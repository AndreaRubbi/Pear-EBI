import os
import re
import warnings

import numpy as np
import pandas as pd

from .._install_helpers import hashrf_binary
from ._exec import (
    PearExecutableError,
    file_fingerprint,
    format_streams,
    raise_for_launch_failure,
    remove_file,
    resolve_binary,
    run_process,
    written_by_this_run,
)

# ──────────────────────────────── RUNNING HASHRF RETURNING CLEANED OUTPUT ─────
# ? From: "A RANDOMIZED ALGORITHM FOR COMPARING SETS OF PHYLOGENETIC TREES"
# ? by Seung-Jin Sul & Tiffani L. Williams
# ? https://code.google.com/archive/p/hashrf/
# ──────────────────────────────────────────────────────────────────────────────


def _parse_hashrf_stdout(text, n_trees, output_file):
    """Last-resort scrape of the RF matrix out of hashrf's stdout.

    Some builds print the matrix to stdout rather than writing the output file.
    Returns a numpy array on success, or None if the text does not look like a
    matrix. Callers warn when this path is taken: a scraped matrix is far less
    trustworthy than one hashrf wrote itself.
    """
    lines = text.splitlines()
    start_idx = None
    for i, ln in enumerate(lines):
        if "Robinson-Foulds" in ln or "matrix format" in ln:
            start_idx = i + 1
            break

    if start_idx is not None:
        cand = []
        for ln in lines[start_idx:]:
            if re.search(r"\d", ln):
                cand.append(ln)
            elif cand:
                break
    else:
        cand = [ln for ln in lines if re.search(r"\d", ln)]

    if len(cand) < n_trees:
        return None
    block = cand[:n_trees]

    try:
        matrix = np.zeros((n_trees, n_trees))
        for i, ln in enumerate(block):
            tokens = [t for t in re.split(r"[,\s]+", ln.strip()) if t != ""]
            vals = [float(t) for t in tokens]
            if len(vals) < (i + 1):
                return None
            matrix[i, : i + 1] = vals[: i + 1]
        df = pd.DataFrame(matrix + matrix.transpose())
        df.to_csv(output_file, header=False, index=False, sep=",")
        return df.values
    except (ValueError, TypeError):
        return None


def _resolve_hashrf():
    """Return the path of the hashrf executable to use."""
    return resolve_binary("hashrf", hashrf_binary(), tool_label="HashRF")


def _read_matrix(output_file, n_trees, label):
    """Read and validate the square distance matrix hashrf wrote.

    Raises ValueError if the file is unusable, so the caller can fold hashrf's own
    output into the message.
    """
    distance_matrix = pd.read_csv(
        output_file, index_col=None, header=None, sep=r"[,;|\s]+", engine="python"
    )

    # Some builds append a trailing delimiter, producing one extra empty column.
    while distance_matrix.shape[1] > n_trees:
        distance_matrix = distance_matrix.iloc[:, :-1]

    if distance_matrix.shape[1] < n_trees:
        raise ValueError(
            f"expected {n_trees} columns in {label} output, got {distance_matrix.shape[1]}. "
            f"This usually means the tree count passed to {label} did not match the "
            f"number of trees actually in the input file."
        )

    numeric = distance_matrix.apply(pd.to_numeric, errors="coerce")
    if numeric.isnull().values.any():
        row, col = (int(v) for v in np.argwhere(pd.isnull(numeric.values))[0])
        raise ValueError(
            f"non-numeric value in {label} output at row {row + 1}, column {col + 1}"
        )

    numeric.to_csv(output_file, header=False, index=False, sep=",")
    return numeric.values


def _finish(run, file, n_trees, output_file, bin_path, label, before):
    """Turn a completed hashrf run into a distance matrix, or raise.

    Shared by the unweighted and weighted entry points. `before` is the output
    file's fingerprint from just before the run.
    """
    raise_for_launch_failure(run, bin_path, tool_label="HashRF")

    if not run.ok:
        # hashrf's exit status is not trustworthy: it writes a complete, correct
        # matrix and still exits 1 on perfectly ordinary input. So a non-zero exit
        # is only fatal if the tool did not actually produce output *in this run*.
        #
        # The `before` fingerprint is what makes this safe. The original code
        # accepted any pre-existing output_file here, which meant a genuinely failed
        # run silently returned a stale matrix from an earlier run -- trimmed to fit
        # n_trees -- and overwrote the user's good CSV with it.
        if written_by_this_run(output_file, before):
            try:
                return _read_matrix(output_file, n_trees, label)
            except (ValueError, pd.errors.ParserError, OSError) as exc:
                raise PearExecutableError(
                    f"{label} exited {run.returncode} and the output it wrote could "
                    f"not be read: {exc}\n"
                    f"  Tree count passed to {label}: {n_trees}\n"
                    f"{run.streams()}"
                ) from exc

        raise PearExecutableError(
            f"{label} failed (exit code {run.returncode}) on {file}.\n"
            f"  Tree count passed to {label}: {n_trees}\n"
            f"{run.streams()}"
        )

    # hashrf reports several fatal conditions on stdout while still exiting 0
    # (see hashrf.cc: "file open error", "RF rate is only for unweighted RF
    # distance", "ullong add overflow"). A zero exit with no output file is
    # therefore a failure, not a read error, and its cause is in stdout.
    if not os.path.exists(output_file):
        parsed = _parse_hashrf_stdout(run.stdout, n_trees, output_file)
        if parsed is not None:
            warnings.warn(
                f"{label} did not write {output_file}; the matrix was recovered by "
                f"parsing its stdout. Treat these distances with caution.",
                RuntimeWarning,
                stacklevel=3,
            )
            return parsed
        raise PearExecutableError(
            f"{label} exited successfully but wrote no output file ({output_file}).\n"
            f"  hashrf reports some fatal errors on stdout while still exiting 0, so "
            f"the cause is most likely below.\n"
            f"  Tree count passed to {label}: {n_trees}\n"
            f"{run.streams()}"
        )

    try:
        return _read_matrix(output_file, n_trees, label)
    except (ValueError, pd.errors.ParserError, OSError) as exc:
        parsed = _parse_hashrf_stdout(run.stdout, n_trees, output_file)
        if parsed is not None:
            warnings.warn(
                f"{label} output file could not be parsed ({exc}); the matrix was "
                f"recovered from its stdout instead. Treat these distances with caution.",
                RuntimeWarning,
                stacklevel=3,
            )
            return parsed
        raise PearExecutableError(
            f"{label} produced output that could not be read: {exc}\n"
            f"  Tree count passed to {label}: {n_trees}\n"
            f"{format_streams(run.stdout, run.stderr)}"
        ) from exc


def hashrf(file, n_trees, output_file):
    """Computes unweighted Robinson Foulds distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (numpy.ndarray): computed distance matrix

    Raises:
        PearExecutableError: if hashrf cannot be run, or fails, or produces
            output that does not match n_trees.
    """
    bin_path = _resolve_hashrf()
    cmd_list = [bin_path, file, str(n_trees), "-p", "matrix", "-o", output_file]
    before = file_fingerprint(output_file)
    run = run_process(cmd_list)
    return _finish(run, file, n_trees, output_file, bin_path, "hashrf", before)


def hashrf_weighted(file, n_trees, output_file):
    """Computes weighted Robinson Foulds distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (numpy.ndarray): computed distance matrix

    Raises:
        PearExecutableError: if hashrf cannot be run, or fails, or produces
            output that does not match n_trees.
    """
    bin_path = _resolve_hashrf()
    cmd_list = [bin_path, file, str(n_trees), "-p", "matrix", "-o", output_file, "-w"]
    before = file_fingerprint(output_file)
    run = run_process(cmd_list)
    # The old code post-processed the output with a `tr -s ' ' | sed` pipeline
    # through a fixed ./tmp_file in the working directory. That was redundant --
    # the pandas reader below already treats runs of whitespace as one separator --
    # and it raced between concurrent runs in the same directory.
    return _finish(run, file, n_trees, output_file, bin_path, "hashrf (weighted)", before)


# Kept as a module-level name because tree_set imports it for cleanup.
__all__ = ["hashrf", "hashrf_weighted", "remove_file", "PearExecutableError"]
