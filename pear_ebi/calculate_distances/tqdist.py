import os
import re
import shutil
import tempfile
import warnings

import numpy as np
import pandas as pd

from .._install_helpers import tqdist_bin_dir
from ._exec import (
    PearExecutableError,
    file_fingerprint,
    raise_for_launch_failure,
    remove_file,
    resolve_binary,
    run_process,
    written_by_this_run,
)

# ──────────────────────────────── RUNNING TQDIST RETURNING CLEANED OUTPUT ─────
# ? From: "tqDist: a library for computing the quartet and triplet distances
# ?        between binary or general trees"
# ? by A. Sand, C. N. S. Pedersen et al - 2014
# ? https://www.birc.au.dk/~cstorm/software/tqdist/
# ──────────────────────────────────────────────────────────────────────────────

# tqDist treats every Newick parse problem as a *warning on stderr* and still
# exits 0 (see NewickParser.cpp: "Parse error! String ended! Continuing anyways...",
# "Parse error! Expected '(' here..."). Its parser is also strictly line-based, so
# a file whose last line has no trailing newline silently loses its last tree.
# Either way the output file ends up short, the untouched rows of the zero-filled
# matrix stay zero, and the result used to be returned as a success. Scanning for
# this marker is what turns that into an error.
_PARSE_ERROR_MARKER = "parse error"


def _resolve(system_name, tool_label):
    bin_dir = tqdist_bin_dir()
    packaged = os.path.join(bin_dir, system_name) if bin_dir else None
    return resolve_binary(system_name, packaged, tool_label=tool_label)


def _read_matrix(output_file, n_trees, run, label):
    """Read tqDist's lower-triangular output into a symmetric matrix.

    Raises PearExecutableError with both of tqDist's streams attached if the file
    is short, unparsable, or if tqDist logged a parse error while exiting 0.
    """
    parse_warned = _PARSE_ERROR_MARKER in run.stderr.lower()

    if not os.path.exists(output_file):
        raise PearExecutableError(
            f"{label} exited successfully but wrote no output file ({output_file}).\n"
            f"  Tree count passed to {label}: {n_trees}\n"
            f"{run.streams()}"
        )

    matrix = np.zeros((n_trees, n_trees))
    rows_seen = 0
    with open(output_file, "r") as out:
        for i, line in enumerate(out):
            if i >= n_trees:
                break
            tokens = [t for t in re.split(r"[,\s]+", line.strip()) if t != ""]
            if not tokens:
                continue
            try:
                vals = [float(t) for t in tokens]
            except ValueError as exc:
                # tqDist writes placeholder text instead of a number when it refuses a
                # pair, and its stderr says why. A taxon-set mismatch is by far the
                # most common cause, so name it up front rather than leaving the user
                # to infer it from a wall of repeated diagnostics.
                if "same set of leaves" in run.stderr or "Leaves doesn't agree" in run.stderr:
                    raise PearExecutableError(
                        f"{label} cannot compare these trees: they do not all have the "
                        f"same set of taxa.\n"
                        f"  Quartet and triplet distances are only defined between trees "
                        f"over the same taxon set.\n"
                        f"  Check that every tree in the input has the same tip labels.\n"
                        f"{run.streams()}"
                    ) from exc
                raise PearExecutableError(
                    f"{label} wrote a non-numeric value on line {i + 1}: "
                    f"{line.strip()!r}\n{run.streams()}"
                ) from exc
            if len(vals) < (i + 1):
                raise PearExecutableError(
                    f"{label} output is malformed: expected at least {i + 1} values on "
                    f"row {i + 1}, got {len(vals)}.\n"
                    f"  Tree count passed to {label}: {n_trees}\n"
                    f"{run.streams()}"
                )
            matrix[i, : i + 1] = vals[: i + 1]
            rows_seen = i + 1

    # A short file means rows of the matrix were never filled and would silently
    # read as zero distances.
    if rows_seen < n_trees:
        raise PearExecutableError(
            f"{label} produced only {rows_seen} of the expected {n_trees} matrix rows, "
            f"so the remaining distances would be zero.\n"
            f"  This usually means the tree count did not match the input file, or "
            f"tqDist could not parse every tree. Note that tqDist reads one tree per "
            f"line and needs a trailing newline after the last tree.\n"
            f"{run.streams()}"
        )

    if parse_warned:
        raise PearExecutableError(
            f"{label} reported a Newick parse error but still exited 0, so the distance "
            f"matrix cannot be trusted.\n"
            f"  tqDist reads one tree per line and requires a newline after the last "
            f"tree.\n{run.streams()}"
        )

    distance_matrix = pd.DataFrame(matrix + matrix.transpose())
    distance_matrix.to_csv(output_file, header=False, index=False)
    return distance_matrix.values


def _run(system_name, label, file, n_trees, output_file):
    """Run a tqDist tool and return the symmetric distance matrix.

    tqDist is pointed at a temporary file rather than at output_file, and the result
    is only written to output_file once it has been validated. Previously tqDist wrote
    straight to output_file and wrote incrementally, so a run that was then refused
    left a truncated matrix behind that looked like a valid result -- and destroyed any
    matrix already at that path.
    """
    bin_path = _resolve(system_name, label)

    with tempfile.TemporaryDirectory(prefix="pear_tqdist_") as tmp:
        staged = os.path.join(tmp, "matrix.csv")
        run = run_process([bin_path, file, staged])
        raise_for_launch_failure(run, bin_path, tool_label=label)
        if not run.ok:
            raise PearExecutableError(
                f"{label} failed (exit code {run.returncode}) on {file}.\n"
                f"  Tree count passed to {label}: {n_trees}\n"
                f"{run.streams()}"
            )
        # _read_matrix validates and rewrites the staged file as a symmetric matrix;
        # it raises before anything reaches output_file if the result is unusable.
        matrix = _read_matrix(staged, n_trees, run, label)
        shutil.copyfile(staged, output_file)

    return matrix


def quartet(file, n_trees, output_file):
    """Computes quartet distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (numpy.ndarray): computed distance matrix

    Raises:
        PearExecutableError: if tqDist cannot be run, fails, or produces a matrix
            that does not match n_trees.
    """
    return _run("all_pairs_quartet_dist", "tqDist (quartet)", file, n_trees, output_file)


def triplet(file, n_trees, output_file):
    """Computes triplet distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (numpy.ndarray): computed distance matrix

    Raises:
        PearExecutableError: if tqDist cannot be run, fails, or produces a matrix
            that does not match n_trees.
    """
    return _run("all_pairs_triplet_dist", "tqDist (triplet)", file, n_trees, output_file)


__all__ = ["quartet", "triplet", "PearExecutableError"]
