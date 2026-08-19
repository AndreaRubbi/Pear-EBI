import os
import re
import warnings

import numpy as np
import pandas as pd

from .._install_helpers import tqdist_bin_dir
from ._exec import (
    PearExecutableError,
    raise_for_launch_failure,
    resolve_binary,
    run_process,
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
    bin_path = _resolve(system_name, label)
    run = run_process([bin_path, file, output_file])
    raise_for_launch_failure(run, bin_path, tool_label=label)
    if not run.ok:
        raise PearExecutableError(
            f"{label} failed (exit code {run.returncode}) on {file}.\n"
            f"  Tree count passed to {label}: {n_trees}\n"
            f"{run.streams()}"
        )
    return _read_matrix(output_file, n_trees, run, label)


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
