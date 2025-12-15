import os
import subprocess
import sys
import shutil
import errno

import numpy as np
import pandas as pd
import re

# Set the value Display variable
os.environ.setdefault("DISPLAY", ":0.0")

# Simple function to perform bash operations


def bash_command(cmd):
    """Executes bash command in subprocess

    Args:
        cmd (str): bash command to be runned in subprocess

    Returns:
        0: returns 0 if everything's alright
    """
    if os.name == "nt":
        try:
            subprocess.run(
                ["C:\\cygwin64\\bin\\bash.exe", "-l", cmd],
                stderr=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
            )
        except:
            try:
                subprocess.run(
                    ["C:\\cygwin\\bin\\bash.exe", "-l", cmd],
                    stderr=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL,
                )
            except:
                sys.exit("Windows platforms need cygwin64/cygwin to run this subprocess")

    elif os.name == "posix":
        subprocess.run(
            ["/bin/bash", "-c", cmd], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL
        )
    else:
        sys.exit("Could not verify OS")
    return 0


def _bin_dir():
    """Return the directory containing packaged tqDist binaries for this platform.
    Uses linux_bin or mac_bin under calculate_distances.
    """
    try:
        here = os.path.dirname(os.path.abspath(__file__))
    except Exception:
        here = "."

    if sys.platform.startswith("linux"):
        cand = os.path.join(here, "linux_bin", "tqDist")
    elif sys.platform == "darwin":
        cand = os.path.join(here, "mac_bin", "tqDist")
    else:
        cand = None

    if cand is not None and os.path.isdir(cand):
        return cand

    # Fallback to original directory layout
    return os.path.join(here, "tqDist", "bin")


def _run_process(cmd_list, *, capture_stdout=False):
    """Run external command and return (returncode, stdout, stderr).

    Provide helpful diagnostics for common failures (missing binary,
    permission denied, exec-format / arch mismatch).
    """
    try:
        completed = subprocess.run(
            cmd_list,
            stdout=subprocess.PIPE if capture_stdout else subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=False,
            text=True,
        )
        return completed.returncode, completed.stdout if capture_stdout else "", completed.stderr
    except FileNotFoundError:
        return 127, "", f"executable not found: {cmd_list[0]}"
    except PermissionError:
        return 126, "", f"permission denied when trying to execute: {cmd_list[0]}"
    except OSError as e:
        if getattr(e, "errno", None) == errno.ENOEXEC or getattr(e, "errno", None) == errno.EFAULT:
            return 8, "", f"exec format error when trying to run: {cmd_list[0]} ({e})"
        return 1, "", f"os error when running {cmd_list[0]}: {e}"


# ──────────────────────────────── RUNNING HASHRF RETURNING CLEANED OUTPUT ─────
# ? From: "tqDist: a library for computing the quartet and triplet distances between binary or general trees"
# ? by A. Sand, C. N. S. Pedersen et al - 2014
# ? https://www.birc.au.dk/~cstorm/software/tqdist/
# ──────────────────────────────────────────────────────────────────────────────


def quartet(file, n_trees, output_file):
    """Computes quartet distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (pandas.DataFrame): computed distance matrix
    """
    bin_dir = _bin_dir()
    # Support two possible layouts: binaries directly under tqDist/ or under tqDist/bin/
    candidate1 = os.path.join(bin_dir, "all_pairs_quartet_dist")
    candidate2 = os.path.join(bin_dir, "bin", "all_pairs_quartet_dist")
    pkg_bin = candidate1 if os.path.exists(candidate1) else candidate2
    system_bin = shutil.which("all_pairs_quartet_dist")

    if system_bin:
        bin_path = system_bin
    elif os.path.exists(pkg_bin):
        # try to make the packaged binary executable
        if not os.access(pkg_bin, os.X_OK):
            try:
                os.chmod(pkg_bin, os.stat(pkg_bin).st_mode | 0o111)
            except Exception:
                pass
        bin_path = pkg_bin
    else:
        sys.exit(
            f"tqDist quartet binary not found. Ensure native tqDist binaries are present or run `pear_ebi._install_helpers.build_tqdist()` to attempt a local build."
        )

    cmd_list = [bin_path, file, output_file]
    rc, _out, err = _run_process(cmd_list)
    if rc == 127:
        sys.exit(
            f"tqDist quartet executable not found at {bin_path}.\n"
            "You can try rebuilding the packaged native tools by running: `python -c 'import pear_ebi._install_helpers as h; print(h.build_tqdist())'` or install a platform wheel that includes prebuilt binaries."
        )
    if rc == 126:
        sys.exit(f"Permission denied when executing {bin_path}. Try `chmod +x {bin_path}` or call `pear_ebi._install_helpers.ensure_native_executables()`.")
    if rc == 8:
        sys.exit(
            f"Could not execute {bin_path}: exec format error (likely architecture mismatch).\n"
            "Install a matching platform wheel or build tqDist locally (requires cmake and a compiler)."
        )
    if rc != 0:
        msg = err.strip() if err else f"tqDist exited with code {rc}"
        sys.exit(f"tqDist (quartet) failed: {msg}")

    # Reads output file and creates numpy array which is used to create a pandas
    # dataframe. The tqDist executables may separate values with spaces, tabs
    # or commas; be liberal when splitting.
    try:
        with open(output_file, "r") as out:
            distance_matrix = np.zeros((n_trees, n_trees))
            for i, line in enumerate(out):
                tokens = re.split(r"[,\s]+", line.strip())
                tokens = [t for t in tokens if t != ""]
                try:
                    vals = [float(t) for t in tokens]
                except ValueError:
                    sys.exit(
                        f"quartet failed parsing line {i}: '{line.strip()}'.\nParsing error. stderr: {err.strip() if err else ''}"
                    )
                if len(vals) < (i + 1):
                    sys.exit(
                        f"quartet failed: expected {i+1} values on row {i}, got {len(vals)}.\nstderr: {err.strip() if err else ''}"
                    )
                distance_matrix[i, : i + 1] = vals[: i + 1]

        distance_matrix_upper = distance_matrix.transpose()
        distance_matrix = pd.DataFrame(distance_matrix + distance_matrix_upper)
        distance_matrix.to_csv(output_file, header=False, index=False)
        return distance_matrix.values
    except Exception as e:
        stderr_msg = err.strip() if err else ""
        sys.exit(f"quartet failed when reading output: {e}.\n{stderr_msg}")


def triplet(file, n_trees, output_file):
    """Computes triplet distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (pandas.DataFrame): computed distance matrix
    """
    bin_dir = _bin_dir()
    candidate1 = os.path.join(bin_dir, "all_pairs_triplet_dist")
    candidate2 = os.path.join(bin_dir, "bin", "all_pairs_triplet_dist")
    pkg_bin = candidate1 if os.path.exists(candidate1) else candidate2
    system_bin = shutil.which("all_pairs_triplet_dist")

    if system_bin:
        bin_path = system_bin
    elif os.path.exists(pkg_bin):
        if not os.access(pkg_bin, os.X_OK):
            try:
                os.chmod(pkg_bin, os.stat(pkg_bin).st_mode | 0o111)
            except Exception:
                pass
        bin_path = pkg_bin
    else:
        sys.exit(
            f"tqDist triplet binary not found. Ensure native tqDist binaries are present or run `pear_ebi._install_helpers.build_tqdist()` to attempt a local build."
        )

    cmd_list = [bin_path, file, output_file]
    rc, _out, err = _run_process(cmd_list)
    if rc == 127:
        sys.exit(
            f"tqDist triplet executable not found at {bin_path}.\n"
            "Try rebuilding the packaged native tools: `python -c 'import pear_ebi._install_helpers as h; print(h.build_tqdist())'` or install a platform wheel that includes prebuilt binaries."
        )
    if rc == 126:
        sys.exit(f"Permission denied when executing {bin_path}. Try `chmod +x {bin_path}` or call `pear_ebi._install_helpers.ensure_native_executables()`." )
    if rc == 8:
        sys.exit(
            f"Could not execute {bin_path}: exec format error (likely architecture mismatch).\n"
            "Install a matching platform wheel or build tqDist locally (requires cmake and a compiler)."
        )
    if rc != 0:
        msg = err.strip() if err else f"tqDist exited with code {rc}"
        sys.exit(f"tqDist (triplet) failed: {msg}")

    # Reads output file and creates numpy array which is used to create a pandas
    # dataframe. Accept commas, tabs or spaces as separators.
    try:
        with open(output_file, "r") as out:
            distance_matrix = np.zeros((n_trees, n_trees))
            for i, line in enumerate(out):
                tokens = re.split(r"[,\s]+", line.strip())
                tokens = [t for t in tokens if t != ""]
                try:
                    vals = [float(t) for t in tokens]
                except ValueError:
                    sys.exit(
                        f"triplet failed parsing line {i}: '{line.strip()}'.\nParsing error. stderr: {err.strip() if err else ''}"
                    )
                if len(vals) < (i + 1):
                    sys.exit(
                        f"triplet failed: expected {i+1} values on row {i}, got {len(vals)}.\nstderr: {err.strip() if err else ''}"
                    )
                distance_matrix[i, : i + 1] = vals[: i + 1]

        distance_matrix_upper = distance_matrix.transpose()
        distance_matrix = pd.DataFrame(distance_matrix + distance_matrix_upper)
        distance_matrix.to_csv(output_file, header=False, index=False)
        return distance_matrix.values
    except Exception as e:
        stderr_msg = err.strip() if err else ""
        sys.exit(f"triplet failed when reading output: {e}.\n{stderr_msg}")
