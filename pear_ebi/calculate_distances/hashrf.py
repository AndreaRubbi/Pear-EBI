import os
import subprocess
import sys
import shutil
import errno

import io
import re
import numpy as np
import pandas as pd

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
        subprocess.run(["/bin/bash", "-c", cmd], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    else:
        sys.exit("Could not verify OS")
    return 0


def _bin_dir():
    """Return the directory containing the packaged hashrf binary for this platform.
    Uses linux_bin or mac_bin under calculate_distances.
    """
    try:
        here = os.path.dirname(os.path.abspath(__file__))
    except Exception:
        here = "."

    if sys.platform.startswith("linux"):
        cand = os.path.join(here, "linux_bin", "HashRF")
    elif sys.platform == "darwin":
        cand = os.path.join(here, "mac_bin", "HashRF")
    else:
        cand = None

    if cand is not None and os.path.isdir(cand):
        return cand

    # Fallback to the original HashRF directory inside the package
    return os.path.join(here, "HashRF")


def _parse_hashrf_stdout(text, n_trees, output_file):
    """Attempt to extract the RF matrix from hashrf stdout text.

    Returns numpy array on success, or None on failure.
    """
    lines = text.splitlines()
    # Try to locate the header that precedes the matrix
    start_idx = None
    for i, ln in enumerate(lines):
        if "Robinson-Foulds" in ln or "matrix format" in ln:
            start_idx = i + 1
            break

    # Collect candidate numeric lines
    if start_idx is not None:
        cand = []
        for ln in lines[start_idx:]:
            if re.search(r"\d", ln):
                cand.append(ln)
            elif cand:
                # stop at first non-numeric after we've started
                break
    else:
        cand = [ln for ln in lines if re.search(r"\d", ln)]

    # If we have at least n_trees candidate lines, take the first n_trees
    if len(cand) >= n_trees:
        block = cand[:n_trees]
    else:
        # fallback: try taking the last n_trees numeric lines
        if len(cand) >= n_trees:
            block = cand[-n_trees:]
        else:
            return None

    # Parse numeric tokens per line and build matrix
    try:
        matrix = np.zeros((n_trees, n_trees))
        for i, ln in enumerate(block):
            tokens = re.split(r"[,\s]+", ln.strip())
            tokens = [t for t in tokens if t != ""]
            vals = [float(t) for t in tokens]
            if len(vals) < (i + 1):
                return None
            matrix[i, : i + 1] = vals[: i + 1]
        df = pd.DataFrame(matrix + matrix.transpose())
        df.to_csv(output_file, header=False, index=False, sep=",")
        return df.values
    except Exception:
        return None


def _run_process(cmd_list, *, capture_stdout=False):
    """Run an external command and return (returncode, stdout, stderr).

    This centralises error handling to provide informative messages on common
    failure modes (missing binary, permission errors, exec-format/arch
    mismatches, runtime errors).
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
        # Binary not found at path
        return 127, "", f"executable not found: {cmd_list[0]}"
    except PermissionError:
        return 126, "", f"permission denied when trying to execute: {cmd_list[0]}"
    except OSError as e:
        # Exec format errors and other OS-level errors
        if getattr(e, "errno", None) == errno.ENOEXEC or getattr(e, "errno", None) == errno.EFAULT:
            return 8, "", f"exec format error when trying to run: {cmd_list[0]} ({e})"
        return 1, "", f"os error when running {cmd_list[0]}: {e}"


# ──────────────────────────────── RUNNING HASHRF RETURNING CLEANED OUTPUT ─────
# ? From: "A RANDOMIZED ALGORITHM FOR COMPARING SETS OF PHYLOGENETIC TREES"
# ? by Seung-Jin Sul & Tiffani L. Williams
# ? https://code.google.com/archive/p/hashrf/
# ──────────────────────────────────────────────────────────────────────────────


def hashrf(file, n_trees, output_file):
    """Computes unweighted Robison Foulds distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (pandas.DataFrame): computed distance matrix
    """
    # helper parsing function is defined at module level: _parse_hashrf_stdout

    # Prefer a system-wide `hashrf` if present (user may have built it).
    system_bin = shutil.which("hashrf")

    bin_dir = _bin_dir()
    packaged_bin = os.path.join(bin_dir, "hashrf")

    # If packaged binary exists, make sure it's executable
    if os.path.exists(packaged_bin) and not os.access(packaged_bin, os.X_OK):
        try:
            os.chmod(packaged_bin, os.stat(packaged_bin).st_mode | 0o111)
        except Exception:
            # If chmod fails, we'll fall back to system binary or error later
            pass

    if system_bin:
        bin_path = system_bin
    elif os.path.exists(packaged_bin):
        bin_path = packaged_bin
    else:
        sys.exit(
            "hashrf binary not found. Either install HashRF system-wide or ensure the package includes the HashRF binary for your platform."
        )

    # Build command using the chosen binary (avoid shell where possible)
    cmd_list = [bin_path, file, str(n_trees), "-p", "matrix", "-o", output_file]

    # Run and handle common failure modes with helpful messages. Capture
    # stdout as well because some builds of hashrf print output or errors to
    # STDOUT instead of STDERR. If the binary prints the matrix to stdout we
    # will attempt to parse that as a fallback.
    rc, _out, err = _run_process(cmd_list, capture_stdout=True)

    if rc == 127:
        sys.exit(
            f"hashrf executable not found at {bin_path}.\n"
            "Install HashRF system-wide or ensure the package includes a HashRF binary built for your platform. "
            "If you installed from source, make sure the binary is on your PATH."
        )
    if rc == 126:
        sys.exit(
            f"Permission denied when trying to execute hashRF at {bin_path}.\n"
            "Try making it executable: `chmod +x {bin_path}` or import the package and call ``pear_ebi._install_helpers.ensure_native_executables()`` to fix permissions."
        )
    if rc == 8:
        # Exec format / architecture mismatch
        sys.exit(
            f"Could not execute {bin_path}: exec format error (likely an architecture mismatch).\n"
            "This happens when the packaged binary was built for a different CPU/OS.\n"
            "Solutions: install a platform-specific wheel, build the native tools on your machine, or install HashRF system-wide."
        )
    if rc != 0:
        # If hashrf exited non-zero, check whether the output file exists and
        # contains the matrix (preferred). Try that first.
        if os.path.exists(output_file):
            try:
                # accept comma, semicolon, pipe, space or tab as separators
                df = pd.read_csv(output_file, index_col=None, header=None, sep=r"[,;|\s]+", engine="python")
                if df.shape[1] > n_trees:
                    while df.shape[1] > n_trees:
                        try:
                            df.drop(df.columns[-1], axis=1, inplace=True)
                        except Exception:
                            break
                if df.shape[1] >= n_trees:
                    # coerce to numeric and validate
                    df_num = df.apply(pd.to_numeric, errors="coerce")
                    if df_num.isnull().values.any():
                        loc = list(zip(*np.where(pd.isnull(df_num.values))))
                        r, c = loc[0]
                        raise ValueError(f"non-numeric value in hashrf output at row {r+1}, col {c+1}")
                    df_num.to_csv(output_file, header=False, index=False, sep=",")
                    return df_num.values
            except Exception:
                pass

        # Next try parsing stdout if available (silent fallback)
        if _out:
            parsed = _parse_hashrf_stdout(_out, n_trees=n_trees, output_file=output_file)
            if parsed is not None:
                return parsed

        # Only include stderr in error messages (do not expose stdout)
        msg = (err or "").strip() or f"hashrf exited with code {rc}"
        sys.exit(f"hashrf failed: {msg}")

    # Reads output file and creates numpy array which is used to create a pandas
    # dataframe. If reading fails, include any stderr captured above.
    try:
        with open(output_file, "r") as out:
            # Accept comma, space or tab separated output from hashrf
            distance_matrix = pd.read_csv(out, index_col=None, header=None, sep=r"[,\s]+", engine="python")

        # If the file contains more columns than expected, drop extra trailing
        # columns (some builds append an extra delimiter). If it contains fewer
        # columns than expected, that's an error.
        if distance_matrix.shape[1] > n_trees:
            # drop trailing extra columns until expected width
            while distance_matrix.shape[1] > n_trees:
                try:
                    distance_matrix.drop(distance_matrix.columns[-1], axis=1, inplace=True)
                except Exception:
                    break
        if distance_matrix.shape[1] < n_trees:
            raise ValueError(f"expected at least {n_trees} columns in hashrf output, got {distance_matrix.shape[1]}")

        distance_matrix.to_csv(output_file, header=False, index=False, sep=",")
        return distance_matrix.values
    except Exception as e:
        # If reading the output file failed, but the binary printed the
        # matrix to stdout, try to parse stdout as a fallback.
        if _out:
            parsed = _parse_hashrf_stdout(_out, n_trees=n_trees, output_file=output_file)
            if parsed is not None:
                return parsed
        stderr_msg = (err or "").strip()
        sys.exit(f"hashrf failed when reading output: {e}.\n{stderr_msg}")


# HashRF calculating weighted RF distances
def hashrf_weighted(file, n_trees, output_file):
    """Computes weighted Robison Foulds distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (pandas.DataFrame): computed distance matrix
    """
    # Prefer a system-wide `hashrf` if present.
    system_bin = shutil.which("hashrf")

    bin_dir = _bin_dir()
    packaged_bin = os.path.join(bin_dir, "hashrf")
    if os.path.exists(packaged_bin) and not os.access(packaged_bin, os.X_OK):
        try:
            os.chmod(packaged_bin, os.stat(packaged_bin).st_mode | 0o111)
        except Exception:
            pass

    if system_bin:
        bin_path = system_bin
    elif os.path.exists(packaged_bin):
        bin_path = packaged_bin
    else:
        sys.exit(
            "hashrf binary not found. Either install HashRF system-wide or ensure the package includes the HashRF binary for your platform."
        )

    cmd_list = [bin_path, file, str(n_trees), "-p", "matrix", "-o", output_file, "-w"]

    rc, _out, err = _run_process(cmd_list, capture_stdout=True)
    if rc == 127:
        sys.exit(
            f"hashrf executable not found at {bin_path}. Install HashRF system-wide or ensure the package includes a HashRF binary built for your platform."
        )
    if rc == 126:
        sys.exit(f"Permission denied when trying to execute {bin_path}. Try `chmod +x {bin_path}`.")
    if rc == 8:
        sys.exit(
            f"Could not execute {bin_path}: exec format error (likely an architecture mismatch).\n"
            "Install a matching platform wheel or build HashRF locally."
        )
    if rc != 0:
        # If hashrf exited non-zero, check whether the output file exists and
        # contains the matrix (preferred). Try that first.
        if os.path.exists(output_file):
            try:
                # accept comma, semicolon, pipe, space or tab as separators
                df = pd.read_csv(output_file, index_col=None, header=None, sep=r"[,;|\s]+", engine="python")
                if df.shape[1] > n_trees:
                    while df.shape[1] > n_trees:
                        try:
                            df.drop(df.columns[-1], axis=1, inplace=True)
                        except Exception:
                            break
                if df.shape[1] >= n_trees:
                    # coerce to numeric and validate
                    df_num = df.apply(pd.to_numeric, errors="coerce")
                    if df_num.isnull().values.any():
                        # report location of first non-numeric cell
                        loc = list(zip(*np.where(pd.isnull(df_num.values))))
                        r, c = loc[0]
                        raise ValueError(f"non-numeric value in hashrf output at row {r+1}, col {c+1}")
                    df_num.to_csv(output_file, header=False, index=False, sep=",")
                    return df_num.values
            except Exception:
                pass

        # try parsing stdout as fallback if available
        if _out:
            try:
                df = pd.read_csv(io.StringIO(_out), index_col=None, header=None, sep=r"[,;|\s]+", engine="python")
                if df.shape[1] > 0:
                    try:
                        df.drop(df.columns[-1], axis=1, inplace=True)
                    except Exception:
                        pass
                # coerce/validate numeric
                df_num = df.apply(pd.to_numeric, errors="coerce")
                if df_num.isnull().values.any():
                    loc = list(zip(*np.where(pd.isnull(df_num.values))))
                    r, c = loc[0]
                    raise ValueError(f"non-numeric value in hashrf stdout at row {r+1}, col {c+1}")
                df_num.to_csv(output_file, header=False, index=False, sep=",")
                return df_num.values
            except Exception:
                pass
        # Only include stderr in error messages
        msg = (err or "").strip() or f"hashrf exited with code {rc}"
        sys.exit(f"hashrf failed: {msg}")

    # Post-process output whitespace and read result
    try:
        subprocess.run(
            ["/bin/sh", "-c", f"tr -s ' ' < {output_file} | sed 's/^[ \t]*//' > ./tmp_file && cat ./tmp_file > {output_file} && rm ./tmp_file"],
            check=False,
        )
    except Exception:
        # non-fatal: continue to try reading the file
        pass

    try:
        with open(output_file, "r") as out:
            distance_matrix = pd.read_csv(out, index_col=None, header=None, sep=r"[,;|\s]+", engine="python")

        if distance_matrix.shape[1] > n_trees:
            while distance_matrix.shape[1] > n_trees:
                try:
                    distance_matrix.drop(distance_matrix.columns[-1], axis=1, inplace=True)
                except Exception:
                    break
        if distance_matrix.shape[1] < n_trees:
            raise ValueError(f"expected at least {n_trees} columns in hashrf output, got {distance_matrix.shape[1]}")

        # coerce and validate numeric
        distance_matrix_num = distance_matrix.apply(pd.to_numeric, errors="coerce")
        if distance_matrix_num.isnull().values.any():
            loc = list(zip(*np.where(pd.isnull(distance_matrix_num.values))))
            r, c = loc[0]
            raise ValueError(f"non-numeric value in hashrf output at row {r+1}, col {c+1}")

        distance_matrix_num.to_csv(output_file, header=False, index=False, sep=",")
        return distance_matrix_num.values
    except Exception as e:
        # If reading the output file failed, but the binary printed the
        # matrix to stdout, try to parse stdout as a fallback.
        if _out:
            parsed = _parse_hashrf_stdout(_out, n_trees=n_trees, output_file=output_file)
            if parsed is not None:
                return parsed
        stderr_msg = (err or "").strip()
        sys.exit(f"hashrf failed when reading output: {e}.\n{stderr_msg}")
