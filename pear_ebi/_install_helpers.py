import os
import platform
import shutil
import stat
import subprocess
import sys
import warnings

# ─── Platform binary layout ───────────────────────────────────────────────────
# The bundled native tools live under a per-platform directory:
#
#   calculate_distances/linux_bin/HashRF/hashrf
#   calculate_distances/linux_bin/tqDist/bin/{quartet_dist,triplet_dist,...}
#   calculate_distances/mac_bin/...
#
# They used to live at calculate_distances/{HashRF,tqDist} and were moved in
# 847b151. This module is the single place that knows the layout; hashrf.py and
# tqdist.py resolve their binaries through it rather than each re-deriving paths.
# Keep this module free of heavy imports (no numpy/pandas): it is imported from
# pear_ebi/__init__.py at package import time.

_PLATFORM_DIRS = {"linux": "linux_bin", "darwin": "mac_bin"}

# tqDist ships these executables; everything else in the tree (objects, CMake
# caches, sample data, READMEs) is build residue and must not be marked executable.
TQDIST_EXECUTABLES = (
    "all_pairs_quartet_dist",
    "all_pairs_triplet_dist",
    "pairs_quartet_dist",
    "pairs_triplet_dist",
    "quartet_dist",
    "triplet_dist",
    "test_quartet",
    "test_triplet",
)


def calculate_distances_dir(package_root=None):
    """Return the calculate_distances directory inside the installed package."""
    if package_root is None:
        package_root = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(package_root, "calculate_distances")


def platform_bin_dir(package_root=None):
    """Return the per-platform binary root for this OS, or None if unsupported.

    Returns e.g. .../calculate_distances/linux_bin. The directory is returned even
    when it does not exist, so callers can produce a useful message; use
    ``os.path.isdir`` to test. Returns None on platforms we ship no binaries for
    (notably Windows).
    """
    for prefix, dirname in _PLATFORM_DIRS.items():
        if sys.platform.startswith(prefix):
            return os.path.join(calculate_distances_dir(package_root), dirname)
    return None


def describe_platform():
    """Short human-readable platform string, used in error messages.

    Includes the machine architecture because the bundled macOS binaries are
    arm64-only: an Intel Mac fails with an architecture error, not a missing file.
    """
    return f"{sys.platform}/{platform.machine()}"


def hashrf_binary(package_root=None):
    """Path to the bundled hashrf executable for this platform, or None."""
    bin_dir = platform_bin_dir(package_root)
    if bin_dir is None:
        return None
    return os.path.join(bin_dir, "HashRF", "hashrf")


def tqdist_bin_dir(package_root=None):
    """Directory holding the bundled tqDist executables, or None."""
    bin_dir = platform_bin_dir(package_root)
    if bin_dir is None:
        return None
    return os.path.join(bin_dir, "tqDist", "bin")


def native_executables(package_root=None):
    """Every bundled native executable that should carry the exec bit."""
    paths = []
    hashrf = hashrf_binary(package_root)
    if hashrf:
        paths.append(hashrf)
    tq_dir = tqdist_bin_dir(package_root)
    if tq_dir:
        paths.extend(os.path.join(tq_dir, name) for name in TQDIST_EXECUTABLES)
    return paths


# ─── Permissions ──────────────────────────────────────────────────────────────


def ensure_native_executables(package_root=None):
    """Ensure the bundled native executables are executable on POSIX systems.

    Safe to call at import time: no-ops on non-POSIX platforms and on platforms we
    ship no binaries for, and skips files that are already executable or absent.

    Returns the list of paths whose mode was changed. Raises nothing; problems are
    reported through ``warnings.warn`` so a read-only or shared install is visible
    without breaking the import.
    """
    if os.name != "posix":
        return []

    bin_dir = platform_bin_dir(package_root)
    if bin_dir is None or not os.path.isdir(bin_dir):
        # No bundled binaries for this platform. hashrf.py/tqdist.py will fall back
        # to a system install and produce a targeted message if there is none.
        return []

    made_executable = []
    for path in native_executables(package_root):
        try:
            if not os.path.exists(path):
                continue
            st = os.stat(path)
            if st.st_mode & (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH):
                continue
            os.chmod(path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            made_executable.append(path)
        except OSError as exc:
            warnings.warn(
                f"Could not mark bundled native tool executable: {path} ({exc}). "
                f"If it later fails with 'Permission denied', run `chmod +x {path}`.",
                RuntimeWarning,
                stacklevel=2,
            )

    return made_executable


# ─── Rebuilding tqDist from source ────────────────────────────────────────────


def build_tqdist(package_root=None, use_cmake_first=True, timeout=300):
    """Attempt to (re)build the bundled tqDist native tools in place.

    Intended for the case where a shipped binary will not run on the current
    platform, for example because it was built for a different architecture
    (the bundled macOS binaries are arm64-only).

    Returns (success: bool, message: str).
    """
    bin_dir = platform_bin_dir(package_root)
    if bin_dir is None:
        return False, (
            f"No bundled tqDist sources for this platform ({describe_platform()}). "
            "Install tqDist system-wide so its executables are on PATH."
        )

    tq_dir = os.path.join(bin_dir, "tqDist")
    if not os.path.isdir(tq_dir):
        return False, f"tqDist source directory not found at {tq_dir}"

    cmake = shutil.which("cmake")
    make = shutil.which("make")

    try:
        if use_cmake_first and cmake:
            build_dir = os.path.join(tq_dir, "build")
            # The shipped tree contains a CMake cache generated on the maintainer's
            # machine, holding absolute paths that exist nowhere else. Clear it so
            # cmake reconfigures instead of failing on a stale path.
            cache_in_source = os.path.join(tq_dir, "CMakeCache.txt")
            files_in_source = os.path.join(tq_dir, "CMakeFiles")
            if os.path.exists(cache_in_source):
                try:
                    os.remove(cache_in_source)
                except OSError as exc:
                    warnings.warn(
                        f"Could not remove stale {cache_in_source} ({exc}); "
                        "cmake may fail to reconfigure.",
                        RuntimeWarning,
                        stacklevel=2,
                    )
            if os.path.isdir(files_in_source):
                try:
                    shutil.rmtree(files_in_source)
                except OSError as exc:
                    warnings.warn(
                        f"Could not remove stale {files_in_source} ({exc}); "
                        "cmake may fail to reconfigure.",
                        RuntimeWarning,
                        stacklevel=2,
                    )

            os.makedirs(build_dir, exist_ok=True)
            # A conservative policy version keeps very new cmake releases, which
            # dropped support for the old compatibility behaviour, from erroring.
            subprocess.check_call(
                [cmake, "..", "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"],
                cwd=build_dir,
                timeout=timeout,
            )
            subprocess.check_call(
                [cmake, "--build", ".", "--config", "Release"],
                cwd=build_dir,
                timeout=timeout,
            )
            return True, f"Built tqDist with cmake in {build_dir}"

        if make:
            subprocess.check_call([make], cwd=tq_dir, timeout=timeout)
            return True, f"Built tqDist with make in {tq_dir}"

        return False, (
            "Neither cmake nor make found on PATH. Install one of them "
            "(conda install -c conda-forge cmake, or brew install cmake) and retry."
        )
    except subprocess.CalledProcessError as exc:
        return False, f"Build command failed with exit {exc.returncode}: {exc}"
    except subprocess.TimeoutExpired:
        return False, f"Build timed out after {timeout}s"
    except FileNotFoundError as exc:
        return False, f"Build tool not found: {exc}"
    except OSError as exc:
        return False, f"Unexpected error while building tqDist: {exc}"
