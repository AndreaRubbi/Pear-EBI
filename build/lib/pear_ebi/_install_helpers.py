import os
import stat
import warnings
import shutil
import subprocess
import sys


def ensure_native_executables():
    """Ensure packaged native helper files are executable on Unix-like systems.

    This is safe to call at import time: it silently no-ops on non-POSIX
    platforms and ignores missing files.
    """
    # Only attempt on POSIX filesystems where exec bits are relevant
    if os.name != "posix":
        return

    pkg_root = os.path.dirname(__file__)
    calc_dir = os.path.join(pkg_root, "calculate_distances")

    candidates = []

    # HashRF main binary (common location)
    candidates.append(os.path.join(calc_dir, "HashRF", "hashrf"))

    # tqDist binaries under tqDist/bin/
    tq_bin_dir = os.path.join(calc_dir, "tqDist", "bin")
    if os.path.isdir(tq_bin_dir):
        for name in os.listdir(tq_bin_dir):
            candidates.append(os.path.join(tq_bin_dir, name))

    # Also ensure any scripts directly under calculate_distances are executable
    # (be conservative and only add regular files)
    if os.path.isdir(calc_dir):
        for name in ("hashrf", "tqdist", "tqDist"):
            path = os.path.join(calc_dir, name)
            if os.path.isfile(path):
                candidates.append(path)

    made_executable = []
    for path in candidates:
        try:
            if not os.path.exists(path):
                continue
            st = os.stat(path)
            # if any of the exec bits already set, skip
            if st.st_mode & (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH):
                continue
            # set user/group/other execute bits
            new_mode = st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
            os.chmod(path, new_mode)
            made_executable.append(path)
        except Exception:
            # Be tolerant: do not fail import/install if chmod fails.
            # Log a non-fatal warning so packagers/users can notice if needed.
            warnings.warn(f"Could not mark native helper executable: {path}")

    return made_executable


def build_tqdist(package_root=None, use_cmake_first=True, timeout=300):
    """Attempt to (re)build the packaged tqDist native tools in-place.

    This is intended to be run by the user when a shipped binary is not
    executable on the current platform (for example because it was built
    for a different architecture). It will try CMake (if available) and
    fall back to invoking `make` inside the tqDist directory.

    Returns (success: bool, message: str).
    """
    if package_root is None:
        package_root = os.path.dirname(__file__)

    tq_dir = os.path.join(package_root, "calculate_distances", "tqDist")
    if not os.path.isdir(tq_dir):
        return False, f"tqDist source directory not found at {tq_dir}"

    # Prefer CMake build when available
    cmake = shutil.which("cmake")
    make = shutil.which("make")

    try:
        if use_cmake_first and cmake:
            build_dir = os.path.join(tq_dir, "build")
            # Remove stale CMake cache/files in the source tree that could
            # reference an absolute cmake path from another machine.
            cache_in_source = os.path.join(tq_dir, "CMakeCache.txt")
            files_in_source = os.path.join(tq_dir, "CMakeFiles")
            if os.path.exists(cache_in_source):
                try:
                    os.remove(cache_in_source)
                except Exception:
                    # if we can't remove, continue and let cmake surface a clearer error
                    pass
            if os.path.isdir(files_in_source):
                try:
                    # remove the directory tree to force a clean reconfigure
                    import shutil as _sh
                    _sh.rmtree(files_in_source)
                except Exception:
                    pass

            os.makedirs(build_dir, exist_ok=True)
            # run: cmake .. && cmake --build . --config Release
            # Pass a conservative policy version to avoid failures with very
            # new cmake versions that removed old compatibility behavior.
            # Some cmake versions accept CMAKE_POLICY_VERSION or
            # CMAKE_POLICY_VERSION_MINIMUM; try the latter if needed.
            try_args = [cmake, "..", "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"]
            subprocess.check_call(try_args, cwd=build_dir, timeout=timeout)
            # Use cmake --build for portability
            subprocess.check_call([cmake, "--build", ".", "--config", "Release"], cwd=build_dir, timeout=timeout)
            return True, f"Built tqDist with cmake in {build_dir}"

        # Fallback to make in tqDist root
        if make:
            subprocess.check_call([make], cwd=tq_dir, timeout=timeout)
            return True, f"Built tqDist with make in {tq_dir}"

        return False, "Neither cmake nor make found on PATH; install CMake (brew install cmake) or make and retry."
    except subprocess.CalledProcessError as e:
        return False, f"Build command failed with exit {e.returncode}: {e}"
    except FileNotFoundError as e:
        return False, f"Build tool not found: {e}"
    except Exception as e:
        return False, f"Unexpected error while building tqDist: {e}"
