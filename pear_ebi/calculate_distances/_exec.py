"""Shared subprocess handling for the bundled native tools (HashRF, tqDist).

hashrf.py and tqdist.py each carried a verbatim copy of ``bash_command`` and
``_run_process``; both live here now.

Two things this module exists to get right:

1. **Both streams are reported.** The previous code captured stdout and then
   deliberately discarded it ("Only include stderr in error messages"). That is
   wrong for these particular tools. HashRF prints
   ``*** Number of trees in the input file: N`` to *stdout*, which is the line that
   actually diagnoses a wrong tree count, and it reports several genuinely fatal
   conditions to stdout while exiting 0 (see hashrf.cc: "file open error",
   "RF rate is only for unweighted RF distance", "ullong add overflow"). Dropping
   stdout turned those into an unexplained "No such file or directory".

2. **Failures raise instead of calling sys.exit().** These are library functions,
   and the documented primary interface is a Jupyter notebook, where ``sys.exit``
   kills the kernel. ``__main__`` catches PearExecutableError and exits non-zero.
"""

import errno
import os
import shutil
import subprocess
import warnings

__all__ = [
    "PearExecutableError",
    "RC_EXEC_FORMAT",
    "RC_NOT_FOUND",
    "RC_PERMISSION",
    "CompletedRun",
    "run_process",
    "resolve_binary",
    "format_streams",
    "remove_file",
]


class PearExecutableError(RuntimeError):
    """A bundled native tool could not be run, or failed while running."""


# Shell-style exit codes, reused as sentinels for exceptions raised before exec.
RC_NOT_FOUND = 127  # binary is not there
RC_PERMISSION = 126  # binary is there but not executable
RC_EXEC_FORMAT = 8  # binary is for a different CPU/OS (ENOEXEC / EBADARCH)

# macOS raises EBADARCH (86) when a Mach-O binary targets another architecture.
# The constant does not exist in Python's errno module on Linux.
_EBADARCH = getattr(errno, "EBADARCH", 86)


class CompletedRun:
    """Result of running a native tool: returncode plus both captured streams."""

    __slots__ = ("returncode", "stdout", "stderr", "argv")

    def __init__(self, returncode, stdout, stderr, argv):
        self.returncode = returncode
        self.stdout = stdout or ""
        self.stderr = stderr or ""
        self.argv = argv

    @property
    def ok(self):
        return self.returncode == 0

    def streams(self):
        """Both streams, labelled, for inclusion in an error message."""
        return format_streams(self.stdout, self.stderr)


def format_streams(stdout, stderr):
    """Render captured output as a labelled block, omitting empty streams.

    Both are included because these tools split their diagnostics across the two
    with no discernible rule.
    """
    parts = []
    out = (stdout or "").strip()
    err = (stderr or "").strip()
    if err:
        parts.append(f"  stderr:\n{_indent(err)}")
    if out:
        parts.append(f"  stdout:\n{_indent(out)}")
    if not parts:
        return "  (the tool produced no output on either stream)"
    return "\n".join(parts)


def _indent(text, prefix="    "):
    return "\n".join(prefix + line for line in text.splitlines())


def run_process(cmd_list, *, timeout=None):
    """Run a native tool, always capturing both streams.

    Never raises for a non-zero exit; inspect the returned CompletedRun. Failures
    that happen before exec are mapped onto the RC_* sentinels so callers can give
    a specific message.
    """
    try:
        completed = subprocess.run(
            cmd_list,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            text=True,
            timeout=timeout,
        )
        return CompletedRun(
            completed.returncode, completed.stdout, completed.stderr, cmd_list
        )
    except FileNotFoundError:
        return CompletedRun(
            RC_NOT_FOUND, "", f"executable not found: {cmd_list[0]}", cmd_list
        )
    except PermissionError as exc:
        return CompletedRun(
            RC_PERMISSION,
            "",
            f"permission denied when trying to execute {cmd_list[0]}: {exc}",
            cmd_list,
        )
    except subprocess.TimeoutExpired as exc:
        return CompletedRun(
            1, exc.stdout or "", f"timed out after {timeout}s", cmd_list
        )
    except OSError as exc:
        code = getattr(exc, "errno", None)
        if code in (errno.ENOEXEC, errno.EFAULT, _EBADARCH):
            return CompletedRun(
                RC_EXEC_FORMAT,
                "",
                f"exec format error when trying to run {cmd_list[0]}: {exc}",
                cmd_list,
            )
        return CompletedRun(1, "", f"os error when running {cmd_list[0]}: {exc}", cmd_list)


def resolve_binary(system_name, packaged_path, *, tool_label):
    """Pick the executable to use, preferring one already on PATH.

    A binary on PATH wins over the bundled copy, which is how a locally built tool
    is honoured -- but it used to happen silently, so a stale or incompatible
    ``hashrf`` on PATH would be used with no indication. That now warns.

    Raises PearExecutableError if neither is available.
    """
    system_bin = shutil.which(system_name)
    if system_bin:
        if packaged_path and os.path.exists(packaged_path):
            warnings.warn(
                f"Using {system_name} found on PATH ({system_bin}) instead of the "
                f"copy bundled with pear_ebi ({packaged_path}). Remove it from PATH "
                f"if you meant to use the bundled one.",
                RuntimeWarning,
                stacklevel=3,
            )
        return system_bin

    if packaged_path and os.path.exists(packaged_path):
        # Last-chance permission repair. ensure_native_executables() normally does
        # this at import; this covers an install whose modes were stripped later.
        if not os.access(packaged_path, os.X_OK):
            try:
                os.chmod(packaged_path, os.stat(packaged_path).st_mode | 0o111)
            except OSError as exc:
                warnings.warn(
                    f"{packaged_path} is not executable and chmod failed ({exc}). "
                    f"Run `chmod +x {packaged_path}`.",
                    RuntimeWarning,
                    stacklevel=3,
                )
        return packaged_path

    from .._install_helpers import describe_platform

    raise PearExecutableError(
        f"{tool_label} executable not found.\n"
        f"  Looked for '{system_name}' on PATH and for the bundled copy at:\n"
        f"    {packaged_path}\n"
        f"  Platform: {describe_platform()}\n"
        f"  pear_ebi ships binaries for Linux x86-64 and macOS arm64 only. On any "
        f"other platform, install {tool_label} system-wide so it is on PATH."
    )


def raise_for_launch_failure(run, bin_path, *, tool_label):
    """Turn a pre-exec failure (RC_* sentinel) into a specific, actionable error."""
    if run.returncode == RC_NOT_FOUND:
        raise PearExecutableError(
            f"{tool_label} executable disappeared at {bin_path}.\n{run.streams()}"
        )

    if run.returncode == RC_PERMISSION:
        raise PearExecutableError(
            f"Permission denied when trying to execute {tool_label} at {bin_path}.\n"
            f"  Fix it with:\n"
            f"    chmod +x {bin_path}\n"
            f"  or from Python:\n"
            f"    import pear_ebi._install_helpers as h; h.ensure_native_executables()\n"
            f"{run.streams()}"
        )

    if run.returncode == RC_EXEC_FORMAT:
        from .._install_helpers import describe_platform

        raise PearExecutableError(
            f"Could not execute {bin_path}: exec format error.\n"
            f"  This means the binary was built for a different CPU or OS.\n"
            f"  Platform: {describe_platform()}\n"
            f"  The bundled Linux binaries are x86-64 and the macOS ones are arm64, "
            f"so an Intel Mac or an ARM Linux box needs a local build:\n"
            f"    import pear_ebi._install_helpers as h; print(h.build_tqdist())\n"
            f"  or install {tool_label} system-wide.\n"
            f"{run.streams()}"
        )


def file_fingerprint(path):
    """Cheap identity of a file: (exists, mtime_ns, size).

    Used to tell "the tool wrote this during the run we just made" from "this is
    left over from an earlier run".
    """
    try:
        st = os.stat(path)
        return (True, st.st_mtime_ns, st.st_size)
    except OSError:
        return (False, 0, 0)


def written_by_this_run(path, before):
    """True if `path` was created or modified since `before` was taken."""
    after = file_fingerprint(path)
    return after[0] and after != before


def remove_file(path):
    """Delete a file, reporting failure instead of hiding it.

    Replaces ``bash_command(f"rm {path}")``, which shelled out with an unquoted
    path (so it broke on spaces), sent both streams to DEVNULL, and returned 0
    unconditionally whether or not the removal worked.
    """
    try:
        os.remove(path)
        return True
    except FileNotFoundError:
        return True
    except OSError as exc:
        warnings.warn(
            f"Could not remove temporary file {path}: {exc}", RuntimeWarning, stacklevel=2
        )
        return False
