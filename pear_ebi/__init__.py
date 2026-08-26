__author__ = "Andrea Rubbi"
__copyright__ = "2023-present Andrea Rubbi and other contributors"
__credits__ = ["Andrea Rubbi", "Lukas Weilguny", "Nick Goldman", "Nicola de Maio"]

__license__ = "MIT"
__version__ = "1.0.1.6"
__maintainer__ = "Andrea Rubbi"
__institute__ = "EMBL-EBI"
__email__ = "andrear@ebi.ac.uk"
__status__ = "Production"


# ─── NumPy version guard ──────────────────────────────────────────────────────
# PEAR requires NumPy < 2.0. The pin is declared in pyproject.toml, but a conda
# environment (or a later `pip install -U numpy`) can still put NumPy 2.x in front
# of us. No first-party code here uses a NumPy-2-removed alias; the incompatibility
# is in the dependency chain -- notably pyDRMetrics 0.0.7, which is unmaintained and
# declares no dependencies at all. Without this check the failure surfaces as an
# obscure traceback from inside a third-party package, so fail early and say why.
def _check_numpy_version() -> None:
    try:
        import numpy
    except ImportError:  # NumPy missing entirely -- let the real import fail later
        return

    version = getattr(numpy, "__version__", "")
    major = version.split(".", 1)[0]
    if major.isdigit() and int(major) >= 2:
        raise ImportError(
            f"pear_ebi requires NumPy < 2.0, but NumPy {version} is installed.\n"
            "NumPy 2.x is not supported: the dependency chain (in particular "
            "pyDRMetrics) has not been ported.\n"
            "Fix it with one of:\n"
            '  pip install "numpy<2"\n'
            '  conda install -c conda-forge "numpy<2"\n'
            "  poetry install        (if you are working from a checkout)"
        )


_check_numpy_version()

# Ensure packaged native helpers are executable where relevant. This runs at
# import time and is intentionally tolerant, but it no longer fails silently:
# a failure here is the usual cause of "Permission denied" from hashrf later on,
# so it is surfaced as a warning.
try:
    # Local import to avoid adding overhead for users who don't need native helpers
    from ._install_helpers import ensure_native_executables

    ensure_native_executables()
except Exception as exc:  # pragma: no cover - depends on install layout
    import warnings

    warnings.warn(
        f"pear_ebi could not verify that its bundled native tools are executable: {exc}. "
        "If hashrf or tqDist later fail with 'Permission denied', run "
        "pear_ebi._install_helpers.ensure_native_executables() or chmod +x the "
        "binaries under pear_ebi/calculate_distances/.",
        RuntimeWarning,
        stacklevel=2,
    )
