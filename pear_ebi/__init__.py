__author__ = "Andrea Rubbi"
__copyright__ = "2023-present Andrea Rubbi and other contributors"
__credits__ = ["Andrea Rubbi", "Lukas Weilguny", "Nick Goldman", "Nicola de Maio"]

__license__ = "MIT"
__version__ = "1.1.2"
__maintainer__ = "Andrea Rubbi"
__institute__ = "EMBL-EBI"
__email__ = "andrear@ebi.ac.uk"
__status__ = "Production"


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
