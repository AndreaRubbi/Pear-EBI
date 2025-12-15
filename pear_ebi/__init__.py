__author__ = "Andrea Rubbi"
__copyright__ = "2023-present Andrea Rubbi and other contributors"
__credits__ = ["Andrea Rubbi", "Lukas Weilguny", "Nick Goldman", "Nicola de Maio"]

__license__ = "MIT"
__version__ = "1.0.1.6"
__maintainer__ = "Andrea Rubbi"
__institute__ = "EMBL-EBI"
__email__ = "andrear@ebi.ac.uk"
__status__ = "Production"
# Ensure packaged native helpers are executable where relevant. This runs at
# import time and is intentionally tolerant (no exception will be raised).
try:
	# Local import to avoid adding overhead for users who don't need native helpers
	from ._install_helpers import ensure_native_executables

	try:
		ensure_native_executables()
	except Exception:
		# Swallow any unexpected errors to avoid breaking imports for downstream code
		pass
except Exception:
	# If the helper module is missing for any reason, keep going silently.
	pass
