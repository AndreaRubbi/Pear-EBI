"""A single place to get a TOML parser across the supported interpreters.

`tomllib` is in the standard library from Python 3.11. On 3.10 it is not, and a bare
`import tomllib` is what made the `test (3.9)` and `test (3.10)` CI jobs fail while
3.11 and 3.12 passed. `tomli` is the same parser under its pre-stdlib name and is
declared in the dev dependency group for `python < 3.11`.
"""

from __future__ import annotations

import sys

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - exercised by the 3.10 CI job
    import tomli as tomllib

__all__ = ["tomllib"]
