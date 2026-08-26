#!/usr/bin/env python
"""Copy the version from pyproject.toml into pear_ebi/__init__.py.

The version is declared in two places: [project] version in pyproject.toml, which is
what the built distribution carries, and __version__ in pear_ebi/__init__.py, which is
what `pear_ebi --version` and the argument parser's banner report. `poetry version`
only updates the first, so a release that forgot the second would ship a package whose
metadata and self-report disagreed.

The release workflow runs this straight after `poetry version <bump>`, and
tests/test_packaging.py::test_version_is_consistent_everywhere fails if they ever drift.

Usage:
    python tools/sync_version.py          # sync, report what changed
    python tools/sync_version.py --check  # exit 1 if out of step, change nothing
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PYPROJECT = REPO / "pyproject.toml"
INIT = REPO / "pear_ebi" / "__init__.py"

# `version = "1.2.3"` in the [project] table. Deliberately anchored to the start of a
# line so it cannot match a dependency's version specifier.
PYPROJECT_VERSION = re.compile(r'(?m)^version\s*=\s*"([^"]+)"')
INIT_VERSION = re.compile(r'(?m)^(__version__\s*=\s*")([^"]+)(")')


def read_pyproject_version() -> str:
    text = PYPROJECT.read_text(encoding="utf-8")
    match = PYPROJECT_VERSION.search(text)
    if not match:
        sys.exit(f"No top-level `version = \"...\"` found in {PYPROJECT}")
    return match.group(1)


def read_init_version(text: str) -> str:
    match = INIT_VERSION.search(text)
    if not match:
        sys.exit(f"No `__version__ = \"...\"` found in {INIT}")
    return match.group(2)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="only report whether the two agree; exit 1 if they do not",
    )
    args = parser.parse_args()

    target = read_pyproject_version()
    # Read as bytes so the file's existing line endings survive the rewrite; several
    # files in this repository are CRLF.
    raw = INIT.read_bytes()
    newline = b"\r\n" if b"\r\n" in raw else b"\n"
    text = raw.replace(b"\r\n", b"\n").decode("utf-8")
    current = read_init_version(text)

    if current == target:
        print(f"version {target}: pyproject.toml and pear_ebi/__init__.py agree")
        return 0

    if args.check:
        print(
            f"version mismatch: pyproject.toml has {target}, "
            f"pear_ebi/__init__.py has {current}\n"
            "Run: python tools/sync_version.py",
            file=sys.stderr,
        )
        return 1

    updated = INIT_VERSION.sub(rf"\g<1>{target}\g<3>", text, count=1)
    body = updated.encode("utf-8")
    if newline == b"\r\n":
        body = body.replace(b"\n", b"\r\n")
    INIT.write_bytes(body)
    print(f"pear_ebi/__init__.py: {current} -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
