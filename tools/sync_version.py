#!/usr/bin/env python
"""Copy the version from pyproject.toml into the other files that carry it.

The version is declared in three places: [project] version in pyproject.toml, which is
what the built distribution carries; __version__ in pear_ebi/__init__.py, which is what
`pear_ebi --version` and the argument parser's banner report; and version: in
CITATION.cff, which is the version people are told to cite. `poetry version` only
updates the first, so a release that forgot the others would ship a package whose
metadata, self-report and citation disagreed.

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
CITATION = REPO / "CITATION.cff"

# `version = "1.2.3"` in the [project] table. Deliberately anchored to the start of a
# line so it cannot match a dependency's version specifier.
PYPROJECT_VERSION = re.compile(r'(?m)^version\s*=\s*"([^"]+)"')
INIT_VERSION = re.compile(r'(?m)^(__version__\s*=\s*")([^"]+)(")')
# `version: 1.2.3` in CITATION.cff, anchored so it cannot match cff-version.
CITATION_VERSION = re.compile(r"(?m)^(version:[ \t]*)(\S+)([ \t]*)$")


def read_pyproject_version() -> str:
    text = PYPROJECT.read_text(encoding="utf-8")
    match = PYPROJECT_VERSION.search(text)
    if not match:
        sys.exit(f'No top-level `version = "..."` found in {PYPROJECT}')
    return match.group(1)


def read_init_version(text: str) -> str:
    match = INIT_VERSION.search(text)
    if not match:
        sys.exit(f'No `__version__ = "..."` found in {INIT}')
    return match.group(2)


def sync(path, pattern, target, *, check):
    """Point one file's version at `target`. Returns (label, changed, ok)."""
    label = path.relative_to(REPO).as_posix()
    # Read as bytes so the file's existing line endings survive the rewrite; several
    # files in this repository are CRLF.
    raw = path.read_bytes()
    newline = b"\r\n" if b"\r\n" in raw else b"\n"
    body = raw.replace(b"\r\n", b"\n").decode("utf-8")

    match = pattern.search(body)
    if not match:
        print(f"{label}: no version line matched {pattern.pattern}", file=sys.stderr)
        return label, False, False

    current = match.group(2)
    if current == target:
        return label, False, True
    if check:
        print(
            f"{label}: has {current}, pyproject.toml has {target}\n"
            "Run: python tools/sync_version.py",
            file=sys.stderr,
        )
        return label, True, False

    updated = pattern.sub(rf"\g<1>{target}\g<3>", body, count=1)
    out = updated.encode("utf-8")
    if newline == b"\r\n":
        out = out.replace(b"\n", b"\r\n")
    path.write_bytes(out)
    print(f"{label}: {current} -> {target}")
    return label, True, True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="only report whether the files agree; exit 1 if they do not",
    )
    args = parser.parse_args()

    target = read_pyproject_version()
    results = [
        sync(INIT, INIT_VERSION, target, check=args.check),
        sync(CITATION, CITATION_VERSION, target, check=args.check),
    ]
    if not all(ok for _, _, ok in results):
        return 1
    if not any(changed for _, changed, _ in results):
        print(
            f"version {target}: pyproject.toml, {results[0][0]} and {results[1][0]} agree"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
