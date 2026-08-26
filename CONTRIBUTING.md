# Contributing to PEAR

## Setting up

Dependencies are managed with [Poetry](https://python-poetry.org/). From a checkout:

```
poetry install --with dev,docs
```

`poetry.lock` is the authoritative record of a known-good environment; it is committed and
should be updated in the same commit as any dependency change (`poetry lock` after editing
`pyproject.toml`, then `poetry check --lock` to confirm they agree).

PEAR requires **NumPy < 2.0** and **Python 3.10 – 3.12**. Both bounds are forced from
outside: 3.13 is excluded because NumPy 1.x publishes no cp313 wheels, and 3.9 was dropped
because the patched releases of the dependency tree (pillow, urllib3, jupyter-server and the
rest) are all published as `requires-python >= 3.10`, so supporting 3.9 meant shipping a
lockfile with known advisories.

## Running the tests

```
poetry run pytest -q
```

Two conventions the suite depends on:

- **Tests must not leave anything behind.** Distance matrices, embeddings and plots are
  written to the current working directory, so any test that runs PEAR must `chdir` into a
  temporary directory first — `TempDirTestCase` in `tests/` does this. CI fails the build if
  a test run leaves the working tree dirty, because those files used to end up in commits.
- The fixtures under `tests/fixtures/` include files that deliberately **lack a trailing
  newline**. That is the input that used to make PEAR miscount trees, so it is the whole point
  of those files. `.pre-commit-config.yaml` excludes the directory from `end-of-file-fixer`;
  do not remove that exclusion or the regression tests will start passing for the wrong reason.

## Documentation

```
poetry run mkdocs serve            # local preview
poetry run mkdocs build --strict   # what CI runs; warnings are errors
```

The API pages under `docs/modules/` are generated from docstrings by `gendoc.py`. Note that
regenerating is **destructive**: it removes `docs/tutorials/` and `docs/docs_assets/`, so
restore them with `git checkout -- docs/` afterwards. `mkgendocs.yml` records this.

`site/` is **not** committed. CI builds the documentation on every push to `pear_ebi` and
deploys it to the `gh-pages` branch, so the published site is always whatever the branch
currently builds to. A committed copy would be a second source of truth that could silently
disagree.

## Releasing

Releases are cut by the **Release** workflow (`.github/workflows/release.yml`), run from the
Actions tab. It is the only supported way to publish: it bumps the version, verifies the
artefacts, publishes to PyPI, tags the commit and redeploys the documentation, and it refuses
to do any of that if a check fails.

Two inputs:

- **version** — either a bump keyword (`patch` / `minor` / `major`) or an explicit
  version (`1.1.1`). Passing the version already in `pyproject.toml` is allowed and
  releases it as it stands.
- **dry_run** — run every check and build every artefact, but publish nothing. Use it first.

### Versions are `MAJOR.MINOR.PATCH`

Three segments, which is what Poetry's bump keywords, PEP 440 and everything else assume. So
from `1.1.0`, `patch` gives `1.1.1`, `minor` gives `1.2.0` and `major` gives `2.0.0`. Check
what a keyword would do without changing anything:

```
poetry version --dry-run --short patch
```

Releases up to and including `1.0.1.6` used a fourth segment. Poetry read that as the patch
level and zeroed what followed, so `patch` on `1.0.1.6` produced `1.0.2.0` rather than
`1.0.1.7` — the reason the scheme was changed. Nothing needs doing about the published
history: `1.1.0` sorts above `1.0.1.6` under PEP 440, so the sequence stays monotonic.

`tools/sync_version.py` propagates the version from `pyproject.toml` to
`pear_ebi/__init__.py` and to `CITATION.cff`; `--check` fails if they drift.

One trap worth knowing: if a `pear_ebi.egg-info/` directory is ever left in the repository
root — `python setup.py` used to create one — it takes precedence over the installed
distribution for anything that reads the metadata while the working directory is the
repository, because `""` is on `sys.path`. The symptom is `importlib.metadata.version`
reporting an old version that no file in the tree contains. Delete the directory; the test
suite fails while one is present.

### What the workflow checks before it publishes

In order, stopping at the first failure: the tree is clean and the tag does not already exist;
`pyproject.toml`, `pear_ebi/__init__.py` and `CITATION.cff` agree on the version
(`tools/sync_version.py`, run straight after the bump because `poetry version` edits only the
first of the three); the lockfile
is consistent; the full test suite passes; the wheel and sdist carry the right version, do not
ship a top-level `test` package that would shadow the stdlib, keep the executable bit on the
bundled HashRF and tqDist binaries, keep the GPL licence texts those binaries require, and stay
under 20 MB; the wheel installs into a clean venv **outside the checkout** and computes a real
distance matrix there; and the documentation builds `--strict` with no reintroduced
`polyfill.io` reference.

Only then does it commit, tag `v<version>`, push, publish and deploy.

It runs the full suite and deploys the documentation itself rather than leaving either to CI,
because the release commit is pushed with the workflow's own `GITHUB_TOKEN` and GitHub does
not start new workflow runs from such a push. So no CI run appears for the release commit --
the release workflow's own log is the record that it was checked.

### One-time PyPI setup

Publishing uses [Trusted Publishing](https://docs.pypi.org/trusted-publishers/) — an OIDC
exchange, so there is no API token to store or rotate. **It must be enabled on PyPI once
before the first release, or the publish step will fail.** On
<https://pypi.org/manage/project/pear-ebi/settings/publishing/>, add a GitHub publisher:

| field | value |
| --- | --- |
| Owner | `AndreaRubbi` |
| Repository | `Pear-EBI` |
| Workflow | `release.yml` |
| Environment | *(leave empty)* |

Until that exists, use `dry_run` — everything except the publish and the push is exercised.

## Code style

`black` and `isort` are configured in `pyproject.toml` and run by `pre-commit`:

```
poetry run pre-commit install      # once
poetry run pre-commit run --all-files
```
