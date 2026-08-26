# Contributing to PEAR

## Setting up

Dependencies are managed with [Poetry](https://python-poetry.org/). From a checkout:

```
poetry install --with dev,docs
```

`poetry.lock` is the authoritative record of a known-good environment; it is committed and
should be updated in the same commit as any dependency change (`poetry lock` after editing
`pyproject.toml`, then `poetry check --lock` to confirm they agree).

PEAR supports **Python 3.10 – 3.13** and runs on NumPy 1.x or 2.x. The floor is 3.10
because the patched releases of the dependency tree (pillow, urllib3, jupyter-server and the
rest) are all published as `requires-python >= 3.10`, so supporting 3.9 — end-of-life since
October 2025 — meant shipping a lockfile with known advisories. 3.14 is excluded only because
it has not been tried.

The scientific stack deliberately carries **no upper bounds**. Those ceilings all followed
from a `numpy<2` pin justified by pyDRMetrics 0.0.7 not being ported; that turned out to be
untrue, and tests now assert the ceilings stay off, because re-adding one silently drops
Python 3.13 with it. Note that `anywidget` is a hard dependency, not an extra: from plotly 6
onwards `FigureWidget` — which every interactive PEAR plot is — raises without it.

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

Only then does it publish to PyPI, commit, tag `v<version>`, push and deploy.

Publishing comes before the tag on purpose. PyPI is the one step that can fail for reasons
outside this repository, and if it did so after tagging, the tag `v<version>` would already
be pushed and a retry would stop at "tag already exists" — the version would be burnt. In
this order a failed publish leaves the repository untouched and the same version can be
retried once the cause is fixed.

It runs the full suite and deploys the documentation itself rather than leaving either to CI,
because the release commit is pushed with the workflow's own `GITHUB_TOKEN` and GitHub does
not start new workflow runs from such a push. So no CI run appears for the release commit --
the release workflow's own log is the record that it was checked.

### One-time PyPI setup — do this before the first release

Publishing uses [Trusted Publishing](https://docs.pypi.org/trusted-publishers/): PyPI
verifies an OIDC token that GitHub mints for this workflow, so there is no API token to
create, paste into a secret, or rotate. Nothing about it can be done from the repository —
it has to be granted once on PyPI, by someone with **Owner** or **Maintainer** rights on
the `pear-ebi` project.

1. Sign in to PyPI as the account that owns `pear-ebi`.
2. Go to <https://pypi.org/manage/project/pear-ebi/settings/publishing/>
   (or: your projects → `pear-ebi` → *Manage* → *Publishing*).
3. Under **Add a new publisher**, choose **GitHub** and fill in:

   | field | value |
   | --- | --- |
   | PyPI Project Name | `pear-ebi` |
   | Owner | `AndreaRubbi` |
   | Repository name | `Pear-EBI` |
   | Workflow name | `release.yml` |
   | Environment name | *(leave empty)* |

4. Press **Add**. It takes effect immediately; there is nothing to copy back.

Two-factor authentication is required on PyPI to manage a project, so if the account does
not have it set up yet, that comes first.

That is the whole grant. From then on, "Run workflow" in the Actions tab is the only step
needed to cut a release: the version is bumped, PyPI receives the new distribution, the
commit is tagged and pushed, and the documentation is redeployed.

Until the publisher exists, run with `dry_run` — everything except the publish, the push
and the deploy is exercised, so the artefacts and the test suite are still fully checked.
A release run without it will fail at the publish step with an OIDC permission error, and
because publishing precedes the tag, nothing will have been written to the repository.

## Code style

`black` and `isort` are configured in `pyproject.toml` and run by `pre-commit`:

```
poetry run pre-commit install      # once
poetry run pre-commit run --all-files
```
