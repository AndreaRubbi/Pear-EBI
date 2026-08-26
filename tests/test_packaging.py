"""Installation and packaging invariants.

These exist because every regression in this area was silent. A bare
find_packages() shipped a top-level `test` package that shadowed the CPython stdlib
module. The wheel grew to 118.8 MB of vendored build byproducts. `numpy < 2` was
declared but nothing checked the resolved version. Removing `tqdm` as "never
imported" broke `import pear_ebi.tree_set`, because pyDRMetrics imports it and
declares no dependencies at all.

The fast tests here run against the *installed* package. The build-based ones are
marked `slow` and skipped unless Poetry is available, since they invoke a real build.
"""

import importlib
import importlib.metadata as md
import os
import pkgutil
import re
import shutil
import subprocess
import sys
import unittest
import zipfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _pyproject():
    from .toml_compat import tomllib

    with open(os.path.join(REPO_ROOT, "pyproject.toml"), "rb") as fh:
        return tomllib.load(fh)


# ─── Metadata consistency ─────────────────────────────────────────────────────


class TestMetadata(unittest.TestCase):
    def test_version_is_consistent_everywhere(self):
        """__version__, pyproject and the installed metadata must agree."""
        import pear_ebi

        declared = _pyproject()["project"]["version"]
        self.assertEqual(pear_ebi.__version__, declared)
        # An editable install's dist-info does not follow a version bump on its own --
        # poetry treats the root project as already current, and pip reports "already
        # satisfied" -- so say what to do rather than just printing two versions.
        self.assertEqual(
            md.version("pear_ebi"),
            declared,
            "the installed metadata is out of step with pyproject.toml; if you just "
            "changed the version, reinstall with `poetry install --only-root` or "
            "`pip install -e . --no-deps --force-reinstall`",
        )

    def test_requires_python_excludes_3_13(self):
        """NumPy 1.x ships no cp313 wheels, so 3.13 cannot work while numpy<2 holds.

        Verified empirically: installing on 3.13 makes pip fall back to building
        scipy from source, which needs a Fortran compiler and fails.
        """
        spec = _pyproject()["project"]["requires-python"]
        self.assertIn(">=3.10", spec)
        self.assertIn("<3.13", spec)

    def test_numpy_pin_is_declared_and_satisfied(self):
        requires = md.requires("pear_ebi") or []
        numpy_reqs = [r for r in requires if r.lower().startswith("numpy")]
        self.assertTrue(numpy_reqs, "numpy is not a declared dependency")
        self.assertTrue(
            any("<2" in r for r in numpy_reqs),
            f"the numpy<2 ceiling is missing from {numpy_reqs}",
        )
        import numpy

        self.assertTrue(
            numpy.__version__.startswith("1."),
            f"numpy {numpy.__version__} is installed but pear_ebi requires < 2.0",
        )

    def test_classifiers_match_requires_python(self):
        """The old classifiers advertised 3.7-3.9, none of which was installable."""
        project = _pyproject()["project"]
        advertised = {
            c.rsplit(" :: ", 1)[-1]
            for c in project["classifiers"]
            if c.startswith("Programming Language :: Python :: 3.")
        }
        self.assertEqual(advertised, {"3.10", "3.11", "3.12"})

    def test_console_script_is_registered_and_resolvable(self):
        entries = [
            ep for ep in md.entry_points(group="console_scripts") if ep.name == "pear_ebi"
        ]
        self.assertTrue(entries, "the pear_ebi console script is not registered")
        self.assertTrue(callable(entries[0].load()))

    def test_notebook_extra_is_declared(self):
        extras = _pyproject()["project"].get("optional-dependencies", {})
        self.assertIn("notebook", extras)
        self.assertTrue(extras["notebook"])

    def test_homepage_points_at_the_real_repository(self):
        """setup.py used to point at the old AndreaRubbi/TreeEmbedding name."""
        urls = _pyproject()["project"]["urls"]
        self.assertIn("Pear-EBI", urls["Homepage"])
        self.assertNotIn("TreeEmbedding", str(urls))


# ─── The installed package ────────────────────────────────────────────────────


class TestInstalledPackage(unittest.TestCase):
    def test_every_submodule_imports(self):
        """Catches a missing transitive dependency anywhere in the package.

        This is what would have caught removing tqdm: pear_ebi never imports it, but
        pyDRMetrics does and declares nothing, so emb_quality -- and therefore
        tree_set -- failed to import.
        """
        import pear_ebi

        skip = {"pear_ebi.test"}  # a stray plotly demo that calls fig.show() on import
        failures = []
        for info in pkgutil.walk_packages(pear_ebi.__path__, prefix="pear_ebi."):
            if info.name in skip or ".linux_bin" in info.name or ".mac_bin" in info.name:
                continue
            try:
                importlib.import_module(info.name)
            except Exception as exc:  # noqa: BLE001 - we want to report any failure
                failures.append(f"{info.name}: {type(exc).__name__}: {exc}")
        self.assertEqual(
            failures, [], "submodules failed to import:\n" + "\n".join(failures)
        )

    def test_no_top_level_test_module_was_installed(self):
        """A bare find_packages() shadowed the CPython stdlib `test` package.

        Only meaningful for a real (non-editable) install. Under `poetry install` the
        project is editable, so pear_ebi.__file__ resolves into the checkout, where a
        sibling tests/ directory is expected. TestBuiltDistribution asserts the same
        invariant against the wheel, which is where it actually matters.
        """
        import pear_ebi

        parent = os.path.dirname(os.path.dirname(pear_ebi.__file__))
        if os.path.exists(os.path.join(parent, "pyproject.toml")):
            self.skipTest("editable install; the wheel test covers this invariant")

        for name in ("test", "tests"):
            candidate = os.path.join(parent, name)
            if os.path.isdir(candidate):
                self.assertFalse(
                    os.path.exists(os.path.join(candidate, "__init__.py")),
                    f"{candidate} was installed and shadows the stdlib `{name}` module",
                )

    def test_native_executables_are_present_and_executable(self):
        from pear_ebi import _install_helpers as ih

        executables = ih.native_executables()
        self.assertEqual(len(executables), 9, "expected hashrf plus 8 tqDist tools")
        for path in executables:
            with self.subTest(binary=os.path.basename(path)):
                self.assertTrue(os.path.exists(path), f"{path} is missing")
                self.assertTrue(os.access(path, os.X_OK), f"{path} is not executable")

    def test_gpl_licence_texts_ship_with_the_binaries(self):
        """HashRF is GPL-2 and tqDist is GPL-3/LGPL-3; the texts must be distributed."""
        from pear_ebi import _install_helpers as ih

        bin_dir = ih.platform_bin_dir()
        required = [
            os.path.join(bin_dir, "HashRF", "COPYING"),
            os.path.join(bin_dir, "tqDist", "COPYING"),
            os.path.join(bin_dir, "tqDist", "COPYING.LESSER"),
        ]
        for path in required:
            with self.subTest(licence=os.path.relpath(path, bin_dir)):
                self.assertTrue(os.path.exists(path), f"{path} is missing")

    def test_gpl_sources_ship_so_the_binaries_can_be_rebuilt(self):
        """Two reasons: the GPL source-provision obligation, and build_tqdist()."""
        from pear_ebi import _install_helpers as ih

        bin_dir = ih.platform_bin_dir()
        self.assertTrue(os.path.exists(os.path.join(bin_dir, "tqDist", "CMakeLists.txt")))
        self.assertTrue(os.path.exists(os.path.join(bin_dir, "HashRF", "configure")))

    def test_no_build_byproducts_ship(self):
        """The wheel used to carry 109 MB of them, 441 files with the maintainer's paths."""
        from pear_ebi import _install_helpers as ih

        bin_dir = os.path.dirname(ih.platform_bin_dir())
        offenders = []
        for dirpath, dirnames, filenames in os.walk(bin_dir):
            if os.path.basename(dirpath) in {"build", "CMakeFiles", ".deps"}:
                offenders.append(os.path.relpath(dirpath, bin_dir))
                dirnames[:] = []
                continue
            for name in filenames:
                if (
                    name.endswith(".o")
                    or name.startswith("._")
                    or name == "CMakeCache.txt"
                ):
                    offenders.append(
                        os.path.relpath(os.path.join(dirpath, name), bin_dir)
                    )
        self.assertEqual(offenders[:10], [], f"{len(offenders)} build byproducts present")

    def test_no_machine_local_paths_in_shipped_text_files(self):
        """239 CMake files used to carry /Users/ar36/... absolute paths."""
        from pear_ebi import _install_helpers as ih

        root = os.path.dirname(ih.platform_bin_dir())
        offenders = []
        for dirpath, _, filenames in os.walk(root):
            for name in filenames:
                path = os.path.join(dirpath, name)
                # mac_bin/HashRF/hashrf embeds linker debug strings and needs a macOS
                # toolchain to strip; it is a documented exception.
                if path.endswith(os.path.join("mac_bin", "HashRF", "hashrf")):
                    continue
                try:
                    with open(path, "rb") as fh:
                        data = fh.read()
                except OSError:
                    continue
                if b"/Users/ar36" in data:
                    offenders.append(os.path.relpath(path, root))
        self.assertEqual(offenders, [], "shipped files contain machine-local paths")


# ─── Build-based invariants ───────────────────────────────────────────────────


@unittest.skipIf(shutil.which("poetry") is None, "poetry not on PATH")
class TestBuiltDistribution(unittest.TestCase):
    """Builds the distribution once and asserts what must be true of the artefacts."""

    wheel = None

    @classmethod
    def setUpClass(cls):
        import tempfile

        cls._tmp = tempfile.TemporaryDirectory(prefix="pear_build_")
        result = subprocess.run(
            ["poetry", "build", "--output", cls._tmp.name],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise unittest.SkipTest(f"poetry build failed:\n{result.stderr[-800:]}")
        wheels = [f for f in os.listdir(cls._tmp.name) if f.endswith(".whl")]
        cls.wheel = os.path.join(cls._tmp.name, wheels[0])

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _names(self):
        with zipfile.ZipFile(self.wheel) as z:
            return z.namelist()

    def test_wheel_has_no_top_level_test_package(self):
        leaked = [n for n in self._names() if n.split("/")[0] in {"test", "tests"}]
        self.assertEqual(leaked, [])

    def test_wheel_preserves_the_executable_bit(self):
        with zipfile.ZipFile(self.wheel) as z:
            for name in [
                "pear_ebi/calculate_distances/linux_bin/HashRF/hashrf",
                "pear_ebi/calculate_distances/linux_bin/tqDist/bin/all_pairs_quartet_dist",
                "pear_ebi/calculate_distances/mac_bin/HashRF/hashrf",
            ]:
                with self.subTest(binary=name):
                    self.assertIn(name, z.namelist())
                    mode = (z.getinfo(name).external_attr >> 16) & 0o777
                    self.assertTrue(mode & 0o111, f"{name} has mode {mode:04o}")

    def test_wheel_stays_small(self):
        """A ratchet. It was 118.8 MB uncompressed, 99.6% vendored build output."""
        with zipfile.ZipFile(self.wheel) as z:
            total_mb = sum(i.file_size for i in z.infolist()) / 1e6
        self.assertLess(total_mb, 20, f"wheel is {total_mb:.1f} MB uncompressed")

    def test_wheel_ships_the_licence_texts(self):
        names = self._names()
        for licence in (
            "pear_ebi/calculate_distances/linux_bin/HashRF/COPYING",
            "pear_ebi/calculate_distances/linux_bin/tqDist/COPYING.LESSER",
        ):
            with self.subTest(licence=licence):
                self.assertIn(licence, names)


class TestHookVersionsMatchTheProject(unittest.TestCase):
    """The formatters used by CI must be the ones the project declares.

    The lint job installs pre-commit standalone, so the `rev:` pinned in
    .pre-commit-config.yaml -- not the version in poetry.lock -- is what decides the
    formatting there. When the two disagree, one reformats what the other just
    formatted and the job can never go green: it sat pinned at black 23.1.0 while the
    project had moved to 26.x. Parsed with a regex rather than PyYAML, which is only in
    the docs dependency group.
    """

    @staticmethod
    def _pinned(tool):
        path = os.path.join(REPO_ROOT, ".pre-commit-config.yaml")
        with open(path, encoding="utf-8") as fh:
            config = fh.read()
        match = re.search(
            r"repo:\s*https://\S*/" + tool + r"\s*\n\s*rev:\s*v?([0-9][^\s#]*)",
            config,
        )
        return match.group(1) if match else None

    def test_black_and_isort_revs_match_the_installed_versions(self):
        for tool in ("black", "isort"):
            with self.subTest(tool=tool):
                pinned = self._pinned(tool)
                self.assertIsNotNone(pinned, f"no {tool} hook found in the config")
                self.assertEqual(
                    pinned,
                    md.version(tool),
                    f"pre-commit pins {tool} {pinned} but the environment has "
                    f"{md.version(tool)}; the two will fight over the formatting",
                )

    def test_blacken_docs_uses_the_same_black(self):
        path = os.path.join(REPO_ROOT, ".pre-commit-config.yaml")
        with open(path, encoding="utf-8") as fh:
            config = fh.read()
        match = re.search(r"additional_dependencies:\s*\[black==([^\]]+)\]", config)
        self.assertIsNotNone(match, "blacken-docs does not pin a black version")
        self.assertEqual(match.group(1), md.version("black"))


if __name__ == "__main__":
    unittest.main()
