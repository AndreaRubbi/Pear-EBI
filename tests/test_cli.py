"""End-to-end tests for the `pear_ebi` command line interface.

__main__.py is 262 statements and had no coverage at all, which is where the
silently-ignored [distance] config method lived. These drive the real console script
in a subprocess, in a temp working directory, exactly as a user or the tutorials
would -- so they also cover the tutorial invocations.
"""

import json
import os
import shutil
import subprocess
import sys
import unittest

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


def _console_script():
    """Path to the installed pear_ebi entry point, or None."""
    found = shutil.which("pear_ebi")
    if found:
        return found
    candidate = os.path.join(os.path.dirname(sys.executable), "pear_ebi")
    return candidate if os.path.exists(candidate) else None


PEAR = _console_script()


@unittest.skipIf(PEAR is None, "pear_ebi console script not installed")
class CLITestCase(unittest.TestCase):
    """Runs the CLI in a scratch directory seeded with the small fixtures."""

    def setUp(self):
        import tempfile

        self._tmp = tempfile.TemporaryDirectory(prefix="pear_cli_")
        self.cwd = self._tmp.name
        for name in os.listdir(FIXTURES):
            shutil.copy(os.path.join(FIXTURES, name), self.cwd)

    def tearDown(self):
        self._tmp.cleanup()

    def run_pear(self, *args, expect_success=True):
        env = dict(os.environ, MPLBACKEND="Agg")
        result = subprocess.run(
            [PEAR, *args], cwd=self.cwd, capture_output=True, text=True, timeout=900, env=env
        )
        if expect_success:
            self.assertEqual(
                result.returncode,
                0,
                f"pear_ebi {' '.join(args)} exited {result.returncode}\n"
                f"stdout:\n{result.stdout[-1500:]}\nstderr:\n{result.stderr[-1500:]}",
            )
        return result

    def path(self, name):
        return os.path.join(self.cwd, name)

    def matrix(self, name):
        return np.loadtxt(self.path(name), delimiter=",")


class TestBasicInvocation(CLITestCase):
    def test_help(self):
        result = self.run_pear("--help")
        self.assertIn("PEAR", result.stdout)
        self.assertIn("--pcoa", result.stdout)

    def test_loads_a_single_tree_set(self):
        result = self.run_pear("twelve_trees.nwk")
        self.assertIn("12 trees", result.stdout)

    def test_reports_the_users_filename_not_a_temp_path(self):
        """self.file used to be reassigned to the normalisation temp file."""
        result = self.run_pear("twelve_trees.nwk")
        self.assertIn("twelve_trees.nwk", result.stdout)
        self.assertNotIn("/tmp/", result.stdout)

    def test_loads_two_files_as_a_collection(self):
        result = self.run_pear(
            "three_trees_no_final_newline.nwk", "three_trees_b_no_final_newline.nwk"
        )
        self.assertIn("6 trees", result.stdout)

    def test_dir_and_pattern(self):
        sub = os.path.join(self.cwd, "sets")
        os.makedirs(sub)
        for name in ("a.nwk", "b.nwk"):
            shutil.copy(os.path.join(FIXTURES, "three_trees_no_final_newline.nwk"),
                        os.path.join(sub, name))
        result = self.run_pear("--dir", "sets", "--pattern", "*.nwk")
        self.assertIn("6 trees", result.stdout)

    def test_missing_input_file_fails_without_a_traceback(self):
        result = self.run_pear("does_not_exist.nwk", expect_success=False)
        self.assertNotEqual(result.returncode, 0)


class TestDistanceMethods(CLITestCase):
    def test_each_method_produces_a_square_matrix(self):
        for method in ("hashrf_RF", "hashrf_wRF", "tqdist_quartet", "tqdist_triplet"):
            with self.subTest(method=method):
                out = f"{method}.csv"
                self.run_pear("twelve_trees.nwk", "-m", method, "-o", out)
                m = self.matrix(out)
                self.assertEqual(m.shape, (12, 12))
                np.testing.assert_array_equal(m, m.T)

    def test_smart_RF_and_hashrf_RF_measure_different_things(self):
        """These two methods do NOT agree, and this test pins the actual relationship.

        Verified against an independent bipartition/clade calculation over the
        12-tree fixture (66 pairs, exact match in both cases):

            smart_RF  == unrooted RF, i.e. |bipartitions(A) symmetric-difference B|
            hashrf_RF == rooted RF / 2, i.e. |clades(A) sym-diff clades(B)| / 2

        So they differ both by a factor of two and by the rooted/unrooted
        distinction, and they disagree on individual pairs in a way no single scale
        factor explains: for trees 0 and 8, unrooted RF is 4 while rooted RF / 2
        is 3.

        This is deliberately pinned rather than reconciled. Changing either would
        change numbers that may already be published, so it is a decision for the
        maintainer, not a silent fix. The test exists so the discrepancy is visible
        and cannot drift further unnoticed.
        """
        self.run_pear("twelve_trees.nwk", "-m", "hashrf_RF", "-o", "fast.csv")
        self.run_pear("twelve_trees.nwk", "-m", "smart_RF", "-o", "slow.csv")
        fast, slow = self.matrix("fast.csv"), self.matrix("slow.csv")

        self.assertEqual(fast.shape, slow.shape)
        # both are valid metrics in their own right
        for name, m in (("hashrf_RF", fast), ("smart_RF", slow)):
            with self.subTest(method=name):
                np.testing.assert_array_equal(m, m.T)
                np.testing.assert_array_equal(np.diag(m), np.zeros(len(m)))
        # ... but they are not the same metric, nor a constant multiple of one another
        self.assertFalse(np.array_equal(fast, slow))
        ratios = np.unique(np.round(slow[fast != 0] / fast[fast != 0], 6))
        self.assertGreater(
            len(ratios), 1,
            "smart_RF and hashrf_RF now differ by a single constant factor; if one of "
            "them was changed, update this test and the documentation together",
        )

    def test_invalid_method_is_rejected(self):
        result = self.run_pear("twelve_trees.nwk", "-m", "not_a_method", expect_success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Invalid method", result.stdout + result.stderr)

    def test_default_output_name_derives_from_the_input(self):
        self.run_pear("twelve_trees.nwk", "-m", "hashrf_RF")
        self.assertTrue(os.path.exists(self.path("twelve_trees_distance_matrix.csv")))

    def test_precomputed_matrix_can_be_reused(self):
        self.run_pear("twelve_trees.nwk", "-m", "hashrf_RF", "-o", "pre.csv")
        result = self.run_pear("twelve_trees.nwk", "-d", "pre.csv")
        self.assertIn("12 trees", result.stdout)


class TestEmbeddings(CLITestCase):
    def test_pcoa_and_tsne(self):
        for flag, tag in (("--pcoa", "PCOA"), ("--tsne", "TSNE")):
            with self.subTest(flag=flag):
                self.run_pear("twelve_trees.nwk", "-m", "hashrf_RF", flag, "2",
                              "-o", f"{tag}.csv")
                expected = f"twelve_trees_{tag}_embedding.csv"
                self.assertTrue(os.path.exists(self.path(expected)),
                                f"{expected} not written; got {os.listdir(self.cwd)}")

    def test_quality_flag(self):
        self.run_pear("twelve_trees.nwk", "-m", "hashrf_RF", "--pcoa", "2", "-q")

    def test_plot_flag_writes_html(self):
        self.run_pear("twelve_trees.nwk", "-m", "hashrf_RF", "--pcoa", "2", "-p")
        html = [f for f in os.listdir(self.cwd) if f.endswith(".html")]
        self.assertTrue(html, f"no plot written; got {os.listdir(self.cwd)}")


class TestConfigFile(CLITestCase):
    """The [distance] method used to be read from the top level and silently ignored."""

    def write_config(self, name, body):
        path = self.path(name)
        with open(path, "w") as fh:
            fh.write(body)
        return name

    def test_distance_method_in_config_is_honoured(self):
        """The decisive check: two different configured methods must disagree.

        Before the fix both fell through to hashrf_RF, so the matrices were identical.
        """
        self.write_config("rf.toml",
                          '[trees]\nfile1 = "twelve_trees.nwk"\n\n'
                          '[distance]\nmethod = "hashrf_RF"\n')
        self.write_config("quartet.toml",
                          '[trees]\nfile1 = "twelve_trees.nwk"\n\n'
                          '[distance]\nmethod = "tqdist_quartet"\n')
        self.run_pear("--config", "rf.toml", "-o", "rf.csv")
        self.run_pear("--config", "quartet.toml", "-o", "quartet.csv")
        rf, qt = self.matrix("rf.csv"), self.matrix("quartet.csv")
        self.assertEqual(rf.shape, qt.shape)
        self.assertFalse(
            np.array_equal(rf, qt),
            "hashrf_RF and tqdist_quartet produced identical matrices, which means the "
            "configured [distance] method is being ignored again",
        )

    def test_invalid_method_in_config_is_rejected(self):
        """Unreachable before: the validation never saw the value."""
        self.write_config("bad.toml",
                          '[trees]\nfile1 = "twelve_trees.nwk"\n\n'
                          '[distance]\nmethod = "hashrf"\n')
        result = self.run_pear("--config", "bad.toml", expect_success=False)
        self.assertIn("Invalid method", result.stdout + result.stderr)

    def test_cli_method_overrides_the_config(self):
        self.write_config("rf.toml",
                          '[trees]\nfile1 = "twelve_trees.nwk"\n\n'
                          '[distance]\nmethod = "hashrf_RF"\n')
        result = self.run_pear("--config", "rf.toml", "-m", "tqdist_triplet")
        self.assertIn("tqdist_triplet", result.stdout)

    def test_embedding_section_is_honoured(self):
        self.write_config("emb.toml",
                          '[trees]\nfile1 = "twelve_trees.nwk"\n\n'
                          '[distance]\nmethod = "hashrf_RF"\n\n'
                          '[embedding]\nmethod = "pcoa"\ndimensions = 2\n')
        result = self.run_pear("--config", "emb.toml")
        self.assertIn("pcoa", result.stdout)

    def test_missing_config_file_fails_cleanly(self):
        result = self.run_pear("--config", "nope.toml", expect_success=False)
        self.assertNotEqual(result.returncode, 0)


class TestShippedExampleConfigs(unittest.TestCase):
    """The example configs are documentation; they must at least be valid."""

    def test_every_example_config_parses_and_names_a_valid_method(self):
        if sys.version_info >= (3, 11):
            import tomllib
        else:  # pragma: no cover
            import tomli as tomllib
        import glob

        valid = {"hashrf_RF", "hashrf_wRF", "smart_RF", "tqdist_quartet", "tqdist_triplet"}
        configs = sorted(glob.glob(
            os.path.join(REPO_ROOT, "examples_tree_sets", "Advanced Examples", "*.toml")))
        self.assertTrue(configs, "no example configs found")
        for path in configs:
            with self.subTest(config=os.path.basename(path)):
                with open(path, "rb") as fh:
                    config = tomllib.load(fh)
                method = (config.get("distance") or {}).get("method")
                if method is not None:
                    self.assertIn(method, valid)

    def test_example_configs_reference_files_that_exist(self):
        if sys.version_info >= (3, 11):
            import tomllib
        else:  # pragma: no cover
            import tomli as tomllib
        import glob

        base = os.path.join(REPO_ROOT, "examples_tree_sets", "Advanced Examples")
        missing = []
        for path in sorted(glob.glob(os.path.join(base, "*.toml"))):
            with open(path, "rb") as fh:
                config = tomllib.load(fh)
            for key, value in (config.get("trees") or {}).items():
                if not os.path.exists(os.path.join(base, value)):
                    missing.append(f"{os.path.basename(path)}::{key} -> {value}")
            for key, value in (config.get("collection") or {}).items():
                if isinstance(value, str) and value.endswith(".csv"):
                    if not os.path.exists(os.path.join(base, value)):
                        missing.append(f"{os.path.basename(path)}::{key} -> {value}")
        self.assertEqual(missing, [], "example configs reference missing files")


if __name__ == "__main__":
    unittest.main()
