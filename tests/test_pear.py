"""Tests for pear_ebi.

Two things changed structurally here. The file is named test_pear.py, not test.py,
because pytest's default python_files patterns are test_*.py and *_test.py -- so the
old name collected zero items and the suite was only reachable via unittest. And the
fixtures are small files in tests/fixtures rather than the 90 MB sets under
examples_tree_sets, so the suite runs in seconds.

All paths are resolved from __file__. The old suite used a cwd-relative
"../examples_tree_sets/..." and only worked if run from inside its own directory.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import warnings

import numpy as np

from pear_ebi.calculate_distances import tqdist
from pear_ebi.calculate_distances._exec import PearExecutableError, run_process
from pear_ebi.calculate_distances.hashrf import hashrf, hashrf_weighted
from pear_ebi.tree_set import (
    _extract_nexus_trees,
    _read_trees,
    _split_newick,
    set_collection,
    tree_set,
)

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")

FIVE = os.path.join(FIXTURES, "five_trees.nwk")
NO_NL_A = os.path.join(FIXTURES, "three_trees_no_final_newline.nwk")
NO_NL_B = os.path.join(FIXTURES, "three_trees_b_no_final_newline.nwk")
NEXUS = os.path.join(FIXTURES, "nexus_three_trees.trees")
TRICKY = os.path.join(FIXTURES, "two_trees_tricky_semicolons.nwk")
# isomap and lle hardcode n_neighbors=5, so they need more than 5 trees
TWELVE = os.path.join(FIXTURES, "twelve_trees.nwk")

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)


class TempDirTestCase(unittest.TestCase):
    """Runs each test in its own temp directory.

    PEAR writes several outputs relative to the working directory, so without this
    the suite scatters CSVs and HTML through the repo -- and the old suite went
    further and `rm -f`'d its own checked-in golden matrix.
    """

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="pear_test_")
        self.tmp = self._tmp.name
        # chdir as well as providing out(). Several defaults are relative to the
        # working directory -- embed() writes ./<name>_<METHOD>_embedding.csv and
        # plot_* write ./<METHOD>_2D.html -- so a test that does not pass an explicit
        # path scatters files into the repository. Caught exactly that way: running
        # the suite left four twelve_trees_*_embedding.csv files at the repo root.
        self._old_cwd = os.getcwd()
        os.chdir(self.tmp)

    def tearDown(self):
        os.chdir(self._old_cwd)
        self._tmp.cleanup()

    def out(self, name):
        return os.path.join(self.tmp, name)


# ─── Fixture integrity ────────────────────────────────────────────────────────


class TestFixtureIntegrity(unittest.TestCase):
    """The no-trailing-newline fixtures are the regression trigger.

    pre-commit's end-of-file-fixer would silently append a newline and make every
    test below pass for the wrong reason, so .pre-commit-config.yaml excludes
    tests/fixtures. This asserts the exclusion is still working.
    """

    def test_regression_fixtures_have_no_trailing_newline(self):
        for path in (NO_NL_A, NO_NL_B):
            with open(path, "rb") as fh:
                fh.seek(-1, os.SEEK_END)
                last = fh.read(1)
            self.assertEqual(
                last,
                b";",
                f"{os.path.basename(path)} must end in ';' with no trailing newline; "
                "it is the regression trigger. Check the pre-commit exclude.",
            )

    def test_wc_l_undercounts_the_fixture(self):
        """Documents the original failure mode: `wc -l` counts newline bytes."""
        with open(NO_NL_A) as fh:
            newline_count = fh.read().count("\n")
        self.assertEqual(newline_count, 2)
        self.assertEqual(len(_read_trees(NO_NL_A)), 3)


# ─── Bug class (a): tree counting ─────────────────────────────────────────────


class TestNewickParsing(unittest.TestCase):
    def test_counts_trees_without_trailing_newline(self):
        self.assertEqual(len(_read_trees(NO_NL_A)), 3)
        self.assertEqual(len(_read_trees(NO_NL_B)), 3)

    def test_counts_trees_with_trailing_newline(self):
        self.assertEqual(len(_read_trees(FIVE)), 5)

    def test_all_trees_on_one_line(self):
        self.assertEqual(len(_split_newick("(A,B,(C,D));(A,C,(B,D));")), 2)

    def test_tree_split_across_lines(self):
        self.assertEqual(len(_split_newick("(A,\nB,\n(C,D));\n")), 1)

    def test_missing_terminator(self):
        self.assertEqual(len(_split_newick("(A,B,(C,D))")), 1)

    def test_semicolon_in_comment_is_not_a_terminator(self):
        self.assertEqual(len(_read_trees(TRICKY)), 2)

    def test_semicolon_in_quoted_label_is_not_a_terminator(self):
        self.assertEqual(len(_split_newick("('weird;name',B,(C,D));\n")), 1)

    def test_escaped_quote_in_label(self):
        self.assertEqual(len(_split_newick("('O''Brien',B,(C,D));\n")), 1)

    def test_space_inside_quoted_label_is_preserved(self):
        (tree,) = _split_newick("('Homo sapiens',B,(C,D));\n")
        self.assertIn("'Homo sapiens'", tree)

    def test_blank_lines_are_not_trees(self):
        self.assertEqual(len(_split_newick("(A,B,(C,D));\n\n\n(A,C,(B,D));\n")), 2)

    def test_empty_input(self):
        self.assertEqual(_split_newick(""), [])
        self.assertEqual(_split_newick("   \n\n "), [])

    def test_every_tree_is_terminated_and_single_line(self):
        for tree in _read_trees(NO_NL_A):
            self.assertTrue(tree.endswith(";"))
            self.assertNotIn("\n", tree)

    def test_nexus_is_not_counted_by_semicolons(self):
        """A naive count(';') says 6 for this file; the answer is 3."""
        with open(NEXUS) as fh:
            self.assertEqual(fh.read().count(";"), 6)
        self.assertEqual(len(_read_trees(NEXUS)), 3)
        self.assertEqual(len(_extract_nexus_trees(open(NEXUS).read())), 3)


class TestTreeSetCounting(TempDirTestCase):
    def test_tree_set_counts_file_without_trailing_newline(self):
        s = tree_set(NO_NL_A, output_file=self.out("d.csv"))
        self.assertEqual(s.n_trees, 3)

    def test_tree_set_keeps_user_path_and_normalises_separately(self):
        """self.file used to be reassigned to the temp file, leaking into output names."""
        s = tree_set(NO_NL_A, output_file=self.out("d.csv"))
        self.assertEqual(s.file, NO_NL_A)
        self.assertNotEqual(s.tool_input(), s.file)
        self.assertEqual(
            s.metadata["SET-ID"].unique().tolist(), ["three_trees_no_final_newline"]
        )

    def test_normalised_copy_is_one_tree_per_line(self):
        s = tree_set(NO_NL_A, output_file=self.out("d.csv"))
        with open(s.tool_input()) as fh:
            content = fh.read()
        self.assertTrue(content.endswith("\n"))
        self.assertEqual(len([ln for ln in content.splitlines() if ln.strip()]), 3)

    def test_empty_file_is_rejected(self):
        empty = self.out("empty.nwk")
        open(empty, "w").close()
        with self.assertRaises(SystemExit):
            tree_set(empty)


# ─── Bug class (a) continued: the reported concatenation bug ──────────────────


class TestCollectionConcatenation(TempDirTestCase):
    def test_collection_counts_both_files_without_trailing_newlines(self):
        c = set_collection(collection=[NO_NL_A, NO_NL_B], output_file=self.out("c.csv"))
        self.assertEqual(c.n_trees, 6)

    def test_concatenated_input_is_not_glued_together(self):
        """The reported bug: the last tree of A ended up on one line with the first of B."""
        trees = _read_trees(NO_NL_A) + _read_trees(NO_NL_B)
        self.assertEqual(len(trees), 6)
        for tree in trees:
            self.assertEqual(tree.count(";"), 1, f"two trees on one line: {tree}")

    def test_distances_across_all_tools(self):
        """hashrf tolerates a bad layout given the right count; tqDist does not."""
        for method in ("hashrf_RF", "tqdist_quartet", "tqdist_triplet"):
            with self.subTest(method=method):
                c = set_collection(
                    collection=[NO_NL_A, NO_NL_B],
                    output_file=self.out(f"c_{method}.csv"),
                )
                c.calculate_distances(method)
                d = np.asarray(c.distance_matrix)
                self.assertEqual(d.shape, (6, 6))
                np.testing.assert_array_equal(d, d.T)
                self.assertTrue((d != 0).any(), f"{method} returned an all-zeros matrix")

    def test_empty_collection_is_rejected(self):
        c = set_collection(collection=[NO_NL_A], output_file=self.out("c.csv"))
        c.collection = []
        with self.assertRaises(SystemExit):
            c.calculate_distances("hashrf_RF")


class TestNonIdentifierFilenames(TempDirTestCase):
    """Filenames that are not valid Python identifiers used to raise SyntaxError."""

    NAMES = ("my-trees.nwk", "1000trees.nwk", "run 2.nwk", "data.set.v2.nwk")

    def test_set_collection_accepts_awkward_filenames(self):
        paths = []
        for name in self.NAMES:
            dest = self.out(name)
            shutil.copy(NO_NL_A, dest)
            paths.append(dest)
        c = set_collection(collection=paths, output_file=self.out("c.csv"))
        self.assertEqual(c.n_trees, 3 * len(self.NAMES))

    def test_concatenate_accepts_awkward_filenames(self):
        dest = self.out("my-trees.nwk")
        shutil.copy(NO_NL_A, dest)
        c = set_collection(collection=[NO_NL_A], output_file=self.out("c.csv"))
        self.assertEqual(c.concatenate([dest]).n_trees, 6)


# ─── Bug class (b): binary permissions ───────────────────────────────────────


class TestNativeToolPermissions(unittest.TestCase):
    def test_install_helper_finds_the_real_binaries(self):
        """These paths were stale, so the helper always returned an empty list."""
        from pear_ebi import _install_helpers as ih

        self.assertTrue(os.path.isdir(ih.platform_bin_dir()))
        self.assertTrue(os.path.exists(ih.hashrf_binary()))
        self.assertTrue(os.path.isdir(ih.tqdist_bin_dir()))
        executables = ih.native_executables()
        self.assertEqual(len(executables), 9)
        for path in executables:
            self.assertTrue(os.path.exists(path), path)

    def test_ensure_native_executables_restores_the_exec_bit(self):
        from pear_ebi import _install_helpers as ih

        target = ih.hashrf_binary()
        original = os.stat(target).st_mode
        try:
            os.chmod(target, original & ~0o111)
            self.assertFalse(os.access(target, os.X_OK))
            repaired = ih.ensure_native_executables()
            self.assertIn(target, repaired)
            self.assertTrue(os.access(target, os.X_OK))
        finally:
            os.chmod(target, original)

    def test_permission_denied_message_names_the_path(self):
        """The old message printed a literal "{bin_path}" -- it was not an f-string."""
        from pear_ebi import _install_helpers as ih
        from pear_ebi.calculate_distances import _exec

        tmp = tempfile.mkdtemp()
        fake = os.path.join(tmp, "hashrf")
        try:
            shutil.copy2(ih.hashrf_binary(), fake)
            os.chmod(fake, 0o444)
            run = run_process([fake, "x", "2"])
            self.assertEqual(run.returncode, 126)
            with self.assertRaises(PearExecutableError) as ctx:
                _exec.raise_for_launch_failure(run, fake, tool_label="HashRF")
            message = str(ctx.exception)
            self.assertIn(fake, message)
            self.assertNotIn("{bin_path}", message)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


# ─── Bug class (c): NumPy version ────────────────────────────────────────────


class TestNumpyVersionGuard(unittest.TestCase):
    def test_installed_numpy_satisfies_the_pin(self):
        import numpy

        self.assertTrue(
            numpy.__version__.startswith("1."),
            f"pear_ebi requires NumPy < 2.0 but {numpy.__version__} is installed",
        )

    def test_guard_rejects_numpy_2(self):
        """Run in a subprocess so a stub numpy cannot leak into this interpreter."""
        code = (
            "import sys, types\n"
            "stub = types.ModuleType('numpy'); stub.__version__ = '2.3.5'\n"
            "sys.modules['numpy'] = stub\n"
            "try:\n"
            "    import pear_ebi\n"
            "except ImportError as exc:\n"
            "    print('GUARD:' + str(exc).splitlines()[0])\n"
            "else:\n"
            "    print('NOGUARD')\n"
        )
        result = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True
        )
        self.assertIn("GUARD:", result.stdout, result.stdout + result.stderr)
        self.assertIn("NumPy < 2.0", result.stdout)


# ─── Error reporting ─────────────────────────────────────────────────────────


class TestErrorReporting(TempDirTestCase):
    def test_hashrf_failure_includes_both_streams(self):
        """The stdout half of the message used to be discarded by design."""
        single = self.out("one.nwk")
        with open(single, "w") as fh:
            fh.write("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")  # no trailing newline
        with self.assertRaises(PearExecutableError) as ctx:
            hashrf(single, 0, self.out("out.csv"))
        message = str(ctx.exception)
        self.assertIn("at least two trees expected", message)  # stderr
        self.assertIn("Number of trees in the input file: 1", message)  # stdout
        self.assertIn("Tree count passed to hashrf: 0", message)

    def test_failed_run_does_not_return_stale_output(self):
        """A leftover matrix from an earlier run must not be reported as success."""
        out = self.out("shared.csv")
        good = np.asarray(hashrf(FIVE, 5, out))
        self.assertEqual(good.shape, (5, 5))

        single = self.out("one.nwk")
        with open(single, "w") as fh:
            fh.write("(A,B,(C,D));")
        with self.assertRaises(PearExecutableError):
            hashrf(single, 0, out)

        # and the earlier matrix is left intact
        np.testing.assert_array_equal(np.loadtxt(out, delimiter=","), good)

    def test_missing_binary_message_is_actionable(self):
        from pear_ebi.calculate_distances import _exec

        with self.assertRaises(PearExecutableError) as ctx:
            _exec.resolve_binary(
                "definitely-not-a-real-tool",
                "/nonexistent/path/to/tool",
                tool_label="Imaginary",
            )
        message = str(ctx.exception)
        self.assertIn("Imaginary", message)
        self.assertIn("/nonexistent/path/to/tool", message)
        self.assertIn("Platform:", message)


# ─── Core pipeline ───────────────────────────────────────────────────────────


class TestPEAR(TempDirTestCase):
    """The original suite, repaired.

    test_compute_embeddings and test_graph used "pca", which has not been a valid
    method since 5a5a3b8 renamed it to "pcoa" in April 2024; both raised. The
    distance test compared hashrf's return value against the very file hashrf had
    just written, so it asserted X == X and could never fail, and then deleted the
    checked-in 4 MB golden matrix it was supposed to be comparing against.
    """

    def setUp(self):
        super().setUp()
        self.Set = tree_set(FIVE, output_file=self.out("set.csv"))
        self.Collection = set_collection(
            collection=[NO_NL_A, NO_NL_B], output_file=self.out("collection.csv")
        )

    def test_init(self):
        self.assertIsInstance(self.Set, tree_set)
        self.assertIsInstance(self.Collection, set_collection)
        self.assertEqual(self.Set.n_trees, 5)
        self.assertEqual(self.Collection.n_trees, 6)

    def test_calculate_distances_all_methods(self):
        for method in ("hashrf_RF", "hashrf_wRF", "tqdist_quartet", "tqdist_triplet"):
            with self.subTest(method=method):
                s = tree_set(FIVE, output_file=self.out(f"{method}.csv"))
                s.calculate_distances(method)
                d = np.asarray(s.distance_matrix)
                self.assertEqual(d.shape, (5, 5))
                np.testing.assert_array_equal(d, d.T)
                np.testing.assert_array_equal(np.diag(d), np.zeros(5))

    def test_identical_topologies_have_zero_distance(self):
        """Trees 1 and 5 in the fixture are the same topology written differently."""
        s = tree_set(FIVE, output_file=self.out("rf.csv"))
        s.calculate_distances("hashrf_RF")
        self.assertEqual(np.asarray(s.distance_matrix)[0, 4], 0)

    def test_weighted_rf_needs_branch_lengths(self):
        """Documents that an all-zero weighted matrix is correct without lengths."""
        no_len = np.asarray(hashrf_weighted(NO_NL_A, 3, self.out("w0.csv")))
        self.assertFalse((no_len != 0).any())
        with_len = np.asarray(hashrf_weighted(FIVE, 5, self.out("w1.csv")))
        self.assertTrue((with_len != 0).any())

    def test_compute_embeddings(self):
        # Twelve trees rather than five: isomap and lle hardcode n_neighbors=5
        # (Isomap_e.py, LLE_e.py) and sklearn requires n_neighbors < n_samples, so
        # they cannot run on a 5-tree set at all. Exposing that parameter is a
        # separate change; this fixture is simply large enough to exercise them.
        s = tree_set(TWELVE, output_file=self.out("twelve.csv"))
        s.calculate_distances("hashrf_RF")
        for method in ("pcoa", "tsne", "isomap", "lle"):
            with self.subTest(method=method):
                s.embed(method, 2, output=self.out(f"emb_{method}.csv"))
                attr = getattr(s, f"embedding_{method}2D")
                self.assertIsNotNone(attr)
                self.assertEqual(np.asarray(attr).shape[0], 12)

    def test_embed_rejects_shape_mismatch(self):
        """This used to trim the metadata and carry on with a yellow warning."""
        self.Set.calculate_distances("hashrf_RF")
        self.Set.metadata = self.Set.metadata.iloc[:3]
        with self.assertRaises(SystemExit):
            self.Set.embed("pcoa", 2, output=self.out("bad.csv"))

    def test_graph(self):
        self.Set.calculate_distances("hashrf_RF")
        self.Set.embed("pcoa", 2, output=self.out("emb.csv"))
        self.Set.plot_2D("pcoa", static=True, name_plot=self.out("plot"))

    def test_get_subset(self):
        s = tree_set(FIVE, output_file=self.out("sub.csv"))
        s.get_subset(3, method="sequence")

    def test_get_subset_refuses_more_than_available(self):
        s = tree_set(FIVE, output_file=self.out("sub2.csv"))
        with self.assertRaises(SystemExit):
            s.get_subset(99, method="sequence")


# ─── CLI ─────────────────────────────────────────────────────────────────────


class TestParser(unittest.TestCase):
    def test_config_distance_method_is_read(self):
        """config["method"] was read from the top level, so [distance] was ignored."""
        from collections import defaultdict

        from .toml_compat import tomllib

        toml_text = b'[distance]\nmethod = "tqdist_quartet"\n'
        config = defaultdict(lambda: None, tomllib.loads(toml_text.decode()))
        self.assertIsNone(config["method"], "the old lookup silently yielded None")
        config["distance"] = defaultdict(lambda: None, config["distance"])
        self.assertEqual(config["distance"]["method"], "tqdist_quartet")

    def test_shipped_configs_use_valid_methods(self):
        import glob

        from .toml_compat import tomllib

        valid = {
            "hashrf_RF",
            "hashrf_wRF",
            "smart_RF",
            "tqdist_quartet",
            "tqdist_triplet",
        }
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        pattern = os.path.join(root, "examples_tree_sets", "Advanced Examples", "*.toml")
        configs = sorted(glob.glob(pattern))
        self.assertTrue(configs, "no example configs found")
        for path in configs:
            with self.subTest(config=os.path.basename(path)):
                with open(path, "rb") as fh:
                    config = tomllib.load(fh)
                method = (config.get("distance") or {}).get("method")
                if method is not None:
                    self.assertIn(method, valid)


if __name__ == "__main__":
    unittest.main()


# ─── Embedding attribute initialisation ───────────────────────────────────────


class TestEmbeddingAttributes(TempDirTestCase):
    """Only pcoa and tsne used to be initialised, so isomap and lle crashed."""

    METHODS = ("pcoa", "tsne", "isomap", "lle")

    def test_all_eight_attributes_exist_on_a_fresh_tree_set(self):
        s = tree_set(FIVE, output_file=self.out("a.csv"))
        for method in self.METHODS:
            for dim in ("2D", "3D"):
                with self.subTest(attr=f"embedding_{method}{dim}"):
                    self.assertTrue(hasattr(s, f"embedding_{method}{dim}"))

    def test_all_eight_attributes_exist_on_a_fresh_collection(self):
        c = set_collection(collection=[NO_NL_A], output_file=self.out("b.csv"))
        for method in self.METHODS:
            for dim in ("2D", "3D"):
                with self.subTest(attr=f"embedding_{method}{dim}"):
                    self.assertTrue(hasattr(c, f"embedding_{method}{dim}"))

    def test_plot_works_for_every_method_without_a_prior_embed(self):
        """pcoa and tsne computed on demand; isomap and lle raised AttributeError."""
        s = tree_set(TWELVE, output_file=self.out("c.csv"))
        s.calculate_distances("hashrf_RF")
        for method in self.METHODS:
            for label, plot in (("2D", s.plot_2D), ("3D", s.plot_3D)):
                with self.subTest(method=method, dim=label):
                    plot(method, static=True, name_plot=self.out(f"{method}{label}"))


# ─── Notebook executability ───────────────────────────────────────────────────


class TestNotebooksAreExecutable(unittest.TestCase):
    """Cells that start interactive mode must stay tagged skip-execution.

    Under nbconvert the `!pear_ebi -i` subprocess inherits the kernel's open stdin
    pipe, which never delivers data and never reaches EOF, so input() blocks forever:
    the cell HANGS rather than raising, and the notebook cannot be executed
    unattended at all. nbclient honours the skip-execution tag, so this test guards
    the tag rather than the behaviour.
    """

    NOTEBOOKS = [
        "examples_tree_sets/How to use pear_ebi.ipynb",
        "examples_tree_sets/Advanced Examples/Advanced use of pear_ebi.ipynb",
        "examples_tree_sets/Advanced Examples/"
        "How to use pear_ebi on python and reproduce the examples on the paper.ipynb",
    ]

    def _repo(self, rel):
        return os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), rel
        )

    def test_interactive_cells_are_tagged_skip_execution(self):
        import json

        for rel in self.NOTEBOOKS:
            path = self._repo(rel)
            if not os.path.exists(path):
                continue
            with open(path, encoding="utf-8") as fh:
                nb = json.load(fh)
            for i, cell in enumerate(nb["cells"]):
                if cell.get("cell_type") != "code":
                    continue
                src = "".join(cell.get("source", []))
                if "pear_ebi" in src and (" -i" in src or src.rstrip().endswith("-i")):
                    tags = cell.get("metadata", {}).get("tags", [])
                    with self.subTest(notebook=os.path.basename(rel), cell=i):
                        self.assertIn(
                            "skip-execution",
                            tags,
                            f"cell {i} starts interactive mode and would hang "
                            "nbconvert forever; it must be tagged skip-execution",
                        )

    def test_notebooks_are_valid_json_with_a_kernelspec(self):
        import json

        for rel in self.NOTEBOOKS:
            path = self._repo(rel)
            with self.subTest(notebook=os.path.basename(rel)):
                self.assertTrue(os.path.exists(path), f"{rel} is missing")
                with open(path, encoding="utf-8") as fh:
                    nb = json.load(fh)
                self.assertIn("cells", nb)
                self.assertIn("kernelspec", nb.get("metadata", {}))


# ─── hashrf's split canonicalisation ──────────────────────────────────────────


class TestHashrfSplitOrientation(TempDirTestCase):
    """hashrf does not canonicalise a bipartition to a fixed side of the split.

    Two Newick strings can encode the same UNROOTED tree while nesting opposite sides
    of the internal edge. A correct unrooted implementation normalises the split (for
    example by always taking the side without a reference taxon) and reports 0. The
    bundled binary reports 1, and prints "# of unique BIDs = 2" where there is only
    one distinct bipartition.

    Measured impact: 6 of the first 1770 BEAST pairs deviate from true_unrooted_RF/2,
    and 2 of the 66 pairs in twelve_trees.nwk. hashrf_RF is also the method embed()
    silently falls back to.

    Pinned rather than fixed: the behaviour is inside a vendored third-party binary
    with no flag to change it, and correcting it would change numbers that may already
    be published.
    """

    def test_reordering_within_the_string_is_handled(self):
        """Permuting children of the same nesting is fine -- the defect is narrower."""
        path = self.out("perm.nwk")
        with open(path, "w") as fh:
            fh.write("(A:1,B:1,(C:1,D:1):1);\n(B:1,A:1,(D:1,C:1):1);\n")
        m = np.asarray(hashrf(path, 2, self.out("perm.csv")))
        self.assertEqual(m[0, 1], 0)

    def test_opposite_side_of_the_split_is_not_recognised(self):
        """Same unrooted tree, other side of the internal edge nested -> reports 1."""
        path = self.out("flip.nwk")
        with open(path, "w") as fh:
            fh.write("(A:1,B:1,(C:1,D:1):1);\n(C:1,D:1,(A:1,B:1):1);\n")
        m = np.asarray(hashrf(path, 2, self.out("flip.csv")))
        self.assertEqual(
            m[0, 1],
            1,
            "hashrf now reports 0 for the same unrooted tree written with the "
            "opposite side nested. If the vendored binary was fixed or replaced, "
            "update this test and the documentation of hashrf_RF together.",
        )
