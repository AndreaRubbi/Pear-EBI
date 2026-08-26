"""Tests for the failure modes found by sweeping the CLI adversarially.

Every case here is a defect that was reproduced before being fixed. They fall into two
kinds, and the second kind is the dangerous one:

  - a raw traceback where a message belonged (ugly, but at least visible);
  - a wrong or incomplete result reported as success with exit code 0.
"""

import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

import numpy as np

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
TWELVE = os.path.join(FIXTURES, "twelve_trees.nwk")
THREE_A = os.path.join(FIXTURES, "three_trees_no_final_newline.nwk")
THREE_B = os.path.join(FIXTURES, "three_trees_b_no_final_newline.nwk")


def _cli():
    found = shutil.which("pear_ebi")
    if found:
        return found
    candidate = os.path.join(os.path.dirname(sys.executable), "pear_ebi")
    return candidate if os.path.exists(candidate) else None


PEAR = _cli()


class ScratchTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="pear_rb_")
        self.cwd = self._tmp.name
        self._old = os.getcwd()
        os.chdir(self.cwd)

    def tearDown(self):
        os.chdir(self._old)
        self._tmp.cleanup()

    def run_cli(self, *args):
        env = dict(os.environ, MPLBACKEND="Agg")
        return subprocess.run(
            [PEAR, *args],
            cwd=self.cwd,
            capture_output=True,
            text=True,
            timeout=1800,
            env=env,
        )

    def assertNoTraceback(self, result, label):
        combined = result.stdout + result.stderr
        self.assertNotIn(
            "Traceback (most recent call last)",
            combined,
            f"{label} dumped a raw traceback at the user:\n{combined[-900:]}",
        )


# ─── Metadata ─────────────────────────────────────────────────────────────────


class TestMetadata(ScratchTestCase):
    """--meta was unusable on single-file input, in three different ways."""

    def _meta(self, name, body):
        with open(os.path.join(self.cwd, name), "w") as fh:
            fh.write(body)
        return name

    def setUp(self):
        super().setUp()
        shutil.copy(TWELVE, os.path.join(self.cwd, "t.nwk"))

    def test_csv_without_a_set_id_column_is_accepted(self):
        """Raised KeyError: 'SET-ID'; the column requirement was undocumented."""
        from pear_ebi.tree_set import tree_set

        self._meta("m.csv", "LK\n" + "\n".join(str(-i) for i in range(12)) + "\n")
        s = tree_set("t.nwk", metadata="m.csv", output_file="d.csv")
        self.assertEqual(s.metadata.shape[0], 12)
        self.assertIn("SET-ID", s.metadata.columns)
        self.assertIn("LK", s.metadata.columns)

    def test_user_supplied_set_id_is_preserved(self):
        from pear_ebi.tree_set import tree_set

        self._meta(
            "m.csv", "SET-ID,LK\n" + "\n".join(f"mine,{-i}" for i in range(12)) + "\n"
        )
        s = tree_set("t.nwk", metadata="m.csv", output_file="d.csv")
        self.assertEqual(s.metadata["SET-ID"].unique().tolist(), ["mine"])

    def test_blank_rows_count_as_rows(self):
        """The class docstring says to leave a blank row for a tree with no info.

        pandas skips blank lines by default, so those rows vanished and every later
        tree's metadata shifted up by one before the count check failed.
        """
        from pear_ebi.tree_set import tree_set

        body = "LK\n" + "\n".join("" if i % 3 else str(-i) for i in range(12)) + "\n"
        self._meta("m.csv", body)
        s = tree_set("t.nwk", metadata="m.csv", output_file="d.csv")
        self.assertEqual(s.metadata.shape[0], 12)

    def test_trailing_spreadsheet_blanks_are_trimmed_but_only_to_the_tree_count(self):
        from pear_ebi.tree_set import tree_set

        body = "LK\n" + "\n".join(str(-i) for i in range(12)) + "\n\n\n\n"
        self._meta("m.csv", body)
        s = tree_set("t.nwk", metadata="m.csv", output_file="d.csv")
        self.assertEqual(s.metadata.shape[0], 12)

    def test_row_count_mismatch_is_refused_with_both_counts(self):
        """Was a raw "Length of values does not match length of index" ValueError."""
        from pear_ebi.tree_set import tree_set

        self._meta("m.csv", "LK\n1\n2\n3\n")
        with self.assertRaises(SystemExit) as ctx:
            tree_set("t.nwk", metadata="m.csv", output_file="d.csv")
        message = str(ctx.exception)
        self.assertIn("12", message)
        self.assertIn("3", message)

    def test_missing_metadata_file_is_named(self):
        from pear_ebi.tree_set import tree_set

        with self.assertRaises(SystemExit) as ctx:
            tree_set("t.nwk", metadata="absent.csv", output_file="d.csv")
        self.assertIn("absent.csv", str(ctx.exception))

    @unittest.skipIf(PEAR is None, "console script not installed")
    def test_meta_through_the_cli_does_not_traceback(self):
        self._meta("m.csv", "LK\n" + "\n".join(str(-i) for i in range(12)) + "\n")
        r = self.run_cli("t.nwk", "--meta", "m.csv", "-m", "hashrf_RF", "-o", "d.csv")
        self.assertNoTraceback(r, "--meta")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)


# ─── Input validation ─────────────────────────────────────────────────────────


@unittest.skipIf(PEAR is None, "console script not installed")
class TestInputValidation(ScratchTestCase):
    """Each of these used to reach the user as a raw traceback."""

    def setUp(self):
        super().setUp()
        shutil.copy(TWELVE, os.path.join(self.cwd, "t.nwk"))
        os.makedirs(os.path.join(self.cwd, "adir"), exist_ok=True)

    def test_missing_input_file(self):
        r = self.run_cli("absent.nwk")
        self.assertNoTraceback(r, "missing input")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not found", r.stdout + r.stderr)

    def test_directory_as_input(self):
        r = self.run_cli("adir")
        self.assertNoTraceback(r, "directory as input")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("directory", r.stdout + r.stderr)

    def test_dir_flag_pointed_at_a_file(self):
        r = self.run_cli("--dir", "t.nwk")
        self.assertNoTraceback(r, "--dir on a file")
        self.assertNotEqual(r.returncode, 0)

    def test_more_dimensions_than_trees(self):
        r = self.run_cli("t.nwk", "-m", "hashrf_RF", "--pcoa", "50")
        self.assertNoTraceback(r, "--pcoa 50 on 12 trees")
        self.assertNotEqual(r.returncode, 0)

    def test_tsne_on_very_few_trees(self):
        """sklearn requires perplexity < n_samples; perplexity was hardcoded at 3."""
        shutil.copy(THREE_A, os.path.join(self.cwd, "three.nwk"))
        r = self.run_cli("three.nwk", "-m", "hashrf_RF", "--tsne", "2")
        self.assertNoTraceback(r, "--tsne on 3 trees")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_pcoa_and_tsne_together_are_refused(self):
        """Only one embedding is carried through, so --pcoa was silently discarded."""
        r = self.run_cli("t.nwk", "-m", "hashrf_RF", "--pcoa", "2", "--tsne", "2")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("cannot be combined", r.stdout + r.stderr)


class TestPerplexityClamp(unittest.TestCase):
    def test_three_is_preserved_whenever_it_is_legal(self):
        """A clamp, not a rescale: existing embeddings must not change."""
        from pear_ebi.embeddings.tSNE_e import _perplexity_for

        for n in (4, 5, 12, 100, 1001):
            with self.subTest(n_samples=n):
                self.assertEqual(_perplexity_for([0] * n), 3.0)

    def test_small_sample_sizes_are_reduced_below_n(self):
        from pear_ebi.embeddings.tSNE_e import _perplexity_for

        for n in (2, 3):
            with self.subTest(n_samples=n):
                self.assertLess(_perplexity_for([0] * n), n)


# ─── Directory expansion ──────────────────────────────────────────────────────


@unittest.skipIf(PEAR is None, "console script not installed")
class TestDirectoryExpansion(ScratchTestCase):
    def setUp(self):
        super().setUp()
        os.makedirs("sets/nested", exist_ok=True)
        shutil.copy(THREE_A, "sets/b.nwk")
        shutil.copy(THREE_B, "sets/a.nwk")
        with open("sets/notes.csv", "w") as fh:
            fh.write("not,a,tree\n")

    def test_non_tree_files_are_skipped_not_counted_as_one_tree(self):
        """--dir defaults to the pattern "*", so a stray notes.csv was swept in."""
        r = self.run_cli("--dir", "sets", "-m", "hashrf_RF", "-o", "c.csv")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("6 trees", r.stdout)
        self.assertIn("Skipping", r.stdout)

    def test_expansion_is_ordered(self):
        """glob() returns filesystem order, which set SET-IDs and matrix row order."""
        r = self.run_cli("--dir", "sets", "-m", "hashrf_RF", "-o", "c.csv")
        self.assertLess(r.stdout.index("a; Containing"), r.stdout.index("b; Containing"))

    def test_pattern_matching_nothing_warns(self):
        r = self.run_cli("--dir", "sets", "--pattern", "*.zzz")
        self.assertIn("no tree files matched", r.stdout + r.stderr)


class TestNonTreeFilesAreNotTrees(unittest.TestCase):
    def test_text_without_parentheses_is_not_a_tree(self):
        from pear_ebi.tree_set import _read_trees

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "notes.csv")
            with open(path, "w") as fh:
                fh.write("alpha,beta,gamma\n1,2,3\n")
            self.assertEqual(_read_trees(path), [])

    def test_a_real_tree_is_still_read(self):
        from pear_ebi.tree_set import _read_trees

        self.assertEqual(len(_read_trees(THREE_A)), 3)


# ─── Output paths ─────────────────────────────────────────────────────────────


@unittest.skipIf(PEAR is None, "console script not installed")
class TestOutputPaths(ScratchTestCase):
    def setUp(self):
        super().setUp()
        shutil.copy(TWELVE, "setA.nwk")
        shutil.copy(THREE_A, "setB.nwk")

    def test_output_flag_is_honoured_verbatim_for_a_collection(self):
        """A UUID was spliced in, so -o wanted.csv produced wanted_<uuid>.csv."""
        r = self.run_cli("setA.nwk", "setB.nwk", "-m", "hashrf_RF", "-o", "wanted.csv")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(os.path.exists("wanted.csv"), os.listdir(self.cwd))

    def test_plots_from_different_inputs_do_not_overwrite_each_other(self):
        """Plots were named after the method only, e.g. PCOA_2D.html."""
        self.run_cli("setA.nwk", "-m", "hashrf_RF", "--pcoa", "2", "-o", "a.csv")
        self.run_cli("setB.nwk", "-m", "hashrf_RF", "--pcoa", "2", "-o", "b.csv")
        html = sorted(f for f in os.listdir(self.cwd) if f.endswith(".html"))
        self.assertGreaterEqual(len(html), 2, f"plots collided: {html}")
        self.assertTrue(any("setA" in f for f in html), html)
        self.assertTrue(any("setB" in f for f in html), html)


class TestSetIdCollisions(ScratchTestCase):
    def test_same_basename_in_two_directories_stays_distinct(self):
        """self.data is keyed by SET-ID, so the second member overwrote the first."""
        from pear_ebi.tree_set import set_collection

        os.makedirs("run1")
        os.makedirs("run2")
        shutil.copy(THREE_A, "run1/trees.nwk")
        shutil.copy(THREE_B, "run2/trees.nwk")
        c = set_collection(
            collection=["run1/trees.nwk", "run2/trees.nwk"], output_file="c.csv"
        )
        self.assertEqual(c.n_trees, 6)
        self.assertEqual(len(c.data), 2, f"SET-IDs collided: {list(c.data)}")

    def test_distinct_basenames_are_not_qualified(self):
        from pear_ebi.tree_set import set_collection

        shutil.copy(THREE_A, "a.nwk")
        shutil.copy(THREE_B, "b.nwk")
        c = set_collection(collection=["a.nwk", "b.nwk"], output_file="c.csv")
        self.assertEqual(sorted(c.data), ["a", "b"])


# ─── Silent corruption ────────────────────────────────────────────────────────


class TestNoSilentlyIncompleteMatrices(ScratchTestCase):
    def test_smart_RF_refuses_trees_with_different_taxa(self):
        """It printed "Done!", exited 0, and wrote a matrix with empty cells."""
        from pear_ebi.calculate_distances import maple_RF

        with open("disjoint.nwk", "w") as fh:
            fh.write("(A:1,B:1,(C:1,D:1):1);\n(E:1,F:1,(G:1,H:1):1);\n")
        with self.assertRaises(SystemExit) as ctx:
            maple_RF.calculate_distance_matrix("disjoint.nwk", 2, "out.csv")
        self.assertIn("undefined", str(ctx.exception))
        self.assertFalse(os.path.exists("out.csv"), "wrote a matrix it had refused")

    def test_smart_RF_still_works_on_a_shared_taxon_set(self):
        from pear_ebi.calculate_distances import maple_RF

        shutil.copy(TWELVE, "t.nwk")
        m = np.asarray(maple_RF.calculate_distance_matrix("t.nwk", 12, "ok.csv"))
        self.assertEqual(m.shape, (12, 12))
        self.assertFalse(np.isnan(m).any())

    def test_tqdist_leaves_no_corrupt_file_when_it_refuses(self):
        from pear_ebi.calculate_distances import tqdist
        from pear_ebi.calculate_distances._exec import PearExecutableError

        shutil.copy(TWELVE, "t.nwk")
        with self.assertRaises(PearExecutableError):
            tqdist.quartet("t.nwk", 99, "fresh.csv")
        self.assertFalse(os.path.exists("fresh.csv"))

    def test_tqdist_preserves_an_existing_matrix_when_it_refuses(self):
        """tqDist writes incrementally, so it used to clobber the path then fail."""
        from pear_ebi.calculate_distances import tqdist
        from pear_ebi.calculate_distances._exec import PearExecutableError

        shutil.copy(TWELVE, "t.nwk")
        good = np.asarray(tqdist.quartet("t.nwk", 12, "keep.csv"))
        with self.assertRaises(PearExecutableError):
            tqdist.quartet("t.nwk", 99, "keep.csv")
        self.assertTrue(os.path.exists("keep.csv"))
        np.testing.assert_array_equal(np.loadtxt("keep.csv", delimiter=","), good)


# ─── Interactive mode ─────────────────────────────────────────────────────────


@unittest.skipIf(PEAR is None, "console script not installed")
class TestInteractiveMode(ScratchTestCase):
    """The menu ran actions via exec() of a source string with no error handling."""

    def setUp(self):
        super().setUp()
        shutil.copy(TWELVE, "t.nwk")
        shutil.copy(TWELVE, "t2.nwk")

    def drive(self, keystrokes, *args):
        env = dict(os.environ, MPLBACKEND="Agg")
        return subprocess.run(
            [PEAR, *args, "-i"],
            input=keystrokes,
            cwd=self.cwd,
            capture_output=True,
            text=True,
            timeout=1800,
            env=env,
        )

    def test_status_and_exit(self):
        r = self.drive("1\n7\n", "t.nwk")
        self.assertNoTraceback(r, "status/exit")
        self.assertEqual(r.returncode, 0)

    def test_eof_at_the_prompt_exits_cleanly(self):
        """Ctrl-D raised an uncaught EOFError traceback."""
        r = self.drive("", "t.nwk")
        self.assertNoTraceback(r, "EOF at prompt")
        self.assertEqual(r.returncode, 0)

    def test_unknown_option_is_reported_and_the_session_continues(self):
        r = self.drive("99\n1\n7\n", "t.nwk")
        self.assertNoTraceback(r, "unknown option")
        self.assertIn("Operation unavailable", r.stdout)
        self.assertEqual(r.returncode, 0)

    def test_non_numeric_input_is_reported(self):
        r = self.drive("abc\n7\n", "t.nwk")
        self.assertNoTraceback(r, "non-numeric input")
        self.assertEqual(r.returncode, 0)

    def test_a_failing_action_does_not_destroy_the_session(self):
        """Asking for 20 dimensions on 12 trees killed the whole session."""
        r = self.drive("3\n1\n20\nn\n1\n7\n", "t.nwk")
        self.assertNoTraceback(r, "failing action")
        self.assertIn("Cannot embed", r.stdout)
        # the tree set is still loaded afterwards
        self.assertGreaterEqual(r.stdout.count("Tree set containing 12 trees"), 2)
        self.assertEqual(r.returncode, 0)

    def test_adding_a_set_rebinds_the_session_state(self):
        """Relied on exec() writing into module globals; PEP 667 ends that."""
        r = self.drive("5\nt2.nwk\n1\n7\n", "t.nwk")
        self.assertNoTraceback(r, "add set")
        self.assertIn("collection containing 24 trees", r.stdout)

    def test_a_bad_path_at_the_file_prompt_reprompts(self):
        r = self.drive("absent.nwk\nt.nwk\n1\n7\n")
        self.assertNoTraceback(r, "bad path at File: prompt")
        self.assertIn("not found", r.stdout)
        self.assertEqual(r.returncode, 0)


if __name__ == "__main__":
    unittest.main()


# ─── Blank metadata cells reaching the plotting code ──────────────────────────


class TestBlankMetadataPlots(ScratchTestCase):
    """Keeping blank rows (correctly) exposed a crash further downstream.

    graph.plot_embedding did `np.unique(metadata[col])`, which sorts, and sorting NaN
    against strings raises
        TypeError: '<' not supported between instances of 'float' and 'str'
    So the very usage the tree_set docstring recommends -- a blank row for a tree with
    no information -- worked at load time and then died as soon as an embedding was
    computed. A blank is now its own labelled category.
    """

    def setUp(self):
        super().setUp()
        shutil.copy(TWELVE, "t.nwk")

    def _meta(self, name, header, values):
        with open(name, "w") as fh:
            fh.write(header + "\n" + "\n".join(values) + "\n")
        return name

    def _plot(self, meta):
        from pear_ebi.tree_set import tree_set

        s = tree_set("t.nwk", metadata=meta, output_file="d.csv")
        s.calculate_distances("hashrf_RF")
        s.embed("pcoa", 2, output="e.csv")
        return s.plot_2D("pcoa", static=True, name_plot="p")

    def test_text_column_with_one_blank(self):
        self._meta("m.csv", "GENE", ["" if i == 4 else f"g{i}" for i in range(12)])
        self.assertIsNotNone(self._plot("m.csv"))

    def test_text_column_with_several_blanks(self):
        self._meta("m.csv", "GENE", ["" if i % 3 else f"g{i}" for i in range(12)])
        self.assertIsNotNone(self._plot("m.csv"))

    def test_entirely_blank_text_column(self):
        self._meta("m.csv", "GENE", [""] * 12)
        self.assertIsNotNone(self._plot("m.csv"))

    def test_numeric_column_with_a_blank(self):
        """A blank in a numeric column would make both colour-scale bounds NaN."""
        self._meta(
            "m.csv", "SUPPORT", ["" if i == 4 else f"{0.1 * i:.2f}" for i in range(12)]
        )
        self.assertIsNotNone(self._plot("m.csv"))

    def test_column_with_no_blanks_is_unaffected(self):
        self._meta("m.csv", "GENE", [f"g{i}" for i in range(12)])
        self.assertIsNotNone(self._plot("m.csv"))

    @unittest.skipIf(PEAR is None, "console script not installed")
    def test_end_to_end_cli_with_blank_text_metadata(self):
        self._meta("m.csv", "GENE", ["" if i == 4 else f"g{i}" for i in range(12)])
        r = self.run_cli("t.nwk", "--meta", "m.csv", "-m", "hashrf_RF", "--pcoa", "2")
        self.assertNoTraceback(r, "blank text metadata + --pcoa")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)


# ─── Temporary files ──────────────────────────────────────────────────────────


class TestNoStrayTempFiles(ScratchTestCase):
    def test_smart_RF_does_not_leave_the_pickle_in_the_working_directory(self):
        """Trees.pckl was written to the cwd and never removed, once per run.

        Only the pypy3 branch writes it, so a shim standing in for pypy3 is what
        exercises the path on a machine without it.
        """
        shutil.copy(TWELVE, "t.nwk")
        os.makedirs("fakebin", exist_ok=True)
        shim = os.path.join(self.cwd, "fakebin", "pypy3")
        with open(shim, "w") as fh:
            fh.write(f'#!/bin/sh\nexec {sys.executable} "$@"\n')
        os.chmod(shim, 0o755)

        env = dict(
            os.environ,
            MPLBACKEND="Agg",
            PATH=os.path.join(self.cwd, "fakebin") + os.pathsep + os.environ["PATH"],
        )
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "from pear_ebi.calculate_distances import maple_RF;"
                "maple_RF.calculate_distance_matrix('t.nwk', 12, 's.csv')",
            ],
            cwd=self.cwd,
            capture_output=True,
            text=True,
            timeout=1800,
            env=env,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(os.path.exists("s.csv"))
        self.assertFalse(
            os.path.exists("Trees.pckl"),
            f"pickle left behind; cwd holds {sorted(os.listdir(self.cwd))}",
        )

    def test_parallel_and_single_core_smart_RF_agree(self):
        """The worker signature changed to carry the pickle directory."""
        from pear_ebi.calculate_distances import maple_RF

        shutil.copy(TWELVE, "t.nwk")
        single = np.asarray(maple_RF.calculate_distance_matrix("t.nwk", 12, "seq.csv"))

        os.makedirs("fakebin", exist_ok=True)
        shim = os.path.join(self.cwd, "fakebin", "pypy3")
        with open(shim, "w") as fh:
            fh.write(f'#!/bin/sh\nexec {sys.executable} "$@"\n')
        os.chmod(shim, 0o755)
        env = dict(
            os.environ,
            MPLBACKEND="Agg",
            PATH=os.path.join(self.cwd, "fakebin") + os.pathsep + os.environ["PATH"],
        )
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "from pear_ebi.calculate_distances import maple_RF;"
                "maple_RF.calculate_distance_matrix('t.nwk', 12, 'par.csv')",
            ],
            cwd=self.cwd,
            capture_output=True,
            text=True,
            timeout=1800,
            env=env,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        parallel = np.loadtxt("par.csv", delimiter=",")
        np.testing.assert_array_equal(parallel, single)


# ─── Telling the user where output went ───────────────────────────────────────


class TestOutputReporting(ScratchTestCase):
    """A run could write three files and name none of them.

    PEAR derives most output names from the input and resolves them against the
    working directory, so `pear_ebi t.nwk -m hashrf_RF --pcoa 2 -p` printed only
    "hashrf_RF | Done!" and "pcoa | Done!" while creating a distance matrix, an
    embedding and an html plot. Reporting them is also what makes an accidental
    overwrite visible.
    """

    def setUp(self):
        super().setUp()
        shutil.copy(TWELVE, "t.nwk")

    def test_paths_inside_the_cwd_are_shown_relative(self):
        from pear_ebi.tree_set import _report_output

        buf = io.StringIO()
        with redirect_stdout(buf):
            _report_output("Thing", os.path.join(self.cwd, "out.csv"))
        text = buf.getvalue()
        self.assertIn("out.csv", text)
        self.assertNotIn(self.cwd, text, "a path in the cwd should be shown relative")

    def test_paths_outside_the_cwd_are_shown_in_full(self):
        from pear_ebi.tree_set import _report_output

        with tempfile.TemporaryDirectory() as other:
            target = os.path.join(other, "out.csv")
            buf = io.StringIO()
            with redirect_stdout(buf):
                _report_output("Thing", target)
            self.assertIn(other, buf.getvalue())

    @unittest.skipIf(PEAR is None, "console script not installed")
    def test_cli_names_every_file_it_writes(self):
        r = self.run_cli("t.nwk", "-m", "hashrf_RF", "--pcoa", "2", "-p")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

        produced = {f for f in os.listdir(self.cwd) if f != "t.nwk"}
        self.assertTrue(produced, "the run wrote nothing")
        for name in produced:
            with self.subTest(file=name):
                self.assertIn(
                    name,
                    r.stdout,
                    f"{name} was written but never mentioned; "
                    f"output was:\n{r.stdout}",
                )

    @unittest.skipIf(PEAR is None, "console script not installed")
    def test_verbose_reports_the_resolved_configuration(self):
        r = self.run_cli("t.nwk", "-m", "hashrf_RF", "-v")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        for expected in (
            "Resolved configuration",
            "working directory",
            "input files",
            "total trees",
            "hashrf binary",
            "platform",
        ):
            with self.subTest(field=expected):
                self.assertIn(expected, r.stdout)
        # the binary actually in use must be named -- one on PATH silently wins
        from pear_ebi import _install_helpers as ih

        self.assertIn(os.path.basename(ih.hashrf_binary()), r.stdout)

    @unittest.skipIf(PEAR is None, "console script not installed")
    def test_verbose_is_off_by_default(self):
        r = self.run_cli("t.nwk", "-m", "hashrf_RF")
        self.assertNotIn("Resolved configuration", r.stdout)

    @unittest.skipIf(PEAR is None, "console script not installed")
    def test_reported_paths_are_not_wrapped_across_lines(self):
        """A path broken over two lines cannot be copied out of the terminal."""
        r = self.run_cli("t.nwk", "-m", "hashrf_RF", "-v")
        for line in r.stdout.splitlines():
            if "->" in line or "directory :" in line:
                with self.subTest(line=line[:40]):
                    self.assertFalse(
                        line.rstrip().endswith(("/", "-")),
                        f"line looks wrapped: {line!r}",
                    )


class TestCapturedStreamSummary(unittest.TestCase):
    """The native tools repeat one diagnostic per tree pair."""

    def test_duplicate_lines_are_collapsed_with_a_count(self):
        from pear_ebi.calculate_distances._exec import format_streams

        # cycling A,B,C x 20 -- not adjacent duplicates, so neighbour-collapsing fails
        stderr = "\n".join(["alpha", "beta", "gamma"] * 20)
        out = format_streams("", stderr)
        self.assertIn("(x20)", out)
        self.assertEqual(out.count("alpha"), 1)
        self.assertLess(len(out.splitlines()), 10)

    def test_long_output_is_capped_and_says_so(self):
        from pear_ebi.calculate_distances._exec import format_streams

        out = format_streams("", "\n".join(f"distinct line {i}" for i in range(200)))
        self.assertIn("more distinct line(s)", out)
        self.assertLess(len(out.splitlines()), 20)

    def test_short_output_is_untouched(self):
        from pear_ebi.calculate_distances._exec import format_streams

        out = format_streams("", "just one problem")
        self.assertIn("just one problem", out)
        self.assertNotIn("(x", out)


class TestTaxonMismatchMessage(ScratchTestCase):
    def test_tqdist_names_the_taxon_mismatch(self):
        """The real cause was buried under 108 lines of repeated diagnostics."""
        from pear_ebi.calculate_distances import tqdist
        from pear_ebi.calculate_distances._exec import PearExecutableError

        with open("mixed.nwk", "w") as fh:
            fh.write(
                "(A:1,B:1,(C:1,D:1):1);\n(A:1,B:1,(C:1,E:1):1);\n"
                "(E:1,F:1,(G:1,H:1):1);\n"
            )
        with self.assertRaises(PearExecutableError) as ctx:
            tqdist.quartet("mixed.nwk", 3, "out.csv")
        message = str(ctx.exception)
        self.assertIn("same set of taxa", message)
        self.assertLess(len(message.splitlines()), 25, "message is still a wall of text")
