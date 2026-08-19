"""Unit-level tests for the modules the end-to-end tests reach only indirectly.

The CLI tests in test_cli.py drive the real console script in a subprocess, which is
the right way to check the user-facing contract but tells coverage nothing and cannot
assert on internal state. These call the same code in-process.
"""

import io
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest import mock

import numpy as np

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
TWELVE = os.path.join(FIXTURES, "twelve_trees.nwk")
NO_NL_A = os.path.join(FIXTURES, "three_trees_no_final_newline.nwk")


# ─── tree_emb_parser ──────────────────────────────────────────────────────────


class TestParser(unittest.TestCase):
    """parser() returned parse_args() with no argv hook, so it was untestable."""

    def parse(self, *argv):
        from pear_ebi.tree_emb_parser import parser

        return parser(list(argv))

    def test_defaults(self):
        args = self.parse()
        self.assertEqual(args.input, [])
        self.assertIsNone(args.method)
        self.assertIsNone(args.pcoa)
        self.assertIsNone(args.tsne)
        self.assertIsNone(args.config)
        self.assertFalse(args.plot)
        self.assertFalse(args.quality)
        self.assertFalse(args.interactive_mode)

    def test_positional_inputs(self):
        self.assertEqual(self.parse("a.nwk", "b.nwk").input, ["a.nwk", "b.nwk"])

    def test_method_and_embeddings(self):
        args = self.parse("a.nwk", "-m", "tqdist_quartet", "--pcoa", "3")
        self.assertEqual(args.method, "tqdist_quartet")
        self.assertEqual(args.pcoa, 3)

    def test_short_and_long_flags_agree(self):
        for short, long, attr in (
            ("-p", "--plot", "plot"),
            ("-q", "--quality", "quality"),
            ("-i", "--interactive", "interactive_mode"),
        ):
            with self.subTest(flag=long):
                self.assertTrue(getattr(self.parse("a.nwk", short), attr))
                self.assertTrue(getattr(self.parse("a.nwk", long), attr))

    def test_dir_and_pattern_are_long_only(self):
        """The docs advertised -dir and -pattern, which argparse mangles."""
        args = self.parse("--dir", "sets", "--pattern", "*.nwk")
        self.assertEqual(args.dir, "sets")
        self.assertEqual(args.pattern, "*.nwk")

    def test_single_dash_dir_is_not_valid(self):
        """-dir parses as -d ir, silently setting distance_matrix. Documented here."""
        args = self.parse("a.nwk", "-dir")
        self.assertEqual(args.distance_matrix, "ir")
        self.assertIsNone(args.dir)

    def test_version_flag_reports_the_package_version(self):
        """There was no --version at all, which is a reproducibility gap for a
        citable tool: nothing recorded which version produced a result."""
        import pear_ebi

        for flag in ("--version", "-V"):
            with self.subTest(flag=flag):
                buf = io.StringIO()
                with self.assertRaises(SystemExit) as ctx:
                    with redirect_stdout(buf):
                        self.parse(flag)
                self.assertEqual(ctx.exception.code, 0)
                self.assertIn(pear_ebi.__version__, buf.getvalue())

    def test_help_exits_zero(self):
        with self.assertRaises(SystemExit) as ctx:
            with redirect_stdout(io.StringIO()):
                self.parse("--help")
        self.assertEqual(ctx.exception.code, 0)


# ─── __main__.main() in-process ───────────────────────────────────────────────


class TestMainInProcess(unittest.TestCase):
    """Drives main() directly so coverage sees it and state can be inspected."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="pear_main_")
        self.cwd = self._tmp.name
        self._old_cwd = os.getcwd()
        for name in os.listdir(FIXTURES):
            shutil.copy(os.path.join(FIXTURES, name), self.cwd)
        os.chdir(self.cwd)

    def tearDown(self):
        os.chdir(self._old_cwd)
        self._tmp.cleanup()

    def main(self, *argv):
        """Run main() in-process. Returns (code, combined output).

        Library paths still raise SystemExit with a message rather than returning, so
        the message arrives as SystemExit.code and is folded into the output here.
        """
        from pear_ebi.__main__ import main

        from contextlib import redirect_stderr

        with mock.patch.object(sys, "argv", ["pear_ebi", *argv]):
            out, err = io.StringIO(), io.StringIO()
            with redirect_stdout(out), redirect_stderr(err):
                try:
                    code = main()
                except SystemExit as exc:
                    code = exc.code
        text = out.getvalue() + err.getvalue()
        if isinstance(code, str):
            text += code
            code = 1
        return code, text

    def test_no_arguments_is_a_usage_error(self):
        """bare exit() returned 0, so doing nothing looked like success."""
        code, output = self.main()
        self.assertEqual(code, 1)
        self.assertIn("No files specified", output)

    def test_missing_config_returns_nonzero(self):
        code, output = self.main("--config", "absent.toml")
        self.assertEqual(code, 1)
        self.assertIn("not found", output)

    def test_malformed_config_returns_nonzero(self):
        with open("broken.toml", "w") as fh:
            fh.write("not = valid = toml\n")
        code, output = self.main("--config", "broken.toml")
        self.assertEqual(code, 1)
        self.assertIn("Could not parse", output)

    def test_invalid_method_returns_nonzero(self):
        code, output = self.main("twelve_trees.nwk", "-m", "nope")
        self.assertNotEqual(code, 0)
        self.assertIn("Invalid method", output)

    def test_distance_run_writes_a_matrix(self):
        code, _ = self.main("twelve_trees.nwk", "-m", "hashrf_RF", "-o", "d.csv")
        self.assertIn(code, (0, None))
        self.assertTrue(os.path.exists("d.csv"))
        self.assertEqual(np.loadtxt("d.csv", delimiter=",").shape, (12, 12))

    def test_config_distance_method_is_used(self):
        with open("q.toml", "w") as fh:
            fh.write('[trees]\nfile1 = "twelve_trees.nwk"\n\n'
                     '[distance]\nmethod = "tqdist_quartet"\n')
        code, output = self.main("--config", "q.toml", "-o", "q.csv")
        self.assertIn(code, (0, None))
        self.assertIn("tqdist_quartet", output)

    def test_highlight_with_a_single_file(self):
        """A single [trees] file makes a tree_set, which has no .data attribute, so
        this documented feature used to die with AttributeError."""
        with open("h.toml", "w") as fh:
            fh.write('[trees]\nfile1 = "twelve_trees.nwk"\n\n'
                     '[distance]\nmethod = "hashrf_RF"\n\n'
                     '[embedding]\nmethod = "pcoa"\ndimensions = 2\n\n'
                     '[highlight]\ntrace1 = [0, 5]\n')
        code, output = self.main("--config", "h.toml")
        self.assertIn(code, (0, None), output[-800:])
        self.assertIn("2 highlighted", output)

    def test_highlight_with_a_collection(self):
        shutil.copy("twelve_trees.nwk", "second.nwk")
        with open("h2.toml", "w") as fh:
            fh.write('[trees]\nfile1 = "twelve_trees.nwk"\nfile2 = "second.nwk"\n\n'
                     '[distance]\nmethod = "hashrf_RF"\n\n'
                     '[embedding]\nmethod = "pcoa"\ndimensions = 2\n\n'
                     '[highlight]\ntrace1 = [0, 5]\ntrace2 = [1]\n')
        code, output = self.main("--config", "h2.toml")
        self.assertIn(code, (0, None), output[-800:])
        self.assertIn("highlighted", output)

    def test_highlight_indices_outside_the_range_are_ignored(self):
        """template_pear.toml ships trace1 = [0, 1001] for a 1001-tree set, where the
        valid range is 0-1000."""
        with open("h3.toml", "w") as fh:
            fh.write('[trees]\nfile1 = "twelve_trees.nwk"\n\n'
                     '[distance]\nmethod = "hashrf_RF"\n\n'
                     '[embedding]\nmethod = "pcoa"\ndimensions = 2\n\n'
                     '[highlight]\ntrace1 = [0, 9999]\n')
        code, output = self.main("--config", "h3.toml")
        self.assertIn(code, (0, None), output[-800:])
        self.assertIn("1 highlighted", output)

    def test_embedding_run(self):
        code, output = self.main("twelve_trees.nwk", "-m", "hashrf_RF", "--pcoa", "2")
        self.assertIn(code, (0, None))
        self.assertIn("pcoa", output)


# ─── subsample ────────────────────────────────────────────────────────────────


class TestSubsample(unittest.TestCase):
    """random.sample stopped accepting dict views in 3.11, killing this whole path."""

    def test_selects_the_requested_number(self):
        from pear_ebi.subsample import subsample

        trees, idxs = subsample.subsample(str([TWELVE]), 12, 4, subp=False)
        self.assertEqual(len(trees), 4)
        self.assertEqual(len(idxs), 4)
        self.assertEqual(len(set(idxs)), 4, "indices should be distinct")

    def test_refuses_more_than_available(self):
        from pear_ebi.subsample import subsample

        with self.assertRaises(ValueError):
            subsample.subsample(str([TWELVE]), 12, 999, subp=False)

    def test_result_marker_is_machine_readable(self):
        """The caller used to read stdout lines 3 and 4 and eval() them."""
        import json
        import subprocess

        from pear_ebi.subsample import subsample

        script = subsample.__file__
        result = subprocess.run(
            [sys.executable, script, str([TWELVE]), "12", "3"],
            capture_output=True, text=True, timeout=600,
        )
        self.assertEqual(result.returncode, 0, result.stderr[-600:])
        lines = [l for l in result.stdout.splitlines()
                 if l.startswith(subsample.RESULT_MARKER)]
        self.assertEqual(len(lines), 1, "expected exactly one result line")
        payload = json.loads(lines[0][len(subsample.RESULT_MARKER):])
        self.assertEqual(len(payload["trees"]), 3)
        self.assertEqual(len(payload["idxs"]), 3)


# ─── embedding quality ────────────────────────────────────────────────────────


class TestEmbeddingQuality(unittest.TestCase):
    def test_euclidean_distance_shape(self):
        from pear_ebi.embeddings.emb_quality import euclidean_distance

        points = np.array([[0.0, 0.0], [3.0, 4.0], [6.0, 8.0]])
        d = euclidean_distance(points)
        # every ordered pair including self-pairs -- documented, see the report
        self.assertEqual(d.shape, (9,))
        self.assertEqual(sorted(np.round(d, 6))[:3], [0.0, 0.0, 0.0])
        self.assertAlmostEqual(d.max(), 10.0)

    def test_pear_correlation_is_one_for_a_faithful_embedding(self):
        from pear_ebi.embeddings.emb_quality import pear_correlation

        pts = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
        corr = pear_correlation(pts, pts)
        self.assertAlmostEqual(corr[0, 1], 1.0, places=6)


# ─── maple_RF internals ───────────────────────────────────────────────────────


class TestMapleRF(unittest.TestCase):
    """maple_RF is a pure-Python RF implementation and the most testable module here."""

    def test_readNewick_parses_every_tree(self):
        from pear_ebi.calculate_distances import maple_RF

        with open(TWELVE) as fh:
            lines = [l for l in fh.read().splitlines() if l.strip()]
        trees = maple_RF.readNewick(lines)
        self.assertEqual(len(trees), len(lines))

    def test_identical_trees_have_zero_distance(self):
        from pear_ebi.calculate_distances import maple_RF

        with open(TWELVE) as fh:
            first = [fh.readline().strip()]
        t1 = maple_RF.readNewick(first)[0]
        t2 = maple_RF.readNewick(first)[0]
        prep = maple_RF.prepareTreeComparison(t1, rooted=False)
        self.assertEqual(maple_RF.RobinsonFouldsWithDay1985(t2, prep, rooted=False)[0], 0)

    def test_distance_is_symmetric(self):
        from pear_ebi.calculate_distances import maple_RF

        with open(TWELVE) as fh:
            lines = [l.strip() for l in fh.read().splitlines() if l.strip()][:4]
        trees = [maple_RF.readNewick([l])[0] for l in lines]
        for i in range(len(trees)):
            for j in range(i + 1, len(trees)):
                with self.subTest(pair=(i, j)):
                    a = maple_RF.RobinsonFouldsWithDay1985(
                        trees[j], maple_RF.prepareTreeComparison(trees[i], rooted=False),
                        rooted=False)[0]
                    b = maple_RF.RobinsonFouldsWithDay1985(
                        trees[i], maple_RF.prepareTreeComparison(trees[j], rooted=False),
                        rooted=False)[0]
                    self.assertEqual(a, b)

    def test_worker_is_picklable(self):
        """It was a closure assigned to a module global, so macOS spawn broke it."""
        import pickle

        from pear_ebi.calculate_distances import maple_RF

        self.assertTrue(pickle.dumps(maple_RF._rf_worker))


# ─── _install_helpers ─────────────────────────────────────────────────────────


class TestInstallHelpers(unittest.TestCase):
    def test_describe_platform_includes_the_architecture(self):
        """The bundled macOS binaries are arm64-only, so the arch belongs in messages."""
        from pear_ebi import _install_helpers as ih

        described = ih.describe_platform()
        self.assertIn(sys.platform, described)
        self.assertIn("/", described)

    def test_resolution_helpers_agree_with_each_other(self):
        from pear_ebi import _install_helpers as ih

        bin_dir = ih.platform_bin_dir()
        self.assertTrue(ih.hashrf_binary().startswith(bin_dir))
        self.assertTrue(ih.tqdist_bin_dir().startswith(bin_dir))
        self.assertEqual(len(ih.native_executables()), 9)

    def test_ensure_is_idempotent(self):
        from pear_ebi import _install_helpers as ih

        ih.ensure_native_executables()
        self.assertEqual(ih.ensure_native_executables(), [],
                         "a second call should have nothing left to fix")

    def test_build_tqdist_reports_rather_than_raises(self):
        """It returns (ok, message) so callers can surface the reason."""
        from pear_ebi import _install_helpers as ih

        ok, message = ih.build_tqdist(package_root="/nonexistent/package/root")
        self.assertFalse(ok)
        self.assertIsInstance(message, str)
        self.assertTrue(message)


# ─── interactive mode ─────────────────────────────────────────────────────────


class TestInteractiveMode(unittest.TestCase):
    def test_usage_lists_the_available_actions(self):
        from pear_ebi.interactive_mode import interactive

        with redirect_stdout(io.StringIO()) as out:
            interactive.usage()
        text = out.getvalue().lower()
        for expected in ("distance", "embed", "plot"):
            with self.subTest(action=expected):
                self.assertIn(expected, text)

    def test_interact_returns_executable_source_for_each_choice(self):
        """interact() returns strings that main() exec()s; check they at least compile."""
        from pear_ebi.interactive_mode import interactive

        compiled = 0
        for choice in range(1, 8):
            snippet = interactive.interact(str(choice))
            if isinstance(snippet, str) and snippet.strip():
                compile(snippet, "<interactive>", "exec")
                compiled += 1
        self.assertGreater(compiled, 0, "no choice produced compilable source")


if __name__ == "__main__":
    unittest.main()


# ─── native tool execution and error handling ─────────────────────────────────


class TestExecLayer(unittest.TestCase):
    """calculate_distances/_exec.py centralises what both native tools go through."""

    def test_both_streams_are_captured_and_labelled(self):
        from pear_ebi.calculate_distances._exec import run_process

        run = run_process(["/bin/sh", "-c", "echo OUT; echo ERR >&2; exit 3"])
        self.assertEqual(run.returncode, 3)
        self.assertFalse(run.ok)
        self.assertIn("OUT", run.stdout)
        self.assertIn("ERR", run.stderr)
        streams = run.streams()
        self.assertIn("stdout:", streams)
        self.assertIn("stderr:", streams)

    def test_streams_omits_empty_ones(self):
        from pear_ebi.calculate_distances._exec import format_streams

        self.assertNotIn("stdout:", format_streams("", "boom"))
        self.assertIn("no output on either stream", format_streams("", ""))

    def test_missing_binary_maps_to_127(self):
        from pear_ebi.calculate_distances._exec import RC_NOT_FOUND, run_process

        self.assertEqual(run_process(["/nonexistent/tool"]).returncode, RC_NOT_FOUND)

    def test_non_executable_binary_maps_to_126(self):
        from pear_ebi.calculate_distances._exec import RC_PERMISSION, run_process

        with tempfile.NamedTemporaryFile(suffix=".sh", delete=False) as fh:
            fh.write(b"#!/bin/sh\necho hi\n")
            path = fh.name
        try:
            os.chmod(path, 0o444)
            self.assertEqual(run_process([path]).returncode, RC_PERMISSION)
        finally:
            os.chmod(path, 0o644)
            os.unlink(path)

    def test_timeout_is_reported_not_raised(self):
        from pear_ebi.calculate_distances._exec import run_process

        run = run_process(["/bin/sh", "-c", "sleep 30"], timeout=1)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("timed out", run.stderr)

    def test_fingerprint_detects_a_fresh_write(self):
        """This is what separates "hashrf wrote this now" from a stale leftover."""
        from pear_ebi.calculate_distances._exec import (
            file_fingerprint,
            written_by_this_run,
        )

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.csv")
            before = file_fingerprint(path)
            self.assertFalse(written_by_this_run(path, before))
            with open(path, "w") as fh:
                fh.write("1,2\n")
            self.assertTrue(written_by_this_run(path, before))
            # unchanged since the new fingerprint
            self.assertFalse(written_by_this_run(path, file_fingerprint(path)))

    def test_remove_file_is_quiet_about_absent_files(self):
        from pear_ebi.calculate_distances._exec import remove_file

        self.assertTrue(remove_file("/nonexistent/path/file.csv"))

    def test_resolve_binary_prefers_path_but_warns(self):
        """A stale hashrf on PATH used to be picked silently."""
        from pear_ebi.calculate_distances import _exec

        with tempfile.TemporaryDirectory() as tmp:
            packaged = os.path.join(tmp, "mytool")
            with open(packaged, "w") as fh:
                fh.write("#!/bin/sh\n")
            os.chmod(packaged, 0o755)
            with mock.patch.object(_exec.shutil, "which", return_value="/usr/bin/mytool"):
                with self.assertWarns(RuntimeWarning):
                    resolved = _exec.resolve_binary("mytool", packaged, tool_label="MyTool")
            self.assertEqual(resolved, "/usr/bin/mytool")


class TestHashrfStdoutFallback(unittest.TestCase):
    """The last-resort scraper. It must decline rather than invent a matrix."""

    def test_declines_non_matrix_output(self):
        from pear_ebi.calculate_distances.hashrf import _parse_hashrf_stdout

        text = "*** Number of trees in the input file: 1\nFatal error\n"
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(
                _parse_hashrf_stdout(text, 5, os.path.join(tmp, "o.csv"))
            )

    def test_declines_empty_output(self):
        from pear_ebi.calculate_distances.hashrf import _parse_hashrf_stdout

        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(_parse_hashrf_stdout("", 3, os.path.join(tmp, "o.csv")))

    def test_parses_a_real_matrix_block(self):
        from pear_ebi.calculate_distances.hashrf import _parse_hashrf_stdout

        text = "Robinson-Foulds distance (matrix format):\n0 1 2\n1 0 3\n2 3 0\n"
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "o.csv")
            parsed = _parse_hashrf_stdout(text, 3, out)
            self.assertIsNotNone(parsed)
            self.assertEqual(np.asarray(parsed).shape, (3, 3))
            self.assertTrue(os.path.exists(out))


class TestDistanceMatrixLoading(unittest.TestCase):
    def test_rejects_a_non_square_matrix(self):
        from pear_ebi.tree_set import _load_distance_matrix

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "m.csv")
            np.savetxt(path, np.zeros((3, 5)), delimiter=",")
            with self.assertRaises(SystemExit):
                _load_distance_matrix(path, 3, "x.nwk")

    def test_rejects_a_size_mismatch(self):
        from pear_ebi.tree_set import _load_distance_matrix

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "m.csv")
            np.savetxt(path, np.zeros((4, 4)), delimiter=",")
            with self.assertRaises(SystemExit) as ctx:
                _load_distance_matrix(path, 12, "x.nwk")
            self.assertIn("does not match", str(ctx.exception))

    def test_accepts_an_in_memory_array(self):
        """DataFrame/array truthiness used to make this impossible."""
        from pear_ebi.tree_set import _load_distance_matrix

        m = np.zeros((3, 3))
        np.testing.assert_array_equal(_load_distance_matrix(m, 3, "x"), m)

    def test_none_passes_through(self):
        from pear_ebi.tree_set import _load_distance_matrix

        self.assertIsNone(_load_distance_matrix(None, 5, "x"))
