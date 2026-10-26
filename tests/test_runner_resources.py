"""Whole-runner thread planning must not change native results or resume identity."""
from __future__ import annotations

import contextlib
import importlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))

# Reuse fixtures by composition. Importing their modules (rather than exposing
# their TestCase classes here) keeps unittest discovery from repeating tests.
import test_batch_scaling as discovery_fixtures
import test_evidence_batch_scaling as evidence_fixtures


def load_script(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DISCOVERY = load_script("discovery_resource_tests", "batch.py")
EVIDENCE = evidence_fixtures.batch


class RunnerResourceMathTests(unittest.TestCase):
    def plan(self, **overrides):
        options = dict(jobs=4, workers=2, selected_samples=4, timed=True,
                       cpu_affinity=144, thread_budget=None, thread_reserve=0)
        options.update(overrides)
        return importlib.import_module("runner_resources").resource_plan(**options)

    def test_all_live_components_and_external_reserve_are_counted(self):
        plan = self.plan(thread_reserve=8)
        self.assertEqual(plan["effective_jobs"], 4)
        self.assertEqual(plan["native_threads"], 12)
        self.assertEqual(plan["python_threads"], 5)
        self.assertEqual(plan["gnu_time_threads"], 4)
        self.assertEqual(plan["runner_threads"], 21)
        self.assertEqual(plan["planned_total_threads"], 29)
        self.assertEqual(plan["thread_reserve"], 8)

    def test_untimed_discovery_removes_only_time_processes(self):
        plan = self.plan(timed=False)
        self.assertEqual(plan["gnu_time_threads"], 0)
        self.assertEqual(plan["runner_threads"], 17)
        self.assertEqual(plan["planned_total_threads"], 17)

    def test_explicit_budget_accepts_exact_boundary_and_rejects_next_job(self):
        at_limit = self.plan(jobs=27, selected_samples=100, thread_budget=144,
                             thread_reserve=8)
        over_limit = self.plan(jobs=28, selected_samples=100, thread_budget=144,
                               thread_reserve=8)
        self.assertEqual(at_limit["planned_total_threads"], 144)
        self.assertTrue(at_limit["within_budget"])
        self.assertEqual(over_limit["planned_total_threads"], 149)
        self.assertFalse(over_limit["within_budget"])
        self.assertEqual(at_limit["budget_mode"], "reject")
        self.assertEqual(at_limit["comparison_limit"], 144)
        self.assertEqual(at_limit["thread_budget"], 144)

    def test_default_budget_is_warning_against_affinity(self):
        plan = self.plan(cpu_affinity=8)
        self.assertIsNone(plan["thread_budget"])
        self.assertEqual(plan["comparison_limit"], 8)
        self.assertEqual(plan["budget_mode"], "warn")
        self.assertFalse(plan["within_budget"])

    def test_selection_caps_parallelism_without_changing_workers(self):
        plan = self.plan(jobs=100, workers=8, selected_samples=2)
        self.assertEqual(plan["effective_jobs"], 2)
        self.assertEqual(plan["native_threads"], 18)
        self.assertEqual(plan["runner_threads"], 23)

    def test_empty_selection_retains_main_and_external_reserve(self):
        plan = self.plan(selected_samples=0, thread_reserve=8)
        self.assertEqual(plan["effective_jobs"], 0)
        self.assertEqual(plan["native_threads"], 0)
        self.assertEqual(plan["python_threads"], 1)
        self.assertEqual(plan["gnu_time_threads"], 0)
        self.assertEqual(plan["runner_threads"], 1)
        self.assertEqual(plan["planned_total_threads"], 9)

    def test_explicit_limit_is_distinct_from_cpu_affinity(self):
        plan = self.plan(cpu_affinity=8, thread_budget=32)
        self.assertEqual(plan["comparison_limit"], 32)
        self.assertTrue(plan["within_budget"])

    def test_144_workers_does_not_fit_144_total_threads(self):
        for timed, total in ((False, 147), (True, 148)):
            with self.subTest(timed=timed):
                plan = self.plan(jobs=1, workers=144, selected_samples=1,
                                 timed=timed)
                self.assertEqual(plan["runner_threads"], total)
                self.assertFalse(plan["within_budget"])


class RunnerResourceCliTests(unittest.TestCase):
    def fixture(self, kind, names=("small", "larger")):
        if kind == "discovery":
            fixture = discovery_fixtures.BatchScalingTest("runTest")
        else:
            fixture = evidence_fixtures.EvidenceScalingTests("runTest")
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        if kind == "discovery":
            fixture.make_samples(names)
        else:
            fixture.samples(names)
        return fixture

    def arguments(self, kind, fixture, *extra):
        if kind == "discovery":
            script = DISCOVERY
            args = [str(fixture.manifest), "--output-dir", str(fixture.output)]
        else:
            script = EVIDENCE
            args = ["--manifest", str(fixture.manifest), "--source", str(fixture.source),
                    "--catalog", str(fixture.catalog), "--output", str(fixture.output)]
        return script, [str(script.__file__), *args, "--binary", str(fixture.binary),
                        "--jobs", "2", "--threads", "1", "--cpu-budget", "2", *extra]

    def invoke(self, kind, fixture, *extra, affinity=144, no_children=False):
        script, arguments = self.arguments(kind, fixture, *extra)
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch.object(sys, "argv", arguments))
            stack.enter_context(patch.object(script, "available_cpus", return_value=affinity))
            stack.enter_context(contextlib.redirect_stdout(stdout))
            stack.enter_context(contextlib.redirect_stderr(stderr))
            if no_children:
                stack.enter_context(patch("runtime_binary.resolve_binary",
                                          side_effect=AssertionError("binary resolution during preflight")))
                stack.enter_context(patch("subprocess.run",
                                          side_effect=AssertionError("child process during preflight")))
                stack.enter_context(patch("subprocess.Popen",
                                          side_effect=AssertionError("child process during preflight")))
            try:
                code = script.main()
            except SystemExit as exc:
                code = exc.code
        return code, stdout.getvalue(), stderr.getvalue()

    def read_summary(self, kind, fixture):
        return (fixture.summary() if kind == "discovery"
                else fixture.report("summary.json"))

    def test_strict_budget_rejects_before_binary_probe_or_output(self):
        for kind in ("discovery", "evidence"):
            with self.subTest(kind=kind):
                fixture = self.fixture(kind)
                code, stdout, stderr = self.invoke(kind, fixture, "--thread-budget", "1",
                                                   no_children=True)
                self.assertEqual(code, 2, stdout + stderr)
                self.assertIn("thread", stderr.lower())
                self.assertFalse(fixture.output.exists())
                self.assertEqual(fixture.launches(), [])

    def test_plan_only_accepts_missing_binary_without_outputs(self):
        for kind in ("discovery", "evidence"):
            with self.subTest(kind=kind):
                fixture = self.fixture(kind)
                fixture.binary.unlink()
                code, stdout, stderr = self.invoke(kind, fixture, "--resource-plan-only",
                                                   "--thread-budget", "32", "--thread-reserve", "4",
                                                   no_children=True)
                self.assertEqual(code, 0, stdout + stderr)
                plan = json.loads(stdout)
                self.assertEqual(plan["effective_jobs"], 2)
                self.assertEqual(plan["runner_threads"], 7 if kind == "discovery" else 9)
                self.assertEqual(plan["planned_total_threads"], 11 if kind == "discovery" else 13)
                self.assertEqual(plan["budget_mode"], "reject")
                self.assertTrue(plan["within_budget"])
                self.assertFalse(fixture.output.exists())
                self.assertEqual(fixture.launches(), [])

    def test_plan_only_strict_failure_still_prints_plan(self):
        for kind in ("discovery", "evidence"):
            with self.subTest(kind=kind):
                fixture = self.fixture(kind)
                code, stdout, stderr = self.invoke(kind, fixture, "--resource-plan-only",
                                                   "--thread-budget", "1", no_children=True)
                self.assertEqual(code, 2, stdout + stderr)
                self.assertFalse(json.loads(stdout)["within_budget"])
                self.assertFalse(fixture.output.exists())
                self.assertEqual(fixture.launches(), [])

    def test_discovery_resource_capture_counts_time_without_probing_it(self):
        fixture = self.fixture("discovery")
        code, stdout, stderr = self.invoke("discovery", fixture, "--resource-plan-only",
                                           "--record-resources", no_children=True)
        self.assertEqual(code, 0, stdout + stderr)
        plan = json.loads(stdout)
        self.assertEqual(plan["gnu_time_threads"], 2)
        self.assertEqual(plan["runner_threads"], 9)
        self.assertFalse(fixture.output.exists())

    def test_plan_only_warns_without_reducing_jobs_or_workers(self):
        for kind in ("discovery", "evidence"):
            with self.subTest(kind=kind):
                fixture = self.fixture(kind)
                code, stdout, stderr = self.invoke(kind, fixture, "--resource-plan-only",
                                                   affinity=2, no_children=True)
                self.assertEqual(code, 0, stdout + stderr)
                self.assertIn("warning", stderr.lower())
                self.assertIn("thread", stderr.lower())
                plan = json.loads(stdout)
                self.assertEqual(plan["effective_jobs"], 2)
                self.assertEqual(plan["native_threads"], 4)
                self.assertEqual(plan["comparison_limit"], 2)
                self.assertFalse(plan["within_budget"])
                self.assertFalse(fixture.output.exists())

    def test_discovery_default_workers_are_warned_not_silently_reduced(self):
        fixture = self.fixture("discovery", names=("only",))
        stdout, stderr = io.StringIO(), io.StringIO()
        arguments = [str(DISCOVERY.__file__), str(fixture.manifest), "--output-dir", str(fixture.output),
                     "--resource-plan-only"]
        with patch.object(sys, "argv", arguments), patch.object(DISCOVERY, "available_cpus", return_value=144), \
                patch("runtime_binary.resolve_binary", side_effect=AssertionError("binary resolution during preflight")), \
                patch("subprocess.run", side_effect=AssertionError("child during preflight")), \
                contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            self.assertEqual(DISCOVERY.main(), 0)
        plan = json.loads(stdout.getvalue())
        self.assertEqual(plan["native_threads"], 145)
        self.assertEqual(plan["planned_total_threads"], 147)
        self.assertFalse(plan["within_budget"])
        self.assertIn("warning", stderr.getvalue().lower())
        self.assertFalse(fixture.output.exists())

    def test_invalid_flags_reject_before_output(self):
        options = (("--thread-budget", "0"), ("--thread-budget", "-1"),
                   ("--thread-budget", "1.5"), ("--thread-reserve", "-1"),
                   ("--thread-reserve", "1.5"))
        for kind in ("discovery", "evidence"):
            fixture = self.fixture(kind)
            for extra in options:
                with self.subTest(kind=kind, extra=extra):
                    code, stdout, stderr = self.invoke(kind, fixture, *extra, no_children=True)
                    self.assertEqual(code, 2, stdout + stderr)
                    self.assertFalse(fixture.output.exists())
                    self.assertEqual(fixture.launches(), [])

    def test_worker_budget_stays_separate_from_total_budget(self):
        for kind in ("discovery", "evidence"):
            with self.subTest(kind=kind):
                fixture = self.fixture(kind)
                code, stdout, stderr = self.invoke(kind, fixture, "--cpu-budget", "1",
                                                   "--thread-budget", "100", "--resource-plan-only",
                                                   no_children=True)
                self.assertEqual(code, 2, stdout + stderr)
                self.assertIn("cpu-budget", stderr)
                self.assertFalse(fixture.output.exists())

    def test_shard_selection_changes_plan_before_empty_shard_validation(self):
        for kind in ("discovery", "evidence"):
            with self.subTest(kind=kind):
                fixture = self.fixture(kind, names=("only",))
                empty_index = 1 - discovery_fixtures.shard("only", 2)
                (fixture.root / "only.fq").unlink()
                code, stdout, stderr = self.invoke(kind, fixture, "--resource-plan-only",
                                                   "--shard-index", str(empty_index), "--shard-count", "2",
                                                   "--thread-budget", "4", "--thread-reserve", "3",
                                                   no_children=True)
                self.assertEqual(code, 0, stdout + stderr)
                plan = json.loads(stdout)
                self.assertEqual(plan["effective_jobs"], 0)
                self.assertEqual(plan["runner_threads"], 1)
                self.assertEqual(plan["planned_total_threads"], 4)
                self.assertTrue(plan["within_budget"])
                self.assertFalse(fixture.output.exists())
                self.assertEqual(fixture.launches(), [])

    def test_evidence_sample_selection_limits_plan(self):
        fixture = self.fixture("evidence")
        code, stdout, stderr = self.invoke("evidence", fixture, "--samples", "small",
                                           "--resource-plan-only", "--thread-budget", "5",
                                           no_children=True)
        self.assertEqual(code, 0, stdout + stderr)
        plan = json.loads(stdout)
        self.assertEqual(plan["effective_jobs"], 1)
        self.assertEqual(plan["planned_total_threads"], 5)
        self.assertFalse(fixture.output.exists())

    def test_changed_budget_does_not_change_results_commands_or_resume_identity(self):
        for kind in ("discovery", "evidence"):
            with self.subTest(kind=kind):
                fixture = self.fixture(kind)
                fixture.run_batch("--thread-budget", "32", jobs=2)
                summary = self.read_summary(kind, fixture)
                self.assertEqual(summary["resource_plan"]["thread_budget"], 32)
                self.assertEqual(summary["resource_plan"]["effective_jobs"], 2)
                self.assertEqual(len(fixture.launches()), 2)
                originals = {}
                for name in ("small", "larger"):
                    output = fixture.output / name
                    result = output / ("result.json" if kind == "discovery" else "evidence.json")
                    state = json.loads((output / ("state.json" if kind == "discovery" else "status.json")).read_text())
                    key = "fingerprint" if kind == "discovery" else "identity"
                    originals[name] = (result.read_bytes(), state[key])
                    for flag in ("--thread-budget", "--thread-reserve", "--resource-plan-only"):
                        self.assertNotIn(flag, json.dumps(state))
                fixture.run_batch("--resume", "--thread-budget", "48", "--thread-reserve", "4", jobs=2)
                summary = self.read_summary(kind, fixture)
                self.assertEqual(summary["resumed_samples"], 2)
                self.assertEqual(summary["resource_plan"]["thread_budget"], 48)
                # A cache hit must not cause an optimistic, zero-worker plan.
                self.assertEqual(summary["resource_plan"]["effective_jobs"], 2)
                self.assertEqual(summary["resource_plan"]["thread_reserve"], 4)
                self.assertEqual(len(fixture.launches()), 2)
                for name in ("small", "larger"):
                    output = fixture.output / name
                    result = output / ("result.json" if kind == "discovery" else "evidence.json")
                    state = json.loads((output / ("state.json" if kind == "discovery" else "status.json")).read_text())
                    key = "fingerprint" if kind == "discovery" else "identity"
                    self.assertEqual((result.read_bytes(), state[key]), originals[name])
                if kind == "evidence":
                    for filename in ("metrics.json", "run_metadata.json"):
                        self.assertEqual(fixture.report(filename)["resource_plan"], summary["resource_plan"])

    def test_default_affinity_warning_does_not_block_valid_stub_scans(self):
        for kind in ("discovery", "evidence"):
            with self.subTest(kind=kind):
                fixture = self.fixture(kind)
                code, stdout, stderr = self.invoke(kind, fixture, affinity=2)
                self.assertEqual(code, 0, stdout + stderr)
                self.assertIn("warning", stderr.lower())
                summary = self.read_summary(kind, fixture)
                self.assertTrue(summary["all_succeeded"])
                self.assertEqual(summary["resource_plan"]["budget_mode"], "warn")
                self.assertFalse(summary["resource_plan"]["within_budget"])
                self.assertEqual(len(fixture.launches()), 2)

    def test_strict_resume_rejection_preserves_cached_and_corrupt_outputs(self):
        for kind in ("discovery", "evidence"):
            fixture = self.fixture(kind)
            fixture.run_batch("--thread-budget", "32", jobs=2)
            self.assertEqual(len(fixture.launches()), 2)
            # This budget would fit one active sample, but not both selected
            # samples. Neither an intact cache nor one corrupt cache entry may
            # cause preflight to plan only the presumed cache misses.
            one_sample_budget = 4 if kind == "discovery" else 5
            for corrupt in (False, True):
                with self.subTest(kind=kind, corrupt_cache=corrupt):
                    if corrupt:
                        filename = "result.json" if kind == "discovery" else "evidence.json"
                        (fixture.output / "small" / filename).write_text("{}")

                    def snapshot():
                        return {
                            str(path.relative_to(fixture.output)):
                                (path.read_bytes(), path.stat().st_mtime_ns)
                            for path in fixture.output.rglob("*") if path.is_file()
                        }

                    before = snapshot()
                    code, stdout, stderr = self.invoke(
                        kind, fixture, "--resume", "--thread-budget", str(one_sample_budget),
                        no_children=True)
                    self.assertEqual(code, 2, stdout + stderr)
                    self.assertIn("thread-budget", stderr)
                    self.assertEqual(snapshot(), before)
                    self.assertEqual(len(fixture.launches()), 2)


if __name__ == "__main__":
    unittest.main()
