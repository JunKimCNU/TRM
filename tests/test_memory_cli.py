"""Linux CLI memory ceilings: evidence invariance and failed-resource isolation.

The ceiling limits virtual address space, not peak RSS. Resource exhaustion is
tested in subprocesses only, with core dumps disabled; parent/test-runner limits
are never changed. Run after building the release executable, or set TRM_BINARY.
"""

import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

if sys.platform == "linux":
    import resource


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]


def child_limits():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


@unittest.skipUnless(sys.platform == "linux", "RLIMIT_AS contract requires Linux")
class MemoryCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = runtime_binary()
        if not cls.binary.is_file():
            raise RuntimeError(f"Build the release executable first: {cls.binary}")

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="trm-memory-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.counter = 0
        self.fasta = self.root / "reads.fasta"
        flank = "ACGTCCGATGCTAACGTCAGTGCATCAGACGTTCGACGTA"
        sequences = ["TTAGGG" * 25, "CCCTAA" * 25,
                     flank + "TTAGGG" * 15, "CCCTAA" * 15 + flank,
                     "TTAGAC" * 25, "AAAAAAAG" * 19, "N" * 150] * 12
        self.fasta.write_text("".join(f">read{i}\n{sequence}\n" for i, sequence in enumerate(sequences)))

    def run_cli(self, action, *args, success=True, setup=child_limits):
        self.counter += 1
        prefix = self.root / f"result_{self.counter}"
        completed = subprocess.run(
            [str(self.binary), action, *map(str, args), "--output", str(prefix)],
            capture_output=True, text=True, timeout=40, preexec_fn=setup,
        )
        if not success:
            self.assertNotEqual(completed.returncode, 0, "invalid/resource-limited run succeeded")
            self.assertTrue(completed.stderr.strip(), "resource/configuration failure needs a diagnostic")
            for extension in (".json", ".tsv", ".rejected.tsv"):
                self.assertFalse(Path(str(prefix) + extension).exists(), completed.stderr)
            return completed, None
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return completed, json.loads(Path(str(prefix) + ".json").read_text())

    @staticmethod
    def observations(doc):
        summary = {key: value for key, value in doc["summary"].items()
                   if key not in ("elapsed_seconds", "megabases_per_second")}
        return {"summary": summary, **{key: doc[key] for key in
                ("config", "candidates", "rejected_candidates", "evidence_qc")}}

    def test_cap_and_smaller_batch_preserve_all_reported_observations(self):
        _, plain = self.run_cli("scan", self.fasta, "--threads", "1", "--min-reads", "1")
        _, capped = self.run_cli("scan", self.fasta, "--threads", "1", "--min-reads", "1",
                                 "--memory-limit-mib", "128")
        self.assertEqual(self.observations(plain), self.observations(capped))
        self.assertEqual(plain["memory_policy"]["enforcement"], "not_requested")
        policy = capped["memory_policy"]
        self.assertEqual(policy["requested_limit_mib"], 128)
        self.assertEqual(policy["effective_virtual_limit_bytes"], 128 * 1024 * 1024)
        self.assertEqual(policy["enforcement"], "linux_rlimit_as_virtual_address_space")
        self.assertGreater(policy["startup_virtual_bytes"], 0)
        self.assertLess(policy["effective_batch_bases"], policy["requested_batch_bases"])

    def test_tiny_limit_fails_before_opening_missing_input(self):
        completed, _ = self.run_cli("scan", self.root / "missing.fasta", "--threads", "1",
                                    "--memory-limit-mib", "1", success=False)
        self.assertIn("--memory-limit-mib", completed.stderr)
        self.assertNotIn("No such file", completed.stderr)

    def test_zero_negative_and_overflow_are_explicit_errors(self):
        for invalid in ("0", "-1", "18446744073709551615"):
            with self.subTest(value=invalid):
                completed, _ = self.run_cli("scan", self.fasta, "--threads", "1",
                                            "--memory-limit-mib", invalid, success=False)
                self.assertIn("--memory-limit-mib", completed.stderr)

    def test_reclassification_rejects_scan_memory_option(self):
        completed, _ = self.run_cli("reclassify", self.root / "missing.json",
                                    "--memory-limit-mib", "128", success=False)
        self.assertIn("reclassification", completed.stderr)

    def test_existing_stricter_os_limit_is_kept_and_reported(self):
        ceiling = 128 * 1024 * 1024

        def limit_child():
            child_limits()
            _, hard = resource.getrlimit(resource.RLIMIT_AS)
            resource.setrlimit(resource.RLIMIT_AS, (ceiling, hard))

        _, capped = self.run_cli("scan", self.fasta, "--threads", "1",
                                 "--memory-limit-mib", "256", setup=limit_child)
        self.assertEqual(capped["memory_policy"]["requested_limit_mib"], 256)
        self.assertEqual(capped["memory_policy"]["effective_virtual_limit_bytes"], ceiling)

    def test_single_record_exhaustion_is_failure_without_result(self):
        path = self.root / "record_larger_than_entire_address_space.fasta.gz"
        # Small compressed artifact, 128 MiB single FASTA line after decoding.
        # This must exceed a 64 MiB process ceiling before repeat detection runs.
        with gzip.open(path, "wb", compresslevel=1) as stream:
            stream.write(b">too_large\n")
            block = b"N" * (1024 * 1024)
            for _ in range(128):
                stream.write(block)
            stream.write(b"\n")
        completed, _ = self.run_cli("scan", path, "--threads", "1",
                                    "--memory-limit-mib", "64", success=False)
        self.assertIn("memory", completed.stderr.lower())
        self.assertIn("Resource exhaustion is a failed analysis", completed.stderr)
        self.assertFalse(list(self.root.glob("core*")))

    def test_mixed_cap_preserves_arrays_and_quality_evidence(self):
        definitions = self.root / "definitions.json"
        definitions.write_text(json.dumps([{"id": "two_units", "units": ["TTGGGTTAC", "ACTTGGGTTAC"]}]))
        fasta = self.root / "mixed.fasta"
        seq = ("TTGGGTTAC" + "ACTTGGGTTAC") * 12
        fasta.write_text(f">mixed\n{seq}\n>pure\n{'TTGGGTTAC' * 20}\n")
        common = (fasta, "--mixed-definitions", definitions, "--threads", "1")
        _, plain = self.run_cli("mixed", *common)
        _, capped = self.run_cli("mixed", *common, "--memory-limit-mib", "128")
        self.assertEqual(plain["mixed_arrays"], capped["mixed_arrays"])
        self.assertEqual(plain["evidence_qc"], capped["evidence_qc"])
        self.assertEqual(capped["memory_policy"]["effective_virtual_limit_bytes"], 128 * 1024 * 1024)


if __name__ == "__main__":
    unittest.main()
