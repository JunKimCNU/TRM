import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "split_control", Path(__file__).resolve().parents[1] / "scripts" / "split_control.py")
split = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(split)


class SplitControlTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)

    def pair(self, count=100):
        paths = [self.root / "R1.fastq.gz", self.root / "R2.fastq.gz"]
        for mate, path in enumerate(paths, 1):
            with gzip.open(path, "wb") as handle:
                for index in range(count):
                    handle.write(f"@read{index}/{mate}\nACGTN\n+\nIIIII\n".encode())
        return paths

    def ids(self, path, mate):
        return [record[0] for record in split.records(path, mate)]

    def test_partition_pair_consistency_and_reproducibility(self):
        first, second = self.pair()
        outputs = [self.root / "one", self.root / "two"]
        results = [split.split_control(first, second, path, "42") for path in outputs]
        all_ids = set()
        for part in ("calibration", "heldout"):
            ids = self.ids(outputs[0] / f"{part}_R1.fastq.gz", 1)
            self.assertEqual(ids, self.ids(outputs[0] / f"{part}_R2.fastq.gz", 2))
            self.assertTrue(ids)
            self.assertFalse(all_ids.intersection(ids))
            all_ids.update(ids)
            for mate in (1, 2):
                self.assertEqual((outputs[0] / f"{part}_R{mate}.fastq.gz").read_bytes(),
                                 (outputs[1] / f"{part}_R{mate}.fastq.gz").read_bytes())
        self.assertEqual(len(all_ids), 100)
        self.assertEqual(results[0]["total_records"], 200)
        self.assertEqual(results[0]["total_bases"], 1000)
        third = self.root / "three"
        split.split_control(first, second, third, "43")
        self.assertNotEqual(self.ids(outputs[0] / "calibration_R1.fastq.gz", 1),
                            self.ids(third / "calibration_R1.fastq.gz", 1))

    def test_sra_and_casava_id_normalization(self):
        self.assertEqual(split.normalized_id(b"@SRR1.123 1/1", 1),
                         split.normalized_id(b"@SRR1.123 1/2", 2))
        self.assertEqual(split.normalized_id(b"@machine:123 1:N:0:ATGC", 1), b"machine:123")
        with self.assertRaisesRegex(ValueError, "Wrong mate"):
            split.normalized_id(b"@read/2", 1)

    def test_mismatched_ids_fail_with_provenance(self):
        first, second = self.pair(1)
        with gzip.open(second, "wb") as handle:
            handle.write(b"@other/2\nACGTN\n+\nIIIII\n")
        out = self.root / "bad"
        with self.assertRaisesRegex(ValueError, "Mate ID mismatch"):
            split.split_control(first, second, out, "42")
        self.assertEqual(json.loads((out / "provenance.json").read_text())["status"], "failed")
        self.assertFalse(list(out.glob("*.fastq.gz")))

    def test_count_mismatch(self):
        first, second = self.pair(2)
        with gzip.open(second, "wb") as handle:
            handle.write(b"@read0/2\nACGTN\n+\nIIIII\n")
        with self.assertRaisesRegex(ValueError, "count mismatch"):
            split.split_control(first, second, self.root / "bad", "42")

    def test_incomplete_fastq_and_length_mismatch(self):
        for index, malformed in enumerate((b"@read0/2\nACGTN\n+\n", b"@read0/2\nACGTN\n+\nIII\n")):
            first, second = self.pair(1)
            with gzip.open(second, "wb") as handle:
                handle.write(malformed)
            with self.assertRaises(ValueError):
                split.split_control(first, second, self.root / f"bad{index}", "42")

    def test_truncated_gzip_rejected(self):
        first, second = self.pair(100)
        second.write_bytes(second.read_bytes()[:-6])
        with self.assertRaises((EOFError, OSError)):
            split.split_control(first, second, self.root / "bad", "42")

    def test_manifest_verification(self):
        first, second = self.pair(10)
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({"complete": True, "controls": {
            "test": {"files": [split.file_info(first), split.file_info(second)]}}}))
        result = split.split_control(first, second, self.root / "good", "42", manifest=manifest)
        self.assertEqual(result["status"], "complete")
        with gzip.open(second, "ab") as handle:
            handle.write(b"@extra/2\nA\n+\nI\n")
        with self.assertRaisesRegex(ValueError, "differs from manifest"):
            split.split_control(first, second, self.root / "bad", "42", manifest=manifest)
        self.assertFalse((self.root / "bad").exists())

    def test_existing_directory_and_invalid_fraction_rejected(self):
        first, second = self.pair(1)
        with self.assertRaises(FileExistsError):
            split.split_control(first, second, self.root, "42")
        for fraction in (0, 1, -0.1, float("nan")):
            with self.assertRaises(ValueError):
                split.split_control(first, second, self.root / "bad", "42", fraction)


if __name__ == "__main__":
    unittest.main()
