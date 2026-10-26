"""Numeric selector and argument/file-identity tests only; no FASTQ fixtures.

Real-source full EOF, output contents, CRC and failure integration checks are
performed separately against existing deposited sequencing data.
"""

from decimal import localcontext
import hashlib
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "subsample_read_pairs", Path(__file__).resolve().parents[1] / "scripts" / "subsample_read_pairs.py")
subsample = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(subsample)


class OrdinalSelectorTests(unittest.TestCase):
    def test_documented_fixed_width_hash_message(self):
        seed, ordinal = 1729, 1
        message = b"trm-pair-ordinal-subsample-v1\0" + seed.to_bytes(8, "big") + ordinal.to_bytes(8, "big")
        expected = int.from_bytes(hashlib.sha256(message).digest()[:8], "big")
        self.assertEqual(subsample.ordinal_value(seed, ordinal), expected)
        self.assertEqual(subsample.ordinal_value(seed, ordinal), subsample.ordinal_value(seed, ordinal))

    def test_deterministic_for_numeric_ordinals(self):
        first = [subsample.ordinal_value(2718, n) for n in range(1, 300)]
        second = [subsample.ordinal_value(2718, n) for n in range(1, 300)]
        self.assertEqual(first, second)
        self.assertNotEqual(first, [subsample.ordinal_value(31415, n) for n in range(1, 300)])
        self.assertEqual(len(set(first)), len(first))

    def test_nested_subsets_within_seed(self):
        thresholds = [subsample.fraction_spec(value)["threshold"] for value in ("0.0001", "0.001", "0.01", "0.1", "1")]
        values = {n: subsample.ordinal_value(1729, n) for n in range(1, 20_001)}
        subsets = [{n for n, value in values.items() if subsample.selected(value, threshold)} for threshold in thresholds]
        for small, large in zip(subsets, subsets[1:]):
            self.assertLessEqual(small, large)
        self.assertEqual(len(subsets[-1]), 20_000)
        self.assertNotEqual(subsets[2], subsets[3])

    def test_domain_separation(self):
        default = subsample.ordinal_value(1729, 1)
        self.assertNotEqual(default, subsample.ordinal_value(1729, 1, b"another-domain\0"))
        self.assertNotEqual(default, subsample.ordinal_value(1, 1729))
        self.assertNotEqual(subsample.seed_prefix(1), subsample.seed_prefix(256))
        for domain in (b"", "not bytes", None):
            with self.subTest(domain=domain), self.assertRaises(ValueError):
                subsample.seed_prefix(1729, domain)

    def test_exact_threshold_is_independent_of_decimal_precision(self):
        with localcontext() as context:
            context.prec = 3
            self.assertEqual(subsample.fraction_spec("0.1")["threshold"], (1 << 64) // 10)
            self.assertEqual(subsample.fraction_spec("0.0001")["threshold"], (1 << 64) // 10_000)
            self.assertEqual(subsample.fraction_spec("0.99999999999999999999")["threshold"], (1 << 64) - 1)

    def test_boundary_strictly_less_than_threshold(self):
        threshold = subsample.fraction_spec("0.1")["threshold"]
        self.assertTrue(subsample.selected(threshold - 1, threshold))
        self.assertFalse(subsample.selected(threshold, threshold))
        self.assertFalse(subsample.selected(threshold + 1, threshold))
        self.assertFalse(subsample.selected(0, 0))
        self.assertTrue(subsample.selected((1 << 64) - 1, 1 << 64))

    def test_smallest_representable_interval(self):
        # This decimal is slightly larger than 2**-64, and selects exactly zero.
        self.assertEqual(subsample.fraction_spec("0.00000000000000000005421010862427523")["threshold"], 1)
        self.assertTrue(subsample.selected(0, 1))
        self.assertFalse(subsample.selected(1, 1))

    def test_fraction_labels_canonicalize_equivalent_spellings(self):
        self.assertEqual(subsample.fraction_spec("1e-3"), subsample.fraction_spec("0.001000"))
        self.assertEqual(subsample.fraction_spec("0.100000")["label"], "0.1")
        self.assertEqual(subsample.fraction_spec("1.000")["label"], "1")
        self.assertEqual(subsample.fraction_spec("1")["threshold"], 1 << 64)

    def test_invalid_fractions(self):
        for value in ("0", "-0.1", "1.1", "NaN", "sNaN", "Infinity", "-Infinity", "bad", "", "1e-100000", "1e-20", True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                subsample.fraction_spec(value)

    def test_seed_and_ordinal_uint64_limits(self):
        self.assertEqual(subsample.uint64(0), 0)
        self.assertEqual(subsample.uint64(str((1 << 64) - 1)), (1 << 64) - 1)
        self.assertIsInstance(subsample.ordinal_value((1 << 64) - 1, (1 << 64) - 1), int)
        for value in (-1, 1 << 64, True, 1.0, "1.0", "1e3", "-0", ""):
            with self.subTest(value=value), self.assertRaises(ValueError):
                subsample.uint64(value)
        with self.assertRaises(ValueError):
            subsample.ordinal_value(0, 0)

    def test_invalid_selection_values_and_thresholds(self):
        for value in (-1, 1 << 64, True, "3", 3.0):
            with self.subTest(value=value), self.assertRaises(ValueError):
                subsample.selected(value, 10)
        for threshold in (-1, (1 << 64) + 1, False, "10", 10.0):
            with self.subTest(threshold=threshold), self.assertRaises(ValueError):
                subsample.selected(3, threshold)


class ParameterTests(unittest.TestCase):
    def arguments(self, **overrides):
        values = {"sample_id": "archive-123", "fractions": ["0.1", "0.001"], "seeds": [1729, 2718],
                  "expected_pairs": 100, "expected_bases": 30000,
                  "r1_sha256": "A" * 64, "r2_sha256": "b" * 64}
        values.update(overrides)
        return values

    def test_valid_parameters_are_normalized(self):
        fractions, seeds, pairs, bases, hashes = subsample.validate_parameters(**self.arguments())
        self.assertEqual([item["label"] for item in fractions], ["0.001", "0.1"])
        self.assertEqual(seeds, [1729, 2718])
        self.assertEqual((pairs, bases), (100, 30000))
        self.assertEqual(hashes, ["a" * 64, "b" * 64])

    def test_duplicate_empty_or_invalid_arguments(self):
        cases = ({"fractions": []}, {"fractions": ["0.1", "1e-1"]},
                 {"seeds": []}, {"seeds": ["01729", 1729]},
                 {"sample_id": "../outside"}, {"sample_id": "name\tvalue"},
                 {"sample_id": ""}, {"sample_id": "."}, {"sample_id": ".."},
                 {"expected_pairs": 0}, {"expected_bases": 0},
                 {"expected_bases": 199}, {"r1_sha256": "x" * 64},
                 {"r2_sha256": "a" * 63}, {"r1_sha256": None})
        for overrides in cases:
            with self.subTest(overrides=overrides), self.assertRaises(ValueError):
                subsample.validate_parameters(**self.arguments(**overrides))

    def test_parser_retains_exact_decimal_strings(self):
        args = subsample.build_parser().parse_args([
            "--sample-id", "example", "--r1", "first", "--r2", "second", "--output-dir", "new",
            "--fractions", "0.1000000000000000001", "1e-3", "--seeds", "1729", "2718",
            "--expected-pairs", "100", "--expected-bases", "30000",
            "--r1-sha256", "a" * 64, "--r2-sha256", "b" * 64])
        self.assertEqual(args.fractions, ["0.1000000000000000001", "1e-3"])
        self.assertEqual(args.seeds, ["1729", "2718"])


class SourceIdentityTests(unittest.TestCase):
    """Empty files test OS identity only; none are passed to a FASTQ reader."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.first = self.root / "first.empty"
        self.first.touch()

    def test_same_file_rejected(self):
        with self.assertRaisesRegex(ValueError, "distinct files"):
            subsample.validate_source_paths(self.first, self.first)

    def test_symlink_alias_rejected(self):
        alias = self.root / "symlink.empty"
        alias.symlink_to(self.first)
        with self.assertRaisesRegex(ValueError, "distinct files"):
            subsample.validate_source_paths(self.first, alias)

    def test_hardlink_alias_rejected(self):
        alias = self.root / "hardlink.empty"
        os.link(self.first, alias)
        with self.assertRaisesRegex(ValueError, "distinct files"):
            subsample.validate_source_paths(self.first, alias)

    def test_distinct_files_allowed(self):
        second = self.root / "second.empty"
        second.touch()
        paths, snapshots = subsample.validate_source_paths(self.first, second)
        self.assertEqual(paths, [self.first.resolve(), second.resolve()])
        self.assertNotEqual(snapshots[0]["inode"], snapshots[1]["inode"])
        self.assertEqual(snapshots[0]["bytes"], 0)

    def test_directory_rejected(self):
        with self.assertRaisesRegex(ValueError, "regular file"):
            subsample.validate_source_paths(self.first, self.root)


if __name__ == "__main__":
    unittest.main()
