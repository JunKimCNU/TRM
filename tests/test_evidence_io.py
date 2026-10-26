"""Adversarial I/O cases shared by both saved-evidence report overlays."""
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import cross_period_evidence as cross
import reference_localization_evidence as reference
from report_code_dependencies import local_dependency_paths

MODULES = (
    (cross, 'cross-period evidence: '),
    (reference, 'reference localization evidence: '),
)


class EvidenceIOTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.path = self.root / 'evidence.json'

    def test_payload_seal_preserves_nested_seal_and_unicode(self):
        # The canonical bytes are specified independently of json.dumps.
        canonical = '{"nested":{"payload_sha256":"keep"},"표본":"종"}'.encode('utf-8')
        expected = hashlib.sha256(canonical).hexdigest()
        for module, _ in MODULES:
            with self.subTest(module=module.__name__):
                document = {'표본': '종', 'payload_sha256': 'discard',
                            'nested': {'payload_sha256': 'keep'}}
                self.assertEqual(module.evidence_payload_sha256(document), expected)
                document['payload_sha256'] = 'another top-level seal'
                self.assertEqual(module.evidence_payload_sha256(document), expected)
                document['nested']['payload_sha256'] = 'changed'
                self.assertNotEqual(module.evidence_payload_sha256(document), expected)

    def test_invalid_payload_types_and_encoding_preserve_error_causes(self):
        for module, prefix in MODULES:
            for value, cause in ((float('nan'), ValueError),
                                 (object(), TypeError),
                                 ('\ud800', UnicodeEncodeError)):
                with self.subTest(module=module.__name__, cause=cause.__name__):
                    with self.assertRaises(module.EvidenceError) as caught:
                        module.evidence_payload_sha256({'value': value})
                    self.assertEqual(str(caught.exception), prefix + 'non-JSON payload')
                    self.assertIsInstance(caught.exception.__cause__, cause)
            with self.assertRaises(module.EvidenceError) as caught:
                module.evidence_payload_sha256([])
            self.assertEqual(str(caught.exception), prefix + 'document must be an object')
            self.assertIsNone(caught.exception.__cause__)

    def test_strict_json_errors_keep_overlay_type_and_detailed_cause(self):
        cases = (
            (b'{"value":{"x":1,"x":2}}', "duplicate JSON key 'x'"),
            (b'{"value":NaN}', 'non-finite JSON number NaN'),
            (b'{"value":Infinity}', 'non-finite JSON number Infinity'),
            (b'{"value":-Infinity}', 'non-finite JSON number -Infinity'),
        )
        for module, prefix in MODULES:
            for raw, detail in cases:
                with self.subTest(module=module.__name__, raw=raw):
                    self.path.write_bytes(raw)
                    with self.assertRaises(module.EvidenceError) as caught:
                        module._read_json(self.path)
                    self.assertEqual(str(caught.exception), prefix + 'malformed JSON in ' + str(self.path))
                    self.assertIsInstance(caught.exception.__cause__, module.EvidenceError)
                    self.assertEqual(str(caught.exception.__cause__), prefix + detail)
            for raw, cause in ((b'{"x":"\xff"}', UnicodeDecodeError),
                               (b'{"x":', json.JSONDecodeError)):
                self.path.write_bytes(raw)
                with self.assertRaises(module.EvidenceError) as caught:
                    module._read_json(self.path)
                self.assertIsInstance(caught.exception.__cause__, cause)

    def test_raw_json_hash_and_top_level_type_have_distinct_behavior(self):
        raw = b'{ "x": 1 }\n'
        self.path.write_bytes(raw)
        for module, _ in MODULES:
            with self.subTest(module=module.__name__):
                self.assertEqual(module._read_json(self.path),
                                 ({'x': 1}, hashlib.sha256(raw).hexdigest()))
        self.path.write_bytes(b'[]')
        for module, prefix in MODULES:
            with self.assertRaises(module.EvidenceError) as caught:
                module._read_json(self.path)
            self.assertEqual(str(caught.exception), prefix + 'document must be an object')
            self.assertIsNone(caught.exception.__cause__)

    def test_runtime_cap_boundary_regular_file_and_os_error_contracts(self):
        self.path.write_bytes(b'12345678')
        for module, prefix in MODULES:
            with self.subTest(module=module.__name__):
                with mock.patch.object(module, 'MAX_ARTIFACT_BYTES', 8):
                    self.assertEqual(module._read_bytes(str(self.path)), b'12345678')
                with mock.patch.object(module, 'MAX_ARTIFACT_BYTES', 7):
                    with self.assertRaises(module.EvidenceError) as caught:
                        module._read_bytes(self.path)
                self.assertEqual(str(caught.exception), prefix + 'source must be a small regular file: ' + str(self.path))
                self.assertIsNone(caught.exception.__cause__)
                with self.assertRaises(module.EvidenceError) as caught:
                    module._read_bytes(self.root)
                self.assertIn('small regular file', str(caught.exception))
                missing = self.root / 'missing.json'
                with self.assertRaises(module.EvidenceError) as caught:
                    module._read_bytes(missing)
                self.assertEqual(str(caught.exception), prefix + 'cannot read ' + str(missing))
                self.assertIsInstance(caught.exception.__cause__, FileNotFoundError)

    def test_descriptor_substitution_is_rejected_when_path_is_unchanged(self):
        self.path.write_bytes(b'original')
        other = self.root / 'other.json'
        other.write_bytes(b'replaced')
        real_open = Path.open

        def substituted(path, *args, **kwargs):
            return real_open(other if path == self.path else path, *args, **kwargs)

        for module, prefix in MODULES:
            with self.subTest(module=module.__name__), mock.patch.object(Path, 'open', substituted):
                with self.assertRaises(module.EvidenceError) as caught:
                    module._read_bytes(self.path)
                self.assertEqual(str(caught.exception), prefix + 'source changed while reading: ' + str(self.path))
                self.assertIsNone(caught.exception.__cause__)
        self.assertEqual(self.path.read_bytes(), b'original')

    def test_symlink_retarget_is_rejected_even_with_identical_inode(self):
        self.path.write_bytes(b'unchanged')
        alias = self.root / 'same_inode.json'
        os.link(self.path, alias)
        link = self.root / 'link.json'
        real_open = Path.open

        def retargeted(path, *args, **kwargs):
            if path == link:
                link.unlink()
                link.symlink_to(alias)
            return real_open(path, *args, **kwargs)

        for module, prefix in MODULES:
            if link.is_symlink():
                link.unlink()
            link.symlink_to(self.path)
            with self.subTest(module=module.__name__), mock.patch.object(Path, 'open', retargeted):
                with self.assertRaises(module.EvidenceError) as caught:
                    module._read_bytes(link)
                self.assertEqual(str(caught.exception), prefix + 'source changed while reading: ' + str(link))
                self.assertIsNone(caught.exception.__cause__)

    def test_shared_io_participates_in_report_cache_fingerprint(self):
        scripts = Path(__file__).resolve().parents[1] / 'scripts'
        self.assertIn(scripts / 'evidence_io.py', local_dependency_paths(scripts))


if __name__ == '__main__':
    unittest.main()
