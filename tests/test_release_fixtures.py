"""Byte-identity checks for relocated, independently frozen test fixtures."""
import hashlib
import json
from pathlib import Path
import unittest

FIXTURES = Path(__file__).resolve().parent / 'fixtures'


class ReleaseFixtureTests(unittest.TestCase):
    def test_every_relocated_fixture_preserves_its_source_bytes(self):
        manifest = json.loads((FIXTURES / 'provenance.json').read_text())
        self.assertEqual(manifest['schema_version'], 1)
        paths = [row['path'] for row in manifest['files']]
        self.assertEqual(len(paths), len(set(paths)))
        self.assertTrue(paths)
        for row in manifest['files']:
            with self.subTest(fixture=row['path']):
                relative = Path(row['path'])
                self.assertFalse(relative.is_absolute())
                self.assertNotIn('..', relative.parts)
                path = FIXTURES / relative
                self.assertFalse(path.is_symlink())
                content = path.read_bytes()
                self.assertEqual(len(content), row['bytes'])
                self.assertEqual(hashlib.sha256(content).hexdigest(), row['sha256'])
                self.assertEqual(row['transformation'], 'none; byte-identical copy')


if __name__ == '__main__':
    unittest.main()
