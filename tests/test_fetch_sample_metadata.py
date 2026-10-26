import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("fetch_sample_metadata", Path(__file__).parents[1] / "scripts/fetch_sample_metadata.py")
metadata = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(metadata)


def run_response(accession="SRR123", tax_id="6239"):
    values = [accession, "SAMEA1", "PRJ1", "Caenorhabditis elegans", tax_id,
              "WGS", "GENOMIC", "PAIRED", "ILLUMINA", "NovaSeq"]
    return ("\t".join(metadata.FIELDS) + "\n" + "\t".join(values) + "\n").encode()


class MetadataTests(unittest.TestCase):
    def test_run_requires_exact_unique_accession(self):
        self.assertEqual(metadata.validate_run(run_response(), "SRR123")["tax_id"], "6239")
        with self.assertRaisesRegex(ValueError, "exactly requested"):
            metadata.validate_run(run_response("SRR999"), "SRR123")
        raw = run_response()
        with self.assertRaisesRegex(ValueError, "exactly requested"):
            metadata.validate_run(raw + raw.splitlines(keepends=True)[1], "SRR123")
        with self.assertRaisesRegex(ValueError, "positive integer"):
            metadata.validate_run(run_response(tax_id="invalid"), "SRR123")
        self.assertEqual(metadata.validate_run(run_response(tax_id=""), "SRR123")["tax_id"], "")

    def test_ncbi_fallback_requires_exact_run_and_explicit_taxonomy(self):
        raw = b"Run,TaxID,ScientificName,BioSample\nSRR123,6239,Caenorhabditis elegans,SAMEA1\n"
        self.assertEqual(metadata.validate_ncbi_run(raw, "SRR123")["TaxID"], "6239")
        with self.assertRaisesRegex(ValueError, "exactly requested"):
            metadata.validate_ncbi_run(raw, "SRR456")
        with self.assertRaisesRegex(ValueError, "valid TaxID"):
            metadata.validate_ncbi_run(raw.replace(b",6239,", b",,"), "SRR123")

    def test_taxonomy_requires_matching_identity_and_lineage(self):
        response = {"taxId": "6239", "scientificName": "Caenorhabditis elegans", "rank": "species",
                    "lineage": "Eukaryota; Metazoa; Nematoda;", "metagenome": "false"}
        raw = json.dumps(response).encode()
        self.assertEqual(metadata.validate_taxonomy(raw, "6239")["rank"], "species")
        with self.assertRaisesRegex(ValueError, "match requested"):
            metadata.validate_taxonomy(raw, "9606")
        response["merged"] = "[9606]"
        self.assertEqual(metadata.validate_taxonomy(json.dumps(response).encode(), "9606")["taxId"], "6239")
        response["merged"] = "[96060]"
        with self.assertRaisesRegex(ValueError, "match requested"):
            metadata.validate_taxonomy(json.dumps(response).encode(), "9606")
        del response["lineage"]
        with self.assertRaisesRegex(ValueError, "lineage"):
            metadata.validate_taxonomy(json.dumps(response).encode(), "6239")

    def test_cache_reuse_verifies_hash_url_and_response(self):
        calls = []
        def transport(url, timeout):
            calls.append(url)
            return run_response()
        with tempfile.TemporaryDirectory() as root:
            cache = metadata.ResponseCache(root, transport=transport, retries=0)
            first = cache.fetch("run", "SRR123", "https://example.test/run", metadata.validate_run)
            self.assertEqual(first, cache.fetch("run", "SRR123", "https://example.test/run", metadata.validate_run))
            self.assertEqual(len(calls), 1)
            with self.assertRaisesRegex(ValueError, "mismatch"):
                cache.fetch("run", "SRR123", "https://example.test/other", metadata.validate_run)
            path = Path(root) / "run_SRR123.response"
            path.write_bytes(run_response("SRR999"))
            with self.assertRaisesRegex(ValueError, "mismatch"):
                cache.fetch("run", "SRR123", "https://example.test/run", metadata.validate_run)
            self.assertEqual(path.read_bytes(), run_response("SRR999"))

    def test_invalid_response_is_not_cached_and_retry_is_bounded(self):
        calls, sleeps = [], []
        def transport(url, timeout):
            calls.append(url)
            return run_response("SRR999")
        with tempfile.TemporaryDirectory() as root:
            cache = metadata.ResponseCache(root, transport=transport, sleeper=sleeps.append, retries=2)
            with self.assertRaisesRegex(ValueError, "exactly requested"):
                cache.fetch("run", "SRR123", "https://example.test/run", metadata.validate_run)
            self.assertEqual(len(calls), 3)
            self.assertEqual(sleeps, [1, 2])
            self.assertEqual(list(Path(root).iterdir()), [])

    def test_manifest_duplicates_rejected_without_opening_read_paths(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "manifest.tsv"
            path.write_text("sample_id\tpath1\nSRR123\t/never/read.fastq\n")
            self.assertEqual(metadata.load_accessions(path), ["SRR123"])
            path.write_text(path.read_text() + "SRR123\t/also/never/read.fastq\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                metadata.load_accessions(path)

    def test_update_preserves_previous_result(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "result.tsv"
            metadata.write_result(path, b"old\n", False)
            with self.assertRaisesRegex(ValueError, "different result"):
                metadata.write_result(path, b"new\n", False)
            metadata.write_result(path, b"new\n", True)
            self.assertEqual(path.read_bytes(), b"new\n")
            self.assertEqual(next((Path(root) / "history").iterdir()).read_bytes(), b"old\n")


if __name__ == "__main__":
    unittest.main()
