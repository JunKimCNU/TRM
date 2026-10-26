import csv
import importlib.util
import io
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("inventory", Path(__file__).parents[1] / "scripts/inventory_public_wgs.py")
inventory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inventory)


def response(**overrides):
    row = {field: "" for field in inventory.FIELDS}
    row.update(run_accession="SRR123", tax_id="6239", scientific_name="C. elegans", library_strategy="WGS",
               library_source="GENOMIC", instrument_platform="ILLUMINA", fastq_bytes="10;20",
               fastq_ftp="ftp/a.fastq.gz;ftp/b.fastq.gz", fastq_md5="abc;def", base_count="1000")
    row.update(overrides)
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=inventory.FIELDS, delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
    return stream.getvalue().encode()


class InventoryTests(unittest.TestCase):
    def test_count_is_not_row_count_or_species_count(self):
        self.assertEqual(inventory.parse_count(b'{"count":"3622076"}'), 3622076)
        for raw in (b'{"count":-1}', b'[]', b'{"count":null}'):
            with self.assertRaises(ValueError):
                inventory.parse_count(raw)

    def test_group_count_rejects_duplicate_or_malformed(self):
        self.assertEqual(inventory.parse_group_counts(b"source\tcount\ngenomic\t12\n", "source"), {"genomic": 12})
        with self.assertRaises(ValueError):
            inventory.parse_group_counts(b"source\tcount\na\t1\na\t2\n", "source")

    def test_sizes_are_summed_with_unknown_preserved(self):
        self.assertEqual(inventory.byte_sum("10;20"), 30)
        self.assertIsNone(inventory.byte_sum(""))
        self.assertEqual(inventory.byte_sum("0"), 0)
        for value in ("1;", "-1", "1.5", "12;NA"):
            with self.assertRaises(ValueError):
                inventory.byte_sum(value)

    def test_paired_bytes_and_archive_formats_are_separate(self):
        row = inventory.parse_runs(response(sra_bytes="5", submitted_bytes="7;8"), 1, "other_metazoa")[0]
        self.assertEqual(row["fastq_total_bytes"], "30")
        self.assertEqual(row["sra_total_bytes"], "5")
        self.assertEqual(row["submitted_total_bytes"], "15")

    def test_missing_bytes_not_zero(self):
        row = inventory.parse_runs(response(fastq_bytes="", fastq_ftp="", fastq_md5=""), 1, "other_metazoa")[0]
        self.assertEqual(row["fastq_total_bytes"], "")
        summary = inventory.summarize([row], {"wgs_genomic": 100}, {})
        self.assertEqual(summary["pilot_only"]["fastq_bytes_missing_runs"], 1)
        self.assertIsNone(summary["full_archive_data_bytes"])
        self.assertIsNone(summary["full_archive_species_count"])

    def test_wrong_query_data_rejected(self):
        with self.assertRaises(ValueError):
            inventory.parse_runs(response(library_source="TRANSCRIPTOMIC"), 1, "other_metazoa")
        with self.assertRaises(ValueError):
            inventory.parse_runs(response(run_accession="SRX123"), 1, "other_metazoa")

    def test_file_metadata_pairing_checked(self):
        with self.assertRaises(ValueError):
            inventory.parse_runs(response(fastq_bytes="10"), 1, "other_metazoa")
        with self.assertRaises(ValueError):
            inventory.parse_runs(response(fastq_md5="abc"), 1, "other_metazoa")

    def test_cap_and_duplicate_enforced(self):
        raw = response()
        with self.assertRaises(ValueError):
            inventory.parse_runs(raw, 0, "other_metazoa")
        with self.assertRaises(ValueError):
            inventory.parse_runs(raw + raw.splitlines(keepends=True)[1], 2, "other_metazoa")

    def test_taxon_ids_not_automatically_species(self):
        row = inventory.parse_runs(response(), 1, "other_metazoa")[0]
        self.assertEqual(inventory.summarize([row], {}, {})["pilot_only"]["rank_species_tax_ids"], 0)
        row["taxonomy_rank"] = "species"
        self.assertEqual(inventory.summarize([row], {}, {})["pilot_only"]["rank_species_tax_ids"], 1)


if __name__ == "__main__":
    unittest.main()
