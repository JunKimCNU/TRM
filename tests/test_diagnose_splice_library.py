"""Coordinate, strand and alignment-denominator checks for the diagnostic."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('splice_diagnostic', Path(__file__).resolve().parents[1] / 'scripts/diagnose_splice_library.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SpliceDiagnosticTests(unittest.TestCase):
    def test_cigar_consumption_and_anchor_boundaries(self):
        blocks, skips = MODULE.cigar_events('3S12M40N15M2I5M4D6M', 100)
        self.assertEqual(blocks, [(100, 112, 3, 15), (152, 167, 15, 30),
                                  (167, 172, 32, 37), (176, 182, 37, 43)])
        self.assertEqual(skips, [(112, 152, 15, 12, 15)])
        _, skips = MODULE.cigar_events('12M1I40N15M', 0)
        self.assertEqual(skips[0][3], 0)

    def test_canonical_exact_anchors_and_annotation_strand(self):
        genome = 'A' * 20 + 'GT' + 'C' * 36 + 'AG' + 'T' * 20
        read = 'A' * 20 + 'T' * 20
        skip = (20, 60, 20, 20, 20)
        details = MODULE.junction_details('c', genome, read, skip, {('c', 20, 60, '+')})
        self.assertEqual(details, ('GTAG', '+', True, True, True))
        mismatch = read[:19] + 'C' + read[20:]
        self.assertFalse(MODULE.junction_details('c', genome, mismatch, skip, set())[-1])
        short_anchor = (20, 60, 20, 9, 20)
        self.assertFalse(MODULE.junction_details('c', genome, read, short_anchor, set())[-1])
        ambiguous_genome = genome[:19] + 'N' + genome[20:]
        ambiguous_read = read[:19] + 'N' + read[20:]
        self.assertFalse(MODULE.junction_details('c', ambiguous_genome, ambiguous_read, skip, set())[-1])

    def test_reverse_genomic_motif_and_noncanonical(self):
        genome = 'A' * 20 + 'CT' + 'C' * 36 + 'AC' + 'T' * 20
        read = 'A' * 20 + 'T' * 20
        details = MODULE.junction_details('c', genome, read, (20, 60, 20, 20, 20), {('c', 20, 60, '-')})
        self.assertEqual(details, ('CTAC', '-', True, True, True))
        genome = 'A' * 20 + 'AA' + 'C' * 36 + 'AA' + 'T' * 20
        self.assertFalse(MODULE.junction_details('c', genome, read, (20, 60, 20, 20, 20), set())[-1])

    def test_annotation_exon_precedence_and_half_open_introns(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'test.gff'
            path.write_text('c\t.\tgene\t11\t90\t.\t+\t.\tID=g\n'
                            'c\t.\texon\t11\t30\t.\t+\t.\tParent=t\n'
                            'c\t.\texon\t71\t90\t.\t+\t.\tParent=t\n')
            masks, _, introns, space = MODULE.read_annotation(path, {'c': 'A' * 100})
            self.assertEqual(introns, {('c', 30, 70, '+')})
            self.assertEqual(space['all_reference_bases'], {'intergenic': 20, 'gene_nonexonic': 40, 'exonic': 40})
            self.assertEqual(masks['c'][9:12], bytearray([0, 2, 2]))

    def test_sam_denominators_primary_flags_and_base_counts(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            genome = {'c': 'A' * 20 + 'GT' + 'C' * 36 + 'AG' + 'T' * 20}
            masks = {'c': bytearray([2] * 20 + [1] * 40 + [2] * 20)}
            read = 'A' * 20 + 'T' * 20
            rows = [f'r\t0\tc\t1\t60\t20M40N20M\t*\t0\t0\t{read}\t*',
                    f'r\t256\tc\t1\t60\t20M40N20M\t*\t0\t0\t{read}\t*',
                    f'r2\t0\tc\t1\t255\t20M40N20M\t*\t0\t0\t{read}\t*',
                    f'r3\t4\t*\t0\t0\t*\t*\t0\t0\t{read}\t*']
            sam = root / 'test.sam'
            sam.write_text('\n'.join(rows) + '\n')
            summary = MODULE.analyze_sam(sam, genome, masks, {('c', 20, 60, '+')}, root)
            self.assertEqual(summary['counts']['primary_records'], 3)
            self.assertEqual(summary['counts']['primary_mapped_records'], 2)
            self.assertEqual(summary['counts']['high_mapq_primary_mapped_records'], 1)
            self.assertEqual(summary['high_mapq_strong_spliced_fraction'], 1)
            self.assertEqual(summary['high_mapq_aligned_reference_bases']['exonic'], 40)
            self.assertEqual(summary['high_mapq_aligned_reference_bases']['gene_nonexonic'], 0)
            self.assertEqual(summary['distinct_strong_annotation_matched_junctions'], 1)


if __name__ == '__main__':
    unittest.main()
