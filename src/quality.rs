//! Opt-in, explicitly declared Phred+33 diagnostics. Never a detection filter.
use crate::detect::Tract;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct Histogram {
    /// Index is the explicitly declared Phred+33 score, 0 through 93.
    pub phred33_counts: Vec<u64>,
    pub bases: u64,
    pub below_q20_bases: u64,
    pub below_q30_bases: u64,
}

impl Default for Histogram {
    fn default() -> Self {
        Self {
            phred33_counts: vec![0; 94],
            bases: 0,
            below_q20_bases: 0,
            below_q30_bases: 0,
        }
    }
}

impl Histogram {
    fn observe(&mut self, quality: &[u8]) {
        for &byte in quality {
            // The input parser has validated the complete quality alphabet.
            let score = usize::from(byte - 33);
            self.phred33_counts[score] += 1;
            self.bases += 1;
            self.below_q20_bases += u64::from(score < 20);
            self.below_q30_bases += u64::from(score < 30);
        }
    }

    pub fn merge(&mut self, other: Self) {
        for (left, right) in self.phred33_counts.iter_mut().zip(other.phred33_counts) {
            *left += right;
        }
        self.bases += other.bases;
        self.below_q20_bases += other.below_q20_bases;
        self.below_q30_bases += other.below_q30_bases;
    }
}

#[derive(Clone, Default, Debug, Serialize)]
pub struct RecordQuality {
    pub observed_records: u64,
    pub quality_records: u64,
    pub fasta_records_without_quality: u64,
    pub fastq_records_above_retention_limit: u64,
    pub whole_read: Histogram,
}

impl RecordQuality {
    pub fn observe(&mut self, quality: Option<&[u8]>, is_fastq: bool) {
        self.observed_records += 1;
        match quality {
            Some(quality) => {
                self.quality_records += 1;
                self.whole_read.observe(quality);
            }
            None if is_fastq => self.fastq_records_above_retention_limit += 1,
            None => self.fasta_records_without_quality += 1,
        }
    }

    pub fn merge(&mut self, other: Self) {
        self.observed_records += other.observed_records;
        self.quality_records += other.quality_records;
        self.fasta_records_without_quality += other.fasta_records_without_quality;
        self.fastq_records_above_retention_limit += other.fastq_records_above_retention_limit;
        self.whole_read.merge(other.whole_read);
    }
}

#[derive(Clone, Default, Debug, Serialize)]
pub struct CandidateQuality {
    pub supporting_records: RecordQuality,
    pub union_tract_bases: u64,
    pub union_tract_bases_without_quality: u64,
    pub union_tract: Histogram,
}

impl CandidateQuality {
    pub fn observe(&mut self, quality: Option<&[u8]>, is_fastq: bool, tracts: &[Tract]) {
        self.supporting_records.observe(quality, is_fastq);
        // Detector tracts are sorted by start. Match the ordinary repeat-base
        // union so overlapping seeds, strands or tracts never count twice.
        let mut covered_end = 0;
        for tract in tracts {
            let start = covered_end.max(tract.start);
            if start < tract.end {
                let bases = (tract.end - start) as u64;
                self.union_tract_bases += bases;
                if let Some(quality) = quality {
                    self.union_tract.observe(&quality[start..tract.end]);
                } else {
                    self.union_tract_bases_without_quality += bases;
                }
            }
            covered_end = covered_end.max(tract.end);
        }
    }

    pub fn merge(&mut self, other: Self) {
        self.supporting_records.merge(other.supporting_records);
        self.union_tract_bases += other.union_tract_bases;
        self.union_tract_bases_without_quality += other.union_tract_bases_without_quality;
        self.union_tract.merge(other.union_tract);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thresholds_are_literal_phred33_and_all_legal_scores_are_counted() {
        let mut h = Histogram::default();
        h.observe(&(33..=126).collect::<Vec<_>>());
        assert_eq!(h.bases, 94);
        assert_eq!(h.below_q20_bases, 20);
        assert_eq!(h.below_q30_bases, 30);
        assert!(h.phred33_counts.iter().all(|&n| n == 1));
    }

    #[test]
    fn union_and_unavailable_quality_have_separate_denominators() {
        let t = |start, end| Tract {
            start,
            end,
            identity: 1.0,
            strand: Some(true),
        };
        let tracts = vec![t(1, 5), t(3, 7), t(9, 10)];
        let mut a = CandidateQuality::default();
        a.observe(Some(b"IIIIIIIIII"), true, &tracts);
        let mut b = CandidateQuality::default();
        b.observe(None, false, &tracts);
        b.observe(None, true, &tracts);
        a.merge(b);
        assert_eq!(a.union_tract_bases, 21);
        assert_eq!(a.union_tract.bases, 7);
        assert_eq!(a.union_tract_bases_without_quality, 14);
        assert_eq!(a.supporting_records.observed_records, 3);
        assert_eq!(a.supporting_records.quality_records, 1);
        assert_eq!(a.supporting_records.fasta_records_without_quality, 1);
        assert_eq!(a.supporting_records.fastq_records_above_retention_limit, 1);
        assert_eq!(a.supporting_records.whole_read.bases, 10);
    }
}
