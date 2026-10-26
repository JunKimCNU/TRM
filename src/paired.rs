//! Optional evidence for explicitly declared inward-facing (FR) read pairs.
//!
//! Counts describe families already detected in each read. Opposite orientation
//! is compatible with sampling one oriented array; it neither establishes a
//! telomere nor joins the mates into a longer observed tract. A mate without a
//! detected repeat is valid: the fragment can extend into subtelomeric sequence,
//! or that mate may not contain a detectable seed. These annotations must not
//! replace ordinary per-read counts or reject a motif automatically.

use crate::detect::Tract;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct PairStats {
    /// Input pairs with at least one detected tract of this family; counted once.
    /// This is not a count of deduplicated physical DNA molecules.
    pub supporting_pairs: u64,
    pub both_mates_repeat_pairs: u64,
    pub r1_only_repeat_pairs: u64,
    pub r2_only_repeat_pairs: u64,
    /// Both mates support this family, each with one unambiguous orientation.
    /// Orientations are relative to the shared canonical key, not C/G richness.
    pub opposite_strand_pairs: u64,
    pub same_strand_pairs: u64,
    /// Both mates support the family, but at least one has an unknown strand or
    /// a mixture of opposing tracts. Never credited as opposite or same strand.
    pub ambiguous_strand_pairs: u64,
    /// A full-repeat mate has at least one individual tract covering >=90%.
    /// These descriptive subsets overlap the orientation/support categories.
    pub both_full_repeat_pairs: u64,
    pub exactly_one_full_repeat_pairs: u64,
}

/// Conservative per-mate orientation: every detected tract must agree. Repeated
/// or overlapping hits of one orientation do not multiply the observation.
fn orientation(tracts: &[Tract]) -> Option<bool> {
    let strand = tracts.first()?.strand?;
    tracts
        .iter()
        .all(|tract| tract.strand == Some(strand))
        .then_some(strand)
}

fn full_repeat(tracts: &[Tract], length: usize) -> bool {
    length > 0
        && tracts
            .iter()
            .any(|tract| (tract.end - tract.start) as f64 >= length as f64 * 0.9)
}

impl PairStats {
    /// Call once per family per synchronized pair; use an empty slice for an
    /// unsupported mate. Tract coordinates have already been validated by the
    /// detector. No sequence rescanning, mapping or overlap inference occurs.
    pub fn observe(&mut self, r1: &[Tract], r1_len: usize, r2: &[Tract], r2_len: usize) {
        match (!r1.is_empty(), !r2.is_empty()) {
            (false, false) => return,
            (true, false) => self.r1_only_repeat_pairs += 1,
            (false, true) => self.r2_only_repeat_pairs += 1,
            (true, true) => {
                self.both_mates_repeat_pairs += 1;
                match (orientation(r1), orientation(r2)) {
                    (Some(a), Some(b)) if a != b => self.opposite_strand_pairs += 1,
                    (Some(_), Some(_)) => self.same_strand_pairs += 1,
                    _ => self.ambiguous_strand_pairs += 1,
                }
            }
        }
        self.supporting_pairs += 1;
        match (full_repeat(r1, r1_len), full_repeat(r2, r2_len)) {
            (true, true) => self.both_full_repeat_pairs += 1,
            (true, false) | (false, true) => self.exactly_one_full_repeat_pairs += 1,
            _ => {}
        }
    }

    pub fn merge(&mut self, other: Self) {
        self.supporting_pairs += other.supporting_pairs;
        self.both_mates_repeat_pairs += other.both_mates_repeat_pairs;
        self.r1_only_repeat_pairs += other.r1_only_repeat_pairs;
        self.r2_only_repeat_pairs += other.r2_only_repeat_pairs;
        self.opposite_strand_pairs += other.opposite_strand_pairs;
        self.same_strand_pairs += other.same_strand_pairs;
        self.ambiguous_strand_pairs += other.ambiguous_strand_pairs;
        self.both_full_repeat_pairs += other.both_full_repeat_pairs;
        self.exactly_one_full_repeat_pairs += other.exactly_one_full_repeat_pairs;
    }
}

/// SRA/ENA archive spot ID, without stripping any dotted part of the identity.
fn archive_spot(id: &str) -> Option<&str> {
    let (run, spot) = id.split_once('.')?;
    let digits = ["SRR", "ERR", "DRR"]
        .iter()
        .find_map(|prefix| run.strip_prefix(prefix))?;
    (!digits.is_empty()
        && digits.bytes().all(|byte| byte.is_ascii_digit())
        && !spot.is_empty()
        && spot.bytes().all(|byte| byte.is_ascii_digit()))
    .then_some(spot)
}

/// Extract first-token identity and optional standard mate label. Recognized
/// comments are Illumina `1:...` / `2:...`, or exactly `<spot>/1` / `<spot>/2`
/// after an archive `SRR|ERR|DRR<digits>.<same spot>` identity. Other comments
/// carry no inferred label. A slash suffix and comment label must agree.
fn identity(header: &str) -> Result<(&str, Option<u8>), String> {
    let mut words = header.split_whitespace();
    let token = words.next().ok_or("paired FASTQ header is empty")?;
    let token = token.strip_prefix('@').unwrap_or(token);
    let (id, suffix) = if let Some(id) = token.strip_suffix("/1") {
        (id, Some(1))
    } else if let Some(id) = token.strip_suffix("/2") {
        (id, Some(2))
    } else {
        (token, None)
    };
    if id.is_empty() {
        return Err("paired FASTQ read identity is empty".into());
    }
    let comment = words.next().and_then(|word| {
        if word.starts_with("1:") {
            Some(1)
        } else if word.starts_with("2:") {
            Some(2)
        } else if let Some(spot) = archive_spot(id) {
            let (comment_spot, mate) = word.split_once('/')?;
            if comment_spot == spot {
                match mate {
                    "1" => Some(1),
                    "2" => Some(2),
                    _ => None,
                }
            } else {
                None
            }
        } else {
            None
        }
    });
    if suffix.is_some() && comment.is_some() && suffix != comment {
        return Err(format!(
            "conflicting /1 or /2 suffix and comment mate label in header {header:?}"
        ));
    }
    Ok((id, suffix.or(comment)))
}

/// Validate an already synchronized R1/R2 pair. Requires exact first-token IDs
/// after stripping only a recognized `/1` or `/2` suffix. Dots or trailing
/// digits (including SRA spot IDs) are never heuristically removed. An unlabeled
/// mate is allowed when its normalized ID matches, but every explicit label
/// must match the file position. The caller rejects unequal stream lengths.
pub fn validate_pair_ids(r1_header: &str, r2_header: &str) -> Result<(), String> {
    let (r1, label1) = identity(r1_header)?;
    let (r2, label2) = identity(r2_header)?;
    if label1.is_some() && label1 != Some(1) {
        return Err(format!("R1 input has an R2 mate label: {r1_header:?}"));
    }
    if label2.is_some() && label2 != Some(2) {
        return Err(format!("R2 input has an R1 mate label: {r2_header:?}"));
    }
    if r1 != r2 {
        return Err(format!(
            "paired FASTQ IDs differ or are out of order: R1={r1:?}, R2={r2:?}"
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::{discover, Config};
    use crate::motif::{canonical, reverse_complement};

    fn tract(start: usize, end: usize, strand: Option<bool>) -> Tract {
        Tract {
            start,
            end,
            identity: 1.0,
            strand,
        }
    }

    fn assert_partitions(stats: &PairStats) {
        assert_eq!(
            stats.supporting_pairs,
            stats.both_mates_repeat_pairs + stats.r1_only_repeat_pairs + stats.r2_only_repeat_pairs
        );
        assert_eq!(
            stats.both_mates_repeat_pairs,
            stats.opposite_strand_pairs + stats.same_strand_pairs + stats.ambiguous_strand_pairs
        );
    }

    #[test]
    fn reverse_complement_mates_support_one_pair_without_double_counting() {
        let r1 = b"CCCTAA".repeat(25);
        let r2 = reverse_complement(&r1);
        let key = canonical(b"TTAGGG").0;
        let a = discover(&r1, &Config::default());
        let b = discover(&r2, &Config::default());
        let mut stats = PairStats::default();
        stats.observe(&a[&key], r1.len(), &b[&key], r2.len());
        assert_eq!(stats.supporting_pairs, 1);
        assert_eq!(stats.opposite_strand_pairs, 1);
        assert_eq!(stats.both_full_repeat_pairs, 1);
        assert_partitions(&stats);
    }

    #[test]
    fn missing_repeat_in_valid_mate_is_not_orientation_conflict() {
        let mut stats = PairStats::default();
        stats.observe(&[tract(0, 90, Some(true))], 150, &[], 150);
        stats.observe(&[], 150, &[tract(0, 150, Some(false))], 150);
        stats.observe(&[], 150, &[], 150);
        assert_eq!(stats.supporting_pairs, 2);
        assert_eq!(stats.r1_only_repeat_pairs, 1);
        assert_eq!(stats.r2_only_repeat_pairs, 1);
        assert_eq!(stats.both_mates_repeat_pairs, 0);
        assert_eq!(stats.same_strand_pairs, 0);
        assert_eq!(stats.ambiguous_strand_pairs, 0);
        assert_eq!(stats.exactly_one_full_repeat_pairs, 1);
        assert_partitions(&stats);
    }

    #[test]
    fn overlapping_same_orientation_hits_count_once() {
        let mut stats = PairStats::default();
        stats.observe(
            &[tract(0, 100, Some(true)), tract(50, 150, Some(true))],
            150,
            &[tract(0, 120, Some(false))],
            150,
        );
        assert_eq!(stats.supporting_pairs, 1);
        assert_eq!(stats.opposite_strand_pairs, 1);
        // A union covering a full read is not an individual full-repeat tract.
        assert_eq!(stats.both_full_repeat_pairs, 0);
        assert_eq!(stats.exactly_one_full_repeat_pairs, 0);
        assert_partitions(&stats);
    }

    #[test]
    fn same_strand_is_annotation_and_retains_support() {
        let mut stats = PairStats::default();
        stats.observe(
            &[tract(0, 60, Some(false))],
            150,
            &[tract(40, 150, Some(false))],
            150,
        );
        assert_eq!(stats.supporting_pairs, 1);
        assert_eq!(stats.both_mates_repeat_pairs, 1);
        assert_eq!(stats.same_strand_pairs, 1);
        assert_eq!(stats.opposite_strand_pairs, 0);
        assert_partitions(&stats);
    }

    #[test]
    fn opposing_tracts_and_palindromic_families_are_ambiguous() {
        let mut stats = PairStats::default();
        stats.observe(
            &[tract(0, 60, Some(true)), tract(90, 150, Some(false))],
            150,
            &[tract(0, 60, Some(false))],
            150,
        );
        stats.observe(&[tract(0, 150, None)], 150, &[tract(0, 150, None)], 150);
        stats.observe(
            &[tract(0, 60, Some(true)), tract(90, 150, None)],
            150,
            &[tract(0, 60, Some(false))],
            150,
        );
        assert_eq!(stats.ambiguous_strand_pairs, 3);
        assert_eq!(stats.opposite_strand_pairs, 0);
        assert_eq!(stats.same_strand_pairs, 0);
        assert_partitions(&stats);
    }

    #[test]
    fn full_repeat_fraction_uses_each_mate_length() {
        let mut stats = PairStats::default();
        stats.observe(
            &[tract(0, 135, Some(true))],
            150,
            &[tract(0, 89, Some(false))],
            100,
        );
        stats.observe(
            &[tract(0, 135, Some(true))],
            150,
            &[tract(0, 90, Some(false))],
            100,
        );
        assert_eq!(stats.exactly_one_full_repeat_pairs, 1);
        assert_eq!(stats.both_full_repeat_pairs, 1);
        assert_partitions(&stats);
    }

    #[test]
    fn merge_preserves_pair_partitions_and_serialization() {
        let mut a = PairStats::default();
        a.observe(&[tract(0, 24, Some(true))], 150, &[], 150);
        let mut b = PairStats::default();
        b.observe(
            &[tract(0, 150, Some(true))],
            150,
            &[tract(0, 150, Some(false))],
            150,
        );
        a.merge(b);
        assert_eq!(a.supporting_pairs, 2);
        assert_eq!(a.both_full_repeat_pairs, 1);
        assert_partitions(&a);
        let restored: PairStats =
            serde_json::from_value(serde_json::to_value(&a).unwrap()).unwrap();
        assert_eq!(restored, a);
    }

    #[test]
    fn accepts_standard_slash_illumina_and_archive_names() {
        for (r1, r2) in [
            ("read/1", "read/2"),
            ("@read/1", "@read/2"),
            (
                "run:lane:tile:x:y 1:N:0:ACGT",
                "run:lane:tile:x:y 2:N:0:ACGT",
            ),
            ("read/1 1:N:0:ACGT", "read/2 2:N:0:ACGT"),
            ("SRR21149832.5", "SRR21149832.5"),
            ("read length=150", "read length=100"),
        ] {
            assert!(validate_pair_ids(r1, r2).is_ok(), "{r1:?}, {r2:?}");
        }
    }

    #[test]
    fn one_unlabeled_mate_is_accepted_when_identity_matches() {
        for (r1, r2) in [
            ("read/1", "read"),
            ("read", "read/2"),
            ("read 1:N:0:X", "read"),
        ] {
            assert!(validate_pair_ids(r1, r2).is_ok());
        }
    }

    #[test]
    fn rejects_wrong_side_labels_even_when_normalized_ids_match() {
        for (r1, r2) in [
            ("read/1", "read/1"),
            ("read/2", "read/2"),
            ("read/2", "read/1"),
            ("read", "read 1:N:0:X"),
            ("read 2:N:0:X", "read"),
        ] {
            assert!(validate_pair_ids(r1, r2).is_err(), "{r1:?}, {r2:?}");
        }
    }

    #[test]
    fn rejects_conflicting_explicit_labels() {
        assert!(validate_pair_ids("read/1 2:N:0:X", "read/2").is_err());
        assert!(validate_pair_ids("read/1", "read/2 1:N:0:X").is_err());
    }

    #[test]
    fn validates_exact_archive_spot_comment_without_guessing_generic_comments() {
        for prefix in ["SRR", "ERR", "DRR"] {
            let id = format!("{prefix}21149832.17");
            assert!(validate_pair_ids(&format!("{id} 17/1"), &format!("{id} 17/2")).is_ok());
            assert!(validate_pair_ids(&format!("{id} 17/2"), &format!("{id} 17/1")).is_err());
            assert!(validate_pair_ids(&format!("{id} 17/1"), &format!("{id} 17/1")).is_err());
            assert!(validate_pair_ids(&format!("{id}/1 17/2"), &format!("{id}/2 17/2")).is_err());
        }
        // Unrecognized comments are not fabricated into mate metadata.
        assert!(validate_pair_ids("custom.17 17/2", "custom.17 17/1").is_ok());
        assert!(validate_pair_ids("SRR1.17 9/2", "SRR1.17 9/1").is_ok());
    }

    #[test]
    fn rejects_out_of_order_ids_without_guessing_away_digits() {
        for (r1, r2) in [
            ("read1/1", "read2/2"),
            ("SRR21149832.1", "SRR21149832.2"),
            ("", "read"),
            ("@/1", "@/2"),
        ] {
            assert!(validate_pair_ids(r1, r2).is_err(), "{r1:?}, {r2:?}");
        }
    }
}
