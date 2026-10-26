//! Bounded, descriptive annotation of legacy-eligible junction flank events.
//!
//! Neither candidate-period agreement nor family-template agreement establishes
//! a subtelomere, physical chromosome end, or phase-continuous repeat extension.
use crate::detect::Tract;
use crate::evidence::QuadrantCounts;
use crate::motif::{canonical, complement, is_dna, primitive_len};
use serde::{Deserialize, Serialize};

const THRESHOLD_PERCENT: usize = 85;
const MAX_PERIOD: usize = 64;

#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct EventCounts {
    pub events: u64,
    /// Raw read-side and strand quadrants, relative to the canonical motif.
    /// These are event counts, not deduplicated reads or independent molecules.
    pub quadrant_counts: QuadrantCounts,
}

impl EventCounts {
    fn observe(&mut self, left: bool, strand: usize) {
        self.events += 1;
        let mut flags = [false; 3];
        flags[strand] = true;
        self.quadrant_counts.merge(if left {
            QuadrantCounts::flags(flags, [false; 3])
        } else {
            QuadrantCounts::flags([false; 3], flags)
        });
    }

    fn merge(&mut self, other: Self) {
        self.events += other.events;
        self.quadrant_counts.merge(other.quadrant_counts);
    }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct EventClasses {
    pub insufficient_window: EventCounts,
    pub neither: EventCounts,
    pub lag_only: EventCounts,
    pub family_only: EventCounts,
    pub both: EventCounts,
}

impl EventClasses {
    fn merge(&mut self, other: Self) {
        self.insufficient_window.merge(other.insufficient_window);
        self.neither.merge(other.neither);
        self.lag_only.merge(other.lag_only);
        self.family_only.merge(other.family_only);
        self.both.merge(other.both);
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct FlankPeriodAnnotation {
    pub schema_version: u32,
    pub candidate_motif: String,
    pub primitive_period_bases: usize,
    pub required_window_bases: usize,
    pub threshold_percent: usize,
    /// Number of legacy-eligible tract-side events actually annotated.
    pub observed_events: u64,
    /// False if a merged or previously observed junction-bearing contribution
    /// was unmeasured. Empty unmeasured contributions do not invalidate this.
    pub measured_all_events: bool,
    pub classes: EventClasses,
}

impl FlankPeriodAnnotation {
    pub fn new(motif: &[u8]) -> Self {
        assert!(
            !motif.is_empty() && motif.len() <= MAX_PERIOD && is_dna(motif),
            "flank-period annotation requires 1..64 uppercase A/C/G/T motif bases"
        );
        let period = primitive_len(motif);
        Self {
            schema_version: 1,
            candidate_motif: canonical(&motif[..period]).0,
            primitive_period_bases: period,
            required_window_bases: 20.max(3 * period),
            threshold_percent: THRESHOLD_PERCENT,
            observed_events: 0,
            measured_all_events: true,
            classes: EventClasses::default(),
        }
    }

    /// Accept rotations, reverse complements, and nonprimitive spellings of
    /// the same family without allocating in the per-record observation path.
    pub(crate) fn assert_motif(&self, motif: &[u8]) {
        if motif == self.candidate_motif.as_bytes() {
            return;
        }
        assert!(
            !motif.is_empty() && motif.len() <= MAX_PERIOD && is_dna(motif),
            "invalid flank-period annotation motif"
        );
        let period = primitive_len(motif);
        assert_eq!(
            period, self.primitive_period_bases,
            "flank-period motif period mismatch"
        );
        let expected = self.candidate_motif.as_bytes();
        let matches = (0..period).any(|phase| {
            (0..period).all(|i| motif[(i + phase) % period] == expected[i])
                || (0..period)
                    .all(|i| complement(motif[period - 1 - ((i + phase) % period)]) == expected[i])
        });
        assert!(matches, "flank-period motif family mismatch");
    }

    pub(crate) fn observe(
        &mut self,
        sequence: &[u8],
        tract: &Tract,
        tracts: &[Tract],
        left: bool,
        strand: usize,
    ) {
        let flank = adjacent_window(sequence, tract, tracts, left, self.required_window_bases);
        self.observed_events += 1;
        let class = if flank.len() < self.required_window_bases {
            &mut self.classes.insufficient_window
        } else {
            let (lag, family) = agreement_tests(flank, self.candidate_motif.as_bytes());
            match (lag, family) {
                (false, false) => &mut self.classes.neither,
                (true, false) => &mut self.classes.lag_only,
                (false, true) => &mut self.classes.family_only,
                (true, true) => &mut self.classes.both,
            }
        };
        class.observe(left, strand);
    }

    pub(crate) fn merge(&mut self, other: Self) {
        // Incompatible measurements must fail explicitly rather than silently
        // combine thresholds, motif families, periods, or window dimensions.
        assert_eq!(
            self.schema_version, other.schema_version,
            "flank-period schema mismatch"
        );
        assert_eq!(
            self.candidate_motif, other.candidate_motif,
            "flank-period motif family mismatch"
        );
        assert_eq!(
            self.primitive_period_bases, other.primitive_period_bases,
            "flank-period period mismatch"
        );
        assert_eq!(
            self.required_window_bases, other.required_window_bases,
            "flank-period window mismatch"
        );
        assert_eq!(
            self.threshold_percent, other.threshold_percent,
            "flank-period threshold mismatch"
        );
        self.observed_events += other.observed_events;
        self.measured_all_events &= other.measured_all_events;
        self.classes.merge(other.classes);
    }
}

/// Retain physical read order on either side. Stop at the first ambiguous
/// base or same-family detected tract; never skip an obstruction to fill a
/// window. At most required bases are inspected, and no sequence is retained.
fn adjacent_window<'a>(
    sequence: &'a [u8],
    tract: &Tract,
    tracts: &[Tract],
    left: bool,
    required: usize,
) -> &'a [u8] {
    if left {
        let end = tracts
            .iter()
            .filter(|other| other.start >= tract.end)
            .map(|other| other.start)
            .min()
            .unwrap_or(sequence.len())
            .min(sequence.len())
            .min(tract.end.saturating_add(required));
        let flank = &sequence[tract.end..end];
        let length = flank
            .iter()
            .take_while(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            .count();
        &flank[..length]
    } else {
        let start = tracts
            .iter()
            .filter(|other| other.end <= tract.start)
            .map(|other| other.end)
            .max()
            .unwrap_or(0)
            .max(tract.start.saturating_sub(required));
        let flank = &sequence[start..tract.start];
        let offset = flank
            .iter()
            .rposition(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
            .map_or(0, |position| position + 1);
        &flank[offset..]
    }
}

/// Candidate-period autocorrelation and candidate-family template identity
/// are independent statistics. The latter tests all phases in both strands;
/// it does not constrain phase or strand to the neighboring detected tract.
fn agreement_tests(flank: &[u8], motif: &[u8]) -> (bool, bool) {
    let period = motif.len();
    debug_assert!(period > 0 && period <= MAX_PERIOD && flank.len() >= 3 * period);
    let agreements = (period..flank.len())
        .filter(|&i| flank[i] == flank[i - period])
        .count();
    let lag = agreements * 100 >= (flank.len() - period) * THRESHOLD_PERCENT;
    // Fixed stack storage, bounded by period 64 and window 192. Accumulating
    // by phase costs O(window + period^2), with no per-event heap allocation.
    let mut counts = [[0_u16; 4]; MAX_PERIOD];
    for (i, &base) in flank.iter().enumerate() {
        counts[i % period][base_index(base)] += 1;
    }
    let family = (0..period).any(|phase| {
        let forward: usize = (0..period)
            .map(|i| usize::from(counts[i][base_index(motif[(i + phase) % period])]))
            .sum();
        let reverse: usize = (0..period)
            .map(|i| {
                usize::from(
                    counts[i][base_index(complement(motif[period - 1 - ((i + phase) % period)]))],
                )
            })
            .sum();
        forward.max(reverse) * 100 >= flank.len() * THRESHOLD_PERCENT
    });
    (lag, family)
}

fn base_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => unreachable!("annotation windows contain only uppercase A/C/G/T"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::motif::reverse_complement;

    fn tract(start: usize, end: usize) -> Tract {
        Tract {
            start,
            end,
            identity: 1.0,
            strand: Some(true),
        }
    }

    fn primitive_motif(period: usize) -> Vec<u8> {
        let mut motif = vec![b'A'; period];
        motif[period - 1] = b'C';
        motif
    }

    #[test]
    fn every_supported_period_has_bounded_rotation_and_strand_invariant_tests() {
        for period in 1..=64 {
            let motif = primitive_motif(period);
            let annotation = FlankPeriodAnnotation::new(&motif);
            assert_eq!(annotation.primitive_period_bases, period);
            assert_eq!(annotation.required_window_bases, 20.max(3 * period));
            assert!(annotation.required_window_bases <= 192);
            for phase in 0..period {
                let rotated: Vec<_> = (0..period).map(|i| motif[(i + phase) % period]).collect();
                for spelling in [&rotated, &reverse_complement(&rotated)] {
                    annotation.assert_motif(spelling);
                    let flank: Vec<_> = (0..annotation.required_window_bases)
                        .map(|i| spelling[i % period])
                        .collect();
                    assert_eq!(
                        agreement_tests(&flank, annotation.candidate_motif.as_bytes()),
                        (true, true)
                    );
                }
            }
        }
        assert_eq!(
            FlankPeriodAnnotation::new(b"TTAGGGGTTAGGGG").primitive_period_bases,
            7
        );
    }

    #[test]
    fn lag_and_family_tests_form_four_distinct_classes() {
        let motif = b"AACCCCT";
        let both = motif.repeat(3);
        let lag_only = b"ACGTCGT".repeat(3);
        let mut family_only = both.clone();
        family_only[7] = b'T';
        family_only[8] = b'G';
        let neither = b"AGTCGATCGGACTACGATCAG";
        assert_eq!(agreement_tests(&both, motif), (true, true));
        assert_eq!(agreement_tests(&lag_only, motif), (true, false));
        assert_eq!(agreement_tests(&family_only, motif), (false, true));
        assert_eq!(agreement_tests(neither, motif), (false, false));
    }

    #[test]
    fn integer_threshold_includes_eighty_five_percent_exactly() {
        let mut flank = vec![b'A'; 20];
        flank[..3].fill(b'C');
        assert!(agreement_tests(&flank, b"A").1); // 17 / 20
        flank[3] = b'C';
        assert!(!agreement_tests(&flank, b"A").1); // 16 / 20
        let motif = b"AACCGGTTAC";
        let mut flank = motif.repeat(3);
        for base in flank.iter_mut().take(3) {
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
        assert!(agreement_tests(&flank, motif).0); // 17 / 20 lag comparisons
        flank[3] = b'A';
        assert!(!agreement_tests(&flank, motif).0); // 16 / 20
    }

    #[test]
    fn bounded_phase_histogram_matches_direct_template_oracle() {
        let mut state = 918_273_u64;
        let mut base = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        };
        for period in [1, 2, 6, 7, 32, 64] {
            for case in 0..16 {
                let motif: Vec<_> = (0..period).map(|_| base()).collect();
                let length = 20.max(3 * period);
                let mut flank: Vec<_> = (0..length).map(|i| motif[(i + case) % period]).collect();
                for _ in 0..case {
                    let i = usize::from(base()) % length;
                    flank[i] = base();
                }
                let reverse = reverse_complement(&motif);
                let best = [&motif, &reverse]
                    .into_iter()
                    .flat_map(|unit| {
                        (0..period).map(|phase| {
                            flank
                                .iter()
                                .enumerate()
                                .filter(|&(i, &b)| b == unit[(i + phase) % period])
                                .count()
                        })
                    })
                    .max()
                    .unwrap();
                let lag = (period..length)
                    .filter(|&i| flank[i] == flank[i - period])
                    .count();
                assert_eq!(
                    agreement_tests(&flank, &motif),
                    (
                        lag * 100 >= (length - period) * 85,
                        best * 100 >= length * 85
                    )
                );
            }
        }
    }

    #[test]
    fn adjacent_window_never_bridges_ambiguity_or_detected_tracts() {
        let sequence = b"AAAAACGTACGTNACGTACGTAAAA";
        let left = tract(0, 4);
        assert_eq!(
            adjacent_window(sequence, &left, std::slice::from_ref(&left), true, 21),
            b"ACGTACGT"
        );
        let right = tract(20, 24);
        assert_eq!(
            adjacent_window(sequence, &right, std::slice::from_ref(&right), false, 21),
            b"ACGTACG"
        );
        let barrier = tract(7, 9);
        assert_eq!(
            adjacent_window(sequence, &left, &[left.clone(), barrier.clone()], true, 21),
            b"ACG"
        );
        let unambiguous = b"AAAAACGTACGTACGTACGTAAAA";
        assert_eq!(
            adjacent_window(unambiguous, &right, &[right.clone(), barrier], false, 21),
            &unambiguous[9..20]
        );
        assert_eq!(
            adjacent_window(unambiguous, &left, std::slice::from_ref(&left), true, 3),
            &unambiguous[4..7]
        );
    }

    #[test]
    fn class_event_partition_and_raw_quadrants_merge_exactly() {
        let motif = b"AACCCCT";
        let base_tract = tract(0, 28);
        let mut expected = FlankPeriodAnnotation::new(motif);
        let mut groups = Vec::new();
        let mut family_only = motif.repeat(3);
        family_only[7] = b'T';
        family_only[8] = b'G';
        for (i, flank) in [
            motif.repeat(3),
            b"ACGTCGT".repeat(3),
            family_only,
            b"AGTCGATCGGACTACGATCAG".to_vec(),
            motif.repeat(2),
        ]
        .into_iter()
        .enumerate()
        {
            let sequence = [motif.repeat(4), flank].concat();
            let mut part = FlankPeriodAnnotation::new(motif);
            part.observe(
                &sequence,
                &base_tract,
                std::slice::from_ref(&base_tract),
                true,
                i % 3,
            );
            expected.observe(
                &sequence,
                &base_tract,
                std::slice::from_ref(&base_tract),
                true,
                i % 3,
            );
            groups.push(part);
        }
        let mut forward = FlankPeriodAnnotation::new(motif);
        let mut reverse = FlankPeriodAnnotation::new(motif);
        for group in groups.clone() {
            forward.merge(group);
        }
        for group in groups.into_iter().rev() {
            reverse.merge(group);
        }
        assert_eq!(forward, expected);
        assert_eq!(reverse, expected);
        assert_eq!(expected.observed_events, 5);
        for class in [
            &expected.classes.insufficient_window,
            &expected.classes.neither,
            &expected.classes.lag_only,
            &expected.classes.family_only,
            &expected.classes.both,
        ] {
            assert_eq!(class.events, 1);
            let q = &class.quadrant_counts;
            assert_eq!(
                q.left_forward + q.left_reverse + q.right_forward + q.right_reverse + q.ambiguous,
                1
            );
        }
    }

    #[test]
    #[should_panic(expected = "motif family mismatch")]
    fn merging_different_families_fails_explicitly() {
        FlankPeriodAnnotation::new(b"TTAGGGG").merge(FlankPeriodAnnotation::new(b"TTAGGGC"));
    }

    #[test]
    #[should_panic(expected = "window mismatch")]
    fn merging_different_window_dimensions_fails_explicitly() {
        let mut changed = FlankPeriodAnnotation::new(b"TTAGGGG");
        changed.required_window_bases += 1;
        FlankPeriodAnnotation::new(b"TTAGGGG").merge(changed);
    }
}
