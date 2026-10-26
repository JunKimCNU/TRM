//! Exact recognition of explicitly supplied, phase-consistent unit grammars.
//!
//! A grammar is `(unit_1|unit_2|...)+`, with units used literally: individual
//! rotations are NOT added. Reverse-complement units form a separate grammar.
//! Thus a reported core is a contiguous, complete-unit parse, never a join of
//! nearby homogeneous tracts. No homogeneous seed or mismatch bridge is used.
//!
//! The endpoint DP retains the earliest reachable covered start. Later starts at
//! the same endpoint are dominated because every subsequent unit can extend
//! either parse. Tied earliest parses retain marginal minimum/maximum counts;
//! these bounds need not be simultaneously attainable in a single parse.
//! `unit_counts` therefore means guaranteed counts, not an invented best parse.
//! Two components must EACH meet the configured guaranteed count to call mixed.
//!
//! Working memory is O(max_unit_length * unit_count), excluding returned hits;
//! matching takes O(record_length * sum(unit_lengths)) per definition/strand.
//! Maximal containment filtering is linear and never combines separated cores.

use crate::motif::{canonical, is_dna, primitive_len, reverse_complement};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct MixedDefinition {
    pub id: String,
    /// Literal G-strand phases. `MixedScanner::new` sorts units lexically.
    pub units: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct MixedConfig {
    pub min_span: usize,
    pub min_component_copies: usize,
    pub allow_edge_partials: bool,
}

impl Default for MixedConfig {
    fn default() -> Self {
        Self {
            min_span: 24,
            min_component_copies: 2,
            allow_edge_partials: true,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct MixedTract {
    pub definition_id: String,
    /// Zero-based, half-open exact coverage, including compatible edge partials.
    pub start: usize,
    pub end: usize,
    /// Strict complete-unit core of a deterministic valid representative parse;
    /// never includes a partial unit or a gap. When core_bounds_ambiguous is
    /// true, other equally covering parses have different complete-core bounds.
    pub core_start: usize,
    pub core_end: usize,
    /// Some(true): supplied G strand; Some(false): its reverse complement.
    /// None: overlapping evidence admits both orientations.
    pub strand: Option<bool>,
    /// Guaranteed counts, ordered as `scanner.definitions[..].units`.
    pub unit_counts: Vec<usize>,
    /// Marginal upper bounds, same ordering. Bounds are not a joint parse.
    pub unit_counts_max: Vec<usize>,
    pub segmentation_ambiguous: bool,
    pub core_bounds_ambiguous: bool,
    pub orientation_ambiguous: bool,
    pub exact_covered_bases: usize,
    pub partial_edge_bases: usize,
    /// Number of components guaranteed to occur at least once.
    pub distinct_units: usize,
    pub is_mixed: bool,
    /// mixed / pure / mixed_below_component_threshold / composition_ambiguous.
    pub composition: String,
}

#[derive(Clone, Debug)]
pub struct MixedScanner {
    pub definitions: Vec<MixedDefinition>,
    pub config: MixedConfig,
    grammars: Vec<[Vec<Vec<u8>>; 2]>,
}

#[derive(Clone, Debug)]
struct State {
    start: usize,
    core_start: usize,
    end: usize,
    minimum: Vec<usize>,
    maximum: Vec<usize>,
    ambiguous: bool,
    core_bounds_ambiguous: bool,
    extended: bool,
}

impl MixedScanner {
    pub fn new(mut definitions: Vec<MixedDefinition>, config: MixedConfig) -> Result<Self, String> {
        if definitions.is_empty() || definitions.len() > 32 {
            return Err("mixed definitions require 1..=32 groups".into());
        }
        if config.min_span == 0 || config.min_component_copies == 0 {
            return Err("mixed minimum span and component copies must be positive".into());
        }
        let mut ids = HashSet::new();
        let mut grammars = Vec::new();
        for definition in &mut definitions {
            if definition.id.trim().is_empty() || !ids.insert(definition.id.clone()) {
                return Err("mixed definition IDs must be nonempty and unique".into());
            }
            if definition.units.len() < 2 || definition.units.len() > 8 {
                return Err(format!("{}: require 2..=8 units", definition.id));
            }
            for unit in &mut definition.units {
                *unit = unit.to_ascii_uppercase();
            }
            definition.units.sort();
            let mut families = HashSet::new();
            for unit in &definition.units {
                let bytes = unit.as_bytes();
                if bytes.is_empty() || bytes.len() > 64 || !is_dna(bytes) {
                    return Err(format!(
                        "{}: units must be 1..=64 ACGT bases",
                        definition.id
                    ));
                }
                if primitive_len(bytes) != bytes.len() {
                    return Err(format!("{}: nonprimitive unit {unit}", definition.id));
                }
                if !families.insert(canonical(bytes).0) {
                    return Err(format!(
                        "{}: redundant rotation/RC-equivalent unit {unit}",
                        definition.id
                    ));
                }
            }
            // A literal concatenation of other supplied units would give an
            // artificial extra component rather than an independent unit.
            for (index, unit) in definition.units.iter().enumerate() {
                let mut reachable = vec![false; unit.len() + 1];
                reachable[0] = true;
                for start in 0..unit.len() {
                    if !reachable[start] {
                        continue;
                    }
                    for (other_index, other) in definition.units.iter().enumerate() {
                        if index != other_index
                            && unit.as_bytes()[start..].starts_with(other.as_bytes())
                        {
                            reachable[start + other.len()] = true;
                        }
                    }
                }
                if reachable[unit.len()] {
                    return Err(format!(
                        "{}: composite unit {unit} is a concatenation of other units",
                        definition.id
                    ));
                }
            }
            let forward: Vec<Vec<u8>> = definition
                .units
                .iter()
                .map(|s| s.as_bytes().to_vec())
                .collect();
            let reverse = forward.iter().map(|s| reverse_complement(s)).collect();
            grammars.push([forward, reverse]);
        }
        Ok(Self {
            definitions,
            config,
            grammars,
        })
    }

    pub fn scan(&self, seq: &[u8]) -> Vec<MixedTract> {
        let mut hits = Vec::new();
        for (definition, grammar) in self.definitions.iter().zip(&self.grammars) {
            let forward = self.scan_orientation(seq, definition, &grammar[0], true);
            let reverse = self.scan_orientation(seq, definition, &grammar[1], false);
            hits.extend(merge_orientations(
                forward,
                reverse,
                self.config.min_component_copies,
            ));
        }
        hits
    }

    fn scan_orientation(
        &self,
        seq: &[u8],
        definition: &MixedDefinition,
        units: &[Vec<u8>],
        strand: bool,
    ) -> Vec<MixedTract> {
        // Lossless negative prescreen: every reported core contains at least
        // one complete literal unit, including an alternating mixed array.
        // Avoid allocating/advancing the DP ring for the overwhelmingly common
        // WGS records that contain none. This is NOT a tandem-repeat seed.
        if !units.iter().any(|unit| {
            seq.len() >= unit.len()
                && seq
                    .windows(unit.len())
                    .any(|window| window == unit.as_slice())
        }) {
            return Vec::new();
        }
        let max_length = units.iter().map(Vec::len).max().unwrap();
        let ring_length = max_length + 1;
        let mut ring: Vec<Option<State>> = vec![None; ring_length];
        let mut hits = Vec::new();
        let mut left_partial = vec![false; max_length];
        if self.config.allow_edge_partials {
            for length in 1..max_length.min(seq.len() + 1) {
                left_partial[length] = units
                    .iter()
                    .any(|unit| length < unit.len() && unit.ends_with(&seq[..length]));
            }
        }
        for end in 1..=seq.len() {
            let slot = end % ring_length;
            if let Some(old) = ring[slot].take() {
                self.finish_state(old, seq, definition, units, strand, &mut hits);
            }
            let mut best: Option<State> = None;
            for (index, unit) in units.iter().enumerate() {
                if end < unit.len() || seq[end - unit.len()..end] != unit[..] {
                    continue;
                }
                let predecessor = &mut ring[(end - unit.len()) % ring_length];
                if let Some(previous) = predecessor {
                    debug_assert_eq!(previous.end, end - unit.len());
                    previous.extended = true;
                    let mut next = State {
                        start: previous.start,
                        core_start: previous.core_start,
                        end,
                        minimum: previous.minimum.clone(),
                        maximum: previous.maximum.clone(),
                        ambiguous: previous.ambiguous,
                        core_bounds_ambiguous: previous.core_bounds_ambiguous,
                        extended: false,
                    };
                    next.minimum[index] += 1;
                    next.maximum[index] += 1;
                    combine_state(&mut best, next);
                }
                // Starting a new complete core can improve covered start when
                // a compatible clipped unit reaches the left record edge.
                let core_start = end - unit.len();
                let start = if core_start < left_partial.len() && left_partial[core_start] {
                    0
                } else {
                    core_start
                };
                let mut seed = State {
                    start,
                    core_start,
                    end,
                    minimum: vec![0; units.len()],
                    maximum: vec![0; units.len()],
                    ambiguous: false,
                    core_bounds_ambiguous: false,
                    extended: false,
                };
                seed.minimum[index] = 1;
                seed.maximum[index] = 1;
                combine_state(&mut best, seed);
            }
            ring[slot] = best;
        }
        // Drain in endpoint order. Earlier states were already retired.
        for end in seq.len().saturating_sub(max_length).max(1)..=seq.len() {
            if let Some(state) = ring[end % ring_length].take() {
                self.finish_state(state, seq, definition, units, strand, &mut hits);
            }
        }
        hits
    }

    fn finish_state(
        &self,
        state: State,
        seq: &[u8],
        definition: &MixedDefinition,
        units: &[Vec<u8>],
        strand: bool,
        hits: &mut Vec<MixedTract>,
    ) {
        let start = state.start;
        let mut end = state.end;
        let mut right_partial = false;
        if self.config.allow_edge_partials {
            let remainder = seq.len() - end;
            if remainder > 0
                && units
                    .iter()
                    .any(|unit| remainder < unit.len() && unit.starts_with(&seq[end..]))
            {
                end = seq.len();
                right_partial = true;
            }
        }
        if state.extended && !right_partial {
            return;
        }
        if end - start < self.config.min_span {
            return;
        }
        let mut hit = MixedTract {
            definition_id: definition.id.clone(),
            start,
            end,
            core_start: state.core_start,
            core_end: state.end,
            strand: Some(strand),
            unit_counts: state.minimum,
            unit_counts_max: state.maximum,
            segmentation_ambiguous: state.ambiguous,
            core_bounds_ambiguous: state.core_bounds_ambiguous,
            orientation_ambiguous: false,
            exact_covered_bases: end - start,
            partial_edge_bases: state.core_start - start + end - state.end,
            distinct_units: 0,
            is_mixed: false,
            composition: String::new(),
        };
        classify(&mut hit, self.config.min_component_copies);
        // Endpoints arrive in increasing order. Drop only contained tracts of
        // this SAME grammar/orientation, never fill gaps between two tracts.
        if let Some(last) = hits.last_mut() {
            if hit.start == last.start && hit.end == last.end {
                merge_equal_span(last, hit, self.config.min_component_copies);
                return;
            }
        }
        while hits
            .last()
            .is_some_and(|last| hit.start <= last.start && hit.end >= last.end)
        {
            hits.pop();
        }
        if hits
            .last()
            .is_some_and(|last| last.start <= hit.start && last.end >= hit.end)
        {
            return;
        }
        hits.push(hit);
    }
}

fn combine_state(best: &mut Option<State>, next: State) {
    match best {
        None => *best = Some(next),
        Some(current) if next.start < current.start => *best = Some(next),
        Some(current) if next.start == current.start => {
            current.ambiguous = true;
            current.core_bounds_ambiguous |=
                next.core_bounds_ambiguous || current.core_start != next.core_start;
            current.core_start = current.core_start.min(next.core_start);
            for i in 0..current.minimum.len() {
                current.minimum[i] = current.minimum[i].min(next.minimum[i]);
                current.maximum[i] = current.maximum[i].max(next.maximum[i]);
            }
        }
        _ => {}
    }
}

fn merge_equal_span(current: &mut MixedTract, next: MixedTract, minimum_copies: usize) {
    current.segmentation_ambiguous |= next.segmentation_ambiguous
        || current.unit_counts != next.unit_counts
        || (current.core_start, current.core_end) != (next.core_start, next.core_end);
    current.core_bounds_ambiguous |= next.core_bounds_ambiguous
        || (current.core_start, current.core_end) != (next.core_start, next.core_end);
    if (next.core_start, std::cmp::Reverse(next.core_end))
        < (current.core_start, std::cmp::Reverse(current.core_end))
    {
        current.core_start = next.core_start;
        current.core_end = next.core_end;
        current.partial_edge_bases = next.partial_edge_bases;
    }
    if current.strand != next.strand {
        current.strand = None;
        current.orientation_ambiguous = true;
    }
    current.orientation_ambiguous |= next.orientation_ambiguous;
    for i in 0..current.unit_counts.len() {
        current.unit_counts[i] = current.unit_counts[i].min(next.unit_counts[i]);
        current.unit_counts_max[i] = current.unit_counts_max[i].max(next.unit_counts_max[i]);
    }
    classify(current, minimum_copies);
}

fn classify(hit: &mut MixedTract, minimum_copies: usize) {
    hit.distinct_units = hit.unit_counts.iter().filter(|&&n| n > 0).count();
    hit.is_mixed = hit
        .unit_counts
        .iter()
        .filter(|&&n| n >= minimum_copies)
        .count()
        >= 2;
    let possible = hit.unit_counts_max.iter().filter(|&&n| n > 0).count();
    hit.composition = if hit.is_mixed {
        "mixed"
    } else if possible == 1 {
        "pure"
    } else if hit.unit_counts != hit.unit_counts_max {
        "composition_ambiguous"
    } else {
        "mixed_below_component_threshold"
    }
    .into();
}

fn merge_orientations(
    forward: Vec<MixedTract>,
    reverse: Vec<MixedTract>,
    minimum_copies: usize,
) -> Vec<MixedTract> {
    let mut forward = forward.into_iter().peekable();
    let mut reverse = reverse.into_iter().peekable();
    let mut merged: Vec<MixedTract> = Vec::new();
    while forward.peek().is_some() || reverse.peek().is_some() {
        let choose_forward = match (forward.peek(), reverse.peek()) {
            (Some(f), Some(r)) => (f.start, f.end) <= (r.start, r.end),
            (Some(_), None) => true,
            _ => false,
        };
        let next = if choose_forward {
            forward.next().unwrap()
        } else {
            reverse.next().unwrap()
        };
        if let Some(last) = merged.last_mut() {
            if (last.start, last.end) == (next.start, next.end) {
                merge_equal_span(last, next, minimum_copies);
                continue;
            }
        }
        merged.push(next);
    }
    // Conservatively flag every tract in a connected overlap component that
    // contains both orientations. Preserve distinct cores instead of choosing
    // a direction by length or erasing the opposing evidence.
    let mut begin = 0;
    while begin < merged.len() {
        let mut end_index = begin + 1;
        let mut end_base = merged[begin].end;
        let mut has_forward = merged[begin].strand != Some(false);
        let mut has_reverse = merged[begin].strand != Some(true);
        while end_index < merged.len() && merged[end_index].start < end_base {
            end_base = end_base.max(merged[end_index].end);
            has_forward |= merged[end_index].strand != Some(false);
            has_reverse |= merged[end_index].strand != Some(true);
            end_index += 1;
        }
        if has_forward && has_reverse {
            for hit in &mut merged[begin..end_index] {
                hit.strand = None;
                hit.orientation_ambiguous = true;
            }
        }
        begin = end_index;
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scanner(units: &[&str], config: MixedConfig) -> MixedScanner {
        MixedScanner::new(
            vec![MixedDefinition {
                id: "test".into(),
                units: units.iter().map(|s| (*s).into()).collect(),
            }],
            config,
        )
        .unwrap()
    }

    fn standard() -> MixedScanner {
        scanner(&["TTGGGTTAC", "ACTTGGGTTAC"], MixedConfig::default())
    }

    #[test]
    fn alternating_units_need_no_homogeneous_seed() {
        let seq = "TTGGGTTACACTTGGGTTAC".repeat(4);
        let scanner = standard();
        let hits = scanner.scan(seq.as_bytes());
        assert_eq!(hits.len(), 1);
        let hit = &hits[0];
        assert_eq!((hit.start, hit.end), (0, seq.len()));
        assert_eq!((hit.core_start, hit.core_end), (0, seq.len()));
        // A leading 9mer can also be the clipped suffix of an 11mer.
        assert_eq!(hit.unit_counts, [4, 3]);
        assert_eq!(hit.unit_counts_max, [4, 4]);
        assert!(hit.core_bounds_ambiguous);
        assert!(hit.is_mixed);
        assert_eq!(hit.strand, Some(true));
    }

    #[test]
    fn reverse_complement_has_opposite_direction() {
        let seq = "TTGGGTTACACTTGGGTTAC".repeat(3);
        let hits = standard().scan(&reverse_complement(seq.as_bytes()));
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].unit_counts, [3, 2]);
        assert_eq!(hits[0].unit_counts_max, [3, 3]);
        assert_eq!(hits[0].strand, Some(false));
        assert!(hits[0].is_mixed);
    }

    #[test]
    fn clipped_read_partials_do_not_add_complete_copies() {
        let seq = "TTGGGTTACACTTGGGTTAC".repeat(5);
        let clipped = &seq.as_bytes()[3..seq.len() - 4];
        let hits = standard().scan(clipped);
        assert_eq!(hits.len(), 1);
        assert_eq!((hits[0].start, hits[0].end), (0, clipped.len()));
        assert_eq!(hits[0].partial_edge_bases, 13);
        assert_eq!(hits[0].unit_counts, [4, 4]);
        assert_eq!(hits[0].core_start, 6);
        assert_eq!(hits[0].core_end, clipped.len() - 7);
        let without = scanner(
            &["TTGGGTTAC", "ACTTGGGTTAC"],
            MixedConfig {
                allow_edge_partials: false,
                ..MixedConfig::default()
            },
        )
        .scan(clipped);
        assert_eq!(without[0].start, 6);
        assert_eq!(without[0].end, clipped.len() - 7);
        assert_eq!(without[0].partial_edge_bases, 0);
    }

    #[test]
    fn internal_partial_bases_are_not_credited() {
        let core = "TTGGGTTACACTTGGGTTAC".repeat(3);
        let seq = format!("NGGTTACN{core}NACTTGN");
        let hits = standard().scan(seq.as_bytes());
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].start, 8);
        assert_eq!(hits[0].end, 8 + core.len());
        assert_eq!(hits[0].partial_edge_bases, 0);
    }

    #[test]
    fn gaps_and_ambiguous_bases_split_arrays() {
        let block = "TTGGGTTACACTTGGGTTAC".repeat(3);
        for gap in ["N", "NNNN", "CCCCCCC", "T"] {
            let seq = format!("{block}{gap}{block}");
            let hits = standard().scan(seq.as_bytes());
            assert_eq!(hits.len(), 2, "gap {gap}");
            assert_eq!(hits[0].end, block.len());
            assert_eq!(hits[1].start, block.len() + gap.len());
            assert!(hits.iter().all(|h| h.exact_covered_bases == block.len()));
        }
    }

    #[test]
    fn independently_rotated_units_cannot_create_a_mixed_join() {
        let seq = format!("{}{}", "TTGGGTTAC".repeat(3), "CTTGGGTTACA".repeat(3));
        let hits = standard().scan(seq.as_bytes());
        assert!(hits.iter().all(|h| !h.is_mixed));
        assert!(hits.iter().all(|h| h.exact_covered_bases < seq.len()));
    }

    #[test]
    fn pure_units_are_separate_from_mixed_calls() {
        let seq = "TTGGGTTAC".repeat(5);
        let hits = standard().scan(seq.as_bytes());
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].composition, "pure");
        assert_eq!(hits[0].unit_counts, [0, 4]);
        assert_eq!(hits[0].unit_counts_max, [0, 5]);
        assert!(!hits[0].is_mixed);
    }

    #[test]
    fn two_copies_per_component_is_an_explicit_threshold() {
        let seq = format!("N{}ACTTGGGTTACN", "TTGGGTTAC".repeat(4));
        let hits = standard().scan(seq.as_bytes());
        assert_eq!(hits[0].composition, "mixed_below_component_threshold");
        assert!(!hits[0].is_mixed);
        let relaxed = scanner(
            &["TTGGGTTAC", "ACTTGGGTTAC"],
            MixedConfig {
                min_component_copies: 1,
                ..MixedConfig::default()
            },
        )
        .scan(seq.as_bytes());
        assert!(relaxed[0].is_mixed);
    }

    #[test]
    fn held_out_units_order_and_case_are_invariant() {
        let seq = "GTTCAGACCTG".repeat(5);
        let a = scanner(&["GTTC", "AGACCTG"], MixedConfig::default());
        let b = scanner(&["agacctg", "gttc"], MixedConfig::default());
        assert_eq!(a.definitions, b.definitions);
        assert_eq!(a.scan(seq.as_bytes()), b.scan(seq.as_bytes()));
        let hits = a.scan(seq.as_bytes());
        assert_eq!(hits.len(), 1);
        assert!(hits[0].is_mixed);
        assert_eq!(hits[0].unit_counts, [5, 5]);
    }

    #[test]
    fn validation_rejects_nonprimitive_redundant_and_composite_units() {
        for units in [
            vec!["TTAGGG", "TTAGGGTTAGGG"],
            vec!["TTAGGG", "TAGGGT"],
            vec!["TTAGGG", "CCCTAA"],
            vec!["ACG", "TCA", "ACGTCA"],
            vec!["ACN", "TG"],
        ] {
            assert!(MixedScanner::new(
                vec![MixedDefinition {
                    id: "bad".into(),
                    units: units.into_iter().map(String::from).collect()
                }],
                MixedConfig::default()
            )
            .is_err());
        }
        assert!(MixedScanner::new(vec![], MixedConfig::default()).is_err());
        let same = MixedDefinition {
            id: "same".into(),
            units: vec!["ACG".into(), "TCA".into()],
        };
        assert!(MixedScanner::new(vec![same.clone(), same], MixedConfig::default()).is_err());
    }

    #[test]
    fn reverse_complement_symmetric_sequence_has_no_confident_direction() {
        let scanner = scanner(&["AT", "CG"], MixedConfig::default());
        let hits = scanner.scan(b"ATATCGCGATATCGCGATATCGCG");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, None);
        assert!(hits[0].orientation_ambiguous);
        assert!(hits[0].is_mixed);
    }

    #[test]
    fn short_empty_and_non_dna_records_are_not_hits() {
        for seq in [b"".as_slice(), b"ACGT", b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNN"] {
            assert!(standard().scan(seq).is_empty());
        }
    }

    #[test]
    fn noncomposite_ambiguous_parses_keep_guaranteed_count_bounds() {
        let scanner = scanner(&["ACGT", "ACGTAC", "GTACACGT"], MixedConfig::default());
        let sequence = format!("N{}N", "ACGTACGTACACGT".repeat(4));
        let hits = scanner.scan(sequence.as_bytes());
        let hit = hits
            .iter()
            .find(|hit| hit.start == 1 && hit.end == sequence.len() - 1)
            .unwrap();
        assert!(hit.segmentation_ambiguous);
        // b+c=a+b+a; the apparent second component is not guaranteed.
        assert_eq!(hit.unit_counts, [0, 4, 0]);
        assert_eq!(hit.unit_counts_max, [8, 4, 4]);
        assert_eq!(hit.composition, "composition_ambiguous");
        assert!(!hit.is_mixed);
    }

    #[test]
    fn internal_arrays_have_unambiguous_complete_counts() {
        let seq = format!("N{}N", "TTGGGTTACACTTGGGTTAC".repeat(4));
        let hits = standard().scan(seq.as_bytes());
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].unit_counts, [4, 4]);
        assert_eq!(hits[0].unit_counts_max, [4, 4]);
        assert!(!hits[0].segmentation_ambiguous);
        assert!(!hits[0].core_bounds_ambiguous);
    }

    #[test]
    fn million_base_pure_array_yields_one_maximal_hit() {
        let sequence = "TTGGGTTAC".repeat(120_000);
        let hits = standard().scan(sequence.as_bytes());
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].end, sequence.len());
        assert_eq!(hits[0].unit_counts, [0, 119_999]);
        assert_eq!(hits[0].unit_counts_max, [0, 120_000]);
    }
}
