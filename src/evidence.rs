//! Read-boundary evidence is descriptive molecule geometry, not telomere proof.
use crate::detect::Tract;
pub use crate::flank_period::FlankPeriodAnnotation;
use crate::motif::{canonical, is_dna, primitive_len};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::BufReader;

const SIGNATURE_CAPACITY: usize = 128;
const CONTEXT_BASES: usize = 20;
const REPRESENTATIVE_BASES: usize = 80;

/// Stable FNV-1a/128 (non-cryptographic); never Rust's randomized HashMap hash.
/// A read hash collision remains possible. Flank identity instead uses all 40
/// packed DNA bits, with the hash used only for ranking and the packed key for ties.
fn stable_hash(bytes: impl IntoIterator<Item = u8>) -> u128 {
    bytes.into_iter().fold(
        0x6c62_272e_07bb_0142_62b8_2175_6295_c58d_u128,
        |hash, byte| {
            (hash ^ u128::from(byte)).wrapping_mul(0x0000_0000_0100_0000_0000_0000_0000_013b)
        },
    )
}

/// Unlike repeat-motif complementation, read identity preserves IUPAC symbols.
fn read_complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => base,
    }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct QuadrantCounts {
    pub left_forward: u64,
    pub left_reverse: u64,
    pub right_forward: u64,
    pub right_reverse: u64,
    pub ambiguous: u64,
}

impl QuadrantCounts {
    pub(crate) fn flags(left: [bool; 3], right: [bool; 3]) -> Self {
        Self {
            left_forward: u64::from(left[0]),
            left_reverse: u64::from(left[1]),
            right_forward: u64::from(right[0]),
            right_reverse: u64::from(right[1]),
            ambiguous: u64::from(left[2] || right[2]),
        }
    }

    fn reverse_complement(&self) -> Self {
        Self {
            left_forward: self.right_reverse,
            left_reverse: self.right_forward,
            right_forward: self.left_reverse,
            right_reverse: self.left_forward,
            ambiguous: self.ambiguous,
        }
    }

    pub(crate) fn merge(&mut self, other: Self) {
        self.left_forward += other.left_forward;
        self.left_reverse += other.left_reverse;
        self.right_forward += other.right_forward;
        self.right_reverse += other.right_reverse;
        self.ambiguous += other.ambiguous;
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ReadSignature {
    pub signature: String,
    pub records: u64,
    pub quadrant_counts: QuadrantCounts,
    pub normalized_quadrant_counts: QuadrantCounts,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ReadSignatures {
    pub capacity: usize,
    pub observed_records: u64,
    pub distinct_lower_bound: u64,
    pub saturated: bool,
    pub entries: Vec<ReadSignature>,
}

impl Default for ReadSignatures {
    fn default() -> Self {
        Self {
            capacity: SIGNATURE_CAPACITY,
            observed_records: 0,
            distinct_lower_bound: 0,
            saturated: false,
            entries: Vec::new(),
        }
    }
}

impl ReadSignatures {
    fn insert(&mut self, entry: ReadSignature) {
        match self
            .entries
            .binary_search_by(|old| old.signature.cmp(&entry.signature))
        {
            Ok(index) => {
                let old = &mut self.entries[index];
                old.records += entry.records;
                old.quadrant_counts.merge(entry.quadrant_counts);
                old.normalized_quadrant_counts
                    .merge(entry.normalized_quadrant_counts);
            }
            Err(index) => {
                if self.entries.len() == self.capacity {
                    self.saturated = true;
                }
                if index < self.capacity {
                    self.entries.insert(index, entry);
                    self.entries.truncate(self.capacity);
                }
            }
        }
        self.distinct_lower_bound = self.entries.len() as u64 + u64::from(self.saturated);
    }

    fn observe(&mut self, sequence: &[u8], quadrants: QuadrantCounts) {
        self.observed_records += 1;
        let reverse = sequence
            .iter()
            .copied()
            .cmp(sequence.iter().rev().map(|&base| read_complement(base)))
            == std::cmp::Ordering::Greater;
        let hash = if reverse {
            stable_hash(sequence.iter().rev().map(|&base| read_complement(base)))
        } else {
            stable_hash(sequence.iter().copied())
        };
        let normalized = if reverse {
            quadrants.reverse_complement()
        } else {
            quadrants.clone()
        };
        self.insert(ReadSignature {
            signature: format!("{hash:032x}"),
            records: 1,
            quadrant_counts: quadrants,
            normalized_quadrant_counts: normalized,
        });
    }

    fn merge(&mut self, other: Self) {
        self.observed_records += other.observed_records;
        self.saturated |= other.saturated;
        for entry in other.entries {
            self.insert(entry);
        }
        self.distinct_lower_bound = self.entries.len() as u64 + u64::from(self.saturated);
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct FlankContext {
    /// Collision-free 40-bit encoding of physical sequence20: A=0,C=1,G=2,T=3.
    pub signature: String,
    pub rank: String,
    pub sequence: String,
    pub representative_flank: String,
    pub events: u64,
    pub quadrant_counts: QuadrantCounts,
    pub normalized_quadrant_counts: QuadrantCounts,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct FlankContexts {
    pub capacity: usize,
    pub flank_bases: usize,
    pub observed_events: u64,
    pub distinct_lower_bound: u64,
    pub saturated: bool,
    pub entries: Vec<FlankContext>,
}

impl Default for FlankContexts {
    fn default() -> Self {
        Self {
            capacity: SIGNATURE_CAPACITY,
            flank_bases: CONTEXT_BASES,
            observed_events: 0,
            distinct_lower_bound: 0,
            saturated: false,
            entries: Vec::new(),
        }
    }
}

impl FlankContexts {
    fn insert(&mut self, entry: FlankContext) {
        match self.entries.binary_search_by(|old| {
            (&old.rank, &old.signature).cmp(&(&entry.rank, &entry.signature))
        }) {
            Ok(index) => {
                let old = &mut self.entries[index];
                old.events += entry.events;
                old.quadrant_counts.merge(entry.quadrant_counts);
                old.normalized_quadrant_counts
                    .merge(entry.normalized_quadrant_counts);
                if entry.representative_flank.len() > old.representative_flank.len()
                    || (entry.representative_flank.len() == old.representative_flank.len()
                        && entry.representative_flank < old.representative_flank)
                {
                    old.representative_flank = entry.representative_flank;
                }
            }
            Err(index) => {
                if self.entries.len() == self.capacity {
                    self.saturated = true;
                }
                if index < self.capacity {
                    self.entries.insert(index, entry);
                    self.entries.truncate(self.capacity);
                }
            }
        }
        self.distinct_lower_bound = self.entries.len() as u64 + u64::from(self.saturated);
    }

    fn observe(
        &mut self,
        sequence: &[u8],
        tract: &Tract,
        tracts: &[Tract],
        left: bool,
        strand: usize,
    ) {
        let mut context = Vec::with_capacity(REPRESENTATIVE_BASES);
        if left {
            let end = tracts
                .iter()
                .filter(|other| other.start >= tract.end)
                .map(|other| other.start)
                .min()
                .unwrap_or(sequence.len())
                .min(sequence.len())
                .min(tract.end.saturating_add(REPRESENTATIVE_BASES));
            context.extend(
                sequence[tract.end..end]
                    .iter()
                    .copied()
                    .take_while(|base| matches!(base, b'A' | b'C' | b'G' | b'T')),
            );
        } else {
            let start = tracts
                .iter()
                .filter(|other| other.end <= tract.start)
                .map(|other| other.end)
                .max()
                .unwrap_or(0)
                .max(tract.start.saturating_sub(REPRESENTATIVE_BASES));
            context.extend(
                sequence[start..tract.start]
                    .iter()
                    .rev()
                    .copied()
                    .take_while(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
                    .map(read_complement),
            );
        }
        if context.len() < CONTEXT_BASES || !usable_flank(&context[..CONTEXT_BASES]) {
            return;
        }
        self.observed_events += 1;
        let bases = &context[..CONTEXT_BASES];
        let packed = bases.iter().fold(0_u64, |packed, base| {
            (packed << 2)
                | match base {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => unreachable!(),
                }
        });
        let mut flags = [false; 3];
        flags[strand] = true;
        let quadrants = if left {
            QuadrantCounts::flags(flags, [false; 3])
        } else {
            QuadrantCounts::flags([false; 3], flags)
        };
        let normalized = if left {
            quadrants.clone()
        } else {
            quadrants.reverse_complement()
        };
        self.insert(FlankContext {
            signature: format!("{packed:010x}"),
            rank: format!("{:032x}", stable_hash(bases.iter().copied())),
            sequence: String::from_utf8(bases.to_vec()).unwrap(),
            representative_flank: String::from_utf8(context).unwrap(),
            events: 1,
            quadrant_counts: quadrants,
            normalized_quadrant_counts: normalized,
        });
    }

    fn merge(&mut self, other: Self) {
        self.observed_events += other.observed_events;
        self.saturated |= other.saturated;
        for entry in other.entries {
            self.insert(entry);
        }
        self.distinct_lower_bound = self.entries.len() as u64 + u64::from(self.saturated);
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Validation {
    pub schema_version: u32,
    #[serde(default)]
    pub measured: bool,
    pub clean_two_flank_repeat_records: u64,
    pub junction_read_signatures: ReadSignatures,
    pub flank_contexts: FlankContexts,
}

impl Default for Validation {
    fn default() -> Self {
        Self {
            schema_version: 1,
            measured: true,
            clean_two_flank_repeat_records: 0,
            junction_read_signatures: ReadSignatures::default(),
            flank_contexts: FlankContexts::default(),
        }
    }
}

fn unmeasured_validation() -> Validation {
    Validation {
        measured: false,
        ..Validation::default()
    }
}

impl Validation {
    fn merge(&mut self, other: Self) {
        self.measured &= other.measured;
        self.clean_two_flank_repeat_records += other.clean_two_flank_repeat_records;
        self.junction_read_signatures
            .merge(other.junction_read_signatures);
        self.flank_contexts.merge(other.flank_contexts);
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct Options {
    pub flank_bases: usize,
    pub end_slack: usize,
    #[serde(skip_serializing_if = "is_false")]
    pub flank_period_annotation: bool,
}

fn is_false(value: &bool) -> bool {
    !*value
}

impl Default for Options {
    fn default() -> Self {
        Self {
            flank_bases: 20,
            end_slack: 3,
            flank_period_annotation: false,
        }
    }
}

impl Options {
    pub fn validate(&self) -> Result<(), String> {
        if !(20..=10_000).contains(&self.flank_bases) {
            return Err("junction-flank-bases must be 20..10000".into());
        }
        if self.end_slack > 1000 {
            return Err("junction-end-slack must be 0..1000".into());
        }
        Ok(())
    }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct JunctionEvidence {
    pub full_repeat_records: u64,
    pub junction_reads: u64,
    pub left_forward_junction_reads: u64,
    pub left_reverse_junction_reads: u64,
    pub left_ambiguous_junction_reads: u64,
    pub right_forward_junction_reads: u64,
    pub right_reverse_junction_reads: u64,
    pub right_ambiguous_junction_reads: u64,
    pub internal_repeat_reads: u64,
    #[serde(default = "unmeasured_validation")]
    pub validation: Validation,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flank_period_annotation: Option<Box<FlankPeriodAnnotation>>,
}

/// Operational flank screen: contiguous A/C/G/T, >=3 distinct bases, and no
/// `>=85%` lag agreement at periods short enough to contain three repeat copies.
/// Passing does not establish a unique locus or a nonrepetitive genome region.
pub(crate) fn usable_flank(flank: &[u8]) -> bool {
    if !is_dna(flank) || b"ACGT".iter().filter(|base| flank.contains(base)).count() < 3 {
        return false;
    }
    for period in 1..=(flank.len() / 3).min(32) {
        let agreements = (period..flank.len())
            .filter(|&position| flank[position] == flank[position - period])
            .count();
        if agreements * 100 >= (flank.len() - period) * 85 {
            return false;
        }
    }
    true
}

impl JunctionEvidence {
    /// Compatibility entry point: without a candidate motif the new optional
    /// annotation is unmeasured, while all established counters still update.
    #[allow(dead_code)]
    pub fn observe(&mut self, sequence: &[u8], tracts: &[Tract], options: &Options) {
        self.observe_impl(sequence, tracts, options, None);
    }

    pub fn observe_with_motif(
        &mut self,
        sequence: &[u8],
        tracts: &[Tract],
        options: &Options,
        motif: &[u8],
    ) {
        self.observe_impl(sequence, tracts, options, Some(motif));
    }

    /// Initialize a measured empty target explicitly, including targets with
    /// no detected tracts. Earlier unannotated junction observations remain
    /// visibly incomplete rather than being relabeled measured zero.
    pub fn init_flank_period_annotation(&mut self, motif: &[u8]) {
        if let Some(annotation) = &self.flank_period_annotation {
            annotation.assert_motif(motif);
        } else {
            let mut annotation = FlankPeriodAnnotation::new(motif);
            annotation.measured_all_events = self.junction_reads == 0;
            self.flank_period_annotation = Some(Box::new(annotation));
        }
    }

    fn observe_impl(
        &mut self,
        sequence: &[u8],
        tracts: &[Tract],
        options: &Options,
        motif: Option<&[u8]>,
    ) {
        let annotate = options.flank_period_annotation && motif.is_some();
        if annotate {
            self.init_flank_period_annotation(motif.unwrap());
        }
        let full = tracts
            .iter()
            .any(|tract| (tract.end - tract.start) as f64 >= sequence.len() as f64 * 0.9);
        self.full_repeat_records += u64::from(full);
        self.internal_repeat_reads += u64::from(tracts.iter().any(|tract| {
            tract.start > options.end_slack && sequence.len() - tract.end > options.end_slack
        }));
        self.validation.clean_two_flank_repeat_records += u64::from(tracts.iter().any(|tract| {
            tract.start >= options.flank_bases
                && sequence.len() - tract.end >= options.flank_bases
                && usable_flank(&sequence[tract.start - options.flank_bases..tract.start])
                && usable_flank(&sequence[tract.end..tract.end + options.flank_bases])
                && !tracts.iter().any(|other| {
                    (other.start < tract.start && other.end > tract.start - options.flank_bases)
                        || (other.start < tract.end + options.flank_bases && other.end > tract.end)
                })
        }));
        // Even long reads with >=20 bp remaining flank are excluded if >=90%
        // covered by one repeat tract: full-repeat evidence is not directional.
        if full {
            return;
        }
        let mut left = [false; 3];
        let mut right = [false; 3];
        for tract in tracts {
            let strand = match tract.strand {
                Some(true) => 0,
                Some(false) => 1,
                None => 2,
            };
            if tract.start <= options.end_slack
                && sequence.len() - tract.end >= options.flank_bases
                && usable_flank(&sequence[tract.end..tract.end + options.flank_bases])
                && !tracts.iter().any(|other| {
                    other.start < tract.end + options.flank_bases && other.end > tract.end
                })
            {
                left[strand] = true;
                self.annotate_flank_event(sequence, tract, tracts, true, strand, annotate);
                self.validation
                    .flank_contexts
                    .observe(sequence, tract, tracts, true, strand);
            }
            if sequence.len() - tract.end <= options.end_slack
                && tract.start >= options.flank_bases
                && usable_flank(&sequence[tract.start - options.flank_bases..tract.start])
                && !tracts.iter().any(|other| {
                    other.start < tract.start && other.end > tract.start - options.flank_bases
                })
            {
                right[strand] = true;
                self.annotate_flank_event(sequence, tract, tracts, false, strand, annotate);
                self.validation
                    .flank_contexts
                    .observe(sequence, tract, tracts, false, strand);
            }
        }
        self.junction_reads += u64::from(left.iter().chain(&right).any(|flag| *flag));
        self.left_forward_junction_reads += u64::from(left[0]);
        self.left_reverse_junction_reads += u64::from(left[1]);
        self.left_ambiguous_junction_reads += u64::from(left[2]);
        self.right_forward_junction_reads += u64::from(right[0]);
        self.right_reverse_junction_reads += u64::from(right[1]);
        self.right_ambiguous_junction_reads += u64::from(right[2]);
        if left.iter().chain(&right).any(|flag| *flag) {
            self.validation
                .junction_read_signatures
                .observe(sequence, QuadrantCounts::flags(left, right));
        }
    }

    fn annotate_flank_event(
        &mut self,
        sequence: &[u8],
        tract: &Tract,
        tracts: &[Tract],
        left: bool,
        strand: usize,
        annotate: bool,
    ) {
        if let Some(annotation) = &mut self.flank_period_annotation {
            if annotate {
                annotation.observe(sequence, tract, tracts, left, strand);
            } else {
                annotation.measured_all_events = false;
            }
        }
    }

    pub fn merge(&mut self, other: Self) {
        match (
            &mut self.flank_period_annotation,
            other.flank_period_annotation,
        ) {
            (Some(annotation), Some(incoming)) => annotation.merge(*incoming),
            (Some(annotation), None) => {
                annotation.measured_all_events &= other.junction_reads == 0;
            }
            (slot @ None, Some(mut incoming)) => {
                incoming.measured_all_events &= self.junction_reads == 0;
                *slot = Some(incoming);
            }
            (None, None) => {}
        }
        self.full_repeat_records += other.full_repeat_records;
        self.junction_reads += other.junction_reads;
        self.left_forward_junction_reads += other.left_forward_junction_reads;
        self.left_reverse_junction_reads += other.left_reverse_junction_reads;
        self.left_ambiguous_junction_reads += other.left_ambiguous_junction_reads;
        self.right_forward_junction_reads += other.right_forward_junction_reads;
        self.right_reverse_junction_reads += other.right_reverse_junction_reads;
        self.right_ambiguous_junction_reads += other.right_ambiguous_junction_reads;
        self.internal_repeat_reads += other.internal_repeat_reads;
        self.validation.merge(other.validation);
    }
}

#[derive(Clone, Debug)]
pub struct QualitySummary {
    pub records: u64,
    pub bases: u64,
    pub min_record_bases: usize,
    pub max_record_bases: usize,
    pub ambiguous_bases: u64,
    pub records_with_ambiguity: u64,
    exact_lengths: Vec<u64>,
    long_lengths: Vec<u64>,
}

impl Default for QualitySummary {
    fn default() -> Self {
        Self {
            records: 0,
            bases: 0,
            min_record_bases: 0,
            max_record_bases: 0,
            ambiguous_bases: 0,
            records_with_ambiguity: 0,
            exact_lengths: vec![0; 1001],
            long_lengths: vec![0; usize::BITS as usize],
        }
    }
}

impl QualitySummary {
    pub fn observe(&mut self, sequence: &[u8]) {
        let length = sequence.len();
        self.records += 1;
        self.bases += length as u64;
        if self.records == 1 || length < self.min_record_bases {
            self.min_record_bases = length;
        }
        self.max_record_bases = self.max_record_bases.max(length);
        let ambiguous = sequence
            .iter()
            .filter(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
            .count() as u64;
        self.ambiguous_bases += ambiguous;
        self.records_with_ambiguity += u64::from(ambiguous > 0);
        if length <= 1000 {
            self.exact_lengths[length] += 1;
        } else {
            self.long_lengths[length.ilog2() as usize] += 1;
        }
    }

    pub fn merge(&mut self, other: Self) {
        if self.records == 0
            || (other.records > 0 && other.min_record_bases < self.min_record_bases)
        {
            self.min_record_bases = other.min_record_bases;
        }
        self.records += other.records;
        self.bases += other.bases;
        self.max_record_bases = self.max_record_bases.max(other.max_record_bases);
        self.ambiguous_bases += other.ambiguous_bases;
        self.records_with_ambiguity += other.records_with_ambiguity;
        for (left, right) in self.exact_lengths.iter_mut().zip(other.exact_lengths) {
            *left += right;
        }
        for (left, right) in self.long_lengths.iter_mut().zip(other.long_lengths) {
            *left += right;
        }
    }

    pub fn json(&self) -> Value {
        let mut bins = Vec::new();
        for (length, &count) in self.exact_lengths.iter().enumerate() {
            if count > 0 {
                bins.push(serde_json::json!({"lower": length, "upper": length, "records": count}));
            }
        }
        for (exponent, &count) in self.long_lengths.iter().enumerate() {
            if count > 0 {
                let lower = (1usize << exponent).max(1001);
                let upper = 1usize
                    .checked_shl(exponent as u32 + 1)
                    .map(|value| value - 1)
                    .unwrap_or(usize::MAX);
                bins.push(serde_json::json!({"lower": lower, "upper": upper, "records": count}));
            }
        }
        serde_json::json!({
            "records": self.records, "bases": self.bases,
            "min_record_bases": self.min_record_bases, "max_record_bases": self.max_record_bases,
            "ambiguous_bases": self.ambiguous_bases, "records_with_ambiguity": self.records_with_ambiguity,
            "length_histogram": bins,
            "length_histogram_method": "Exact inclusive bins through 1000 bp; bounded inclusive power-of-two bins above 1000 bp.",
            "quality_scores_assessed": false, "duplicates_assessed": false,
        })
    }
}

pub fn load_targets(path: &str) -> Result<Vec<String>, String> {
    let file =
        File::open(path).map_err(|error| format!("cannot open targets {path:?}: {error}"))?;
    let value: Value = serde_json::from_reader(BufReader::new(file))
        .map_err(|error| format!("invalid targets JSON: {error}"))?;
    let raw = if let Some(motifs) = value.get("motifs") {
        motifs
            .as_array()
            .ok_or("targets motifs must be an array")?
            .clone()
    } else {
        value
            .get("candidates")
            .and_then(Value::as_array)
            .ok_or("targets require motifs array or saved candidates array")?
            .iter()
            .map(|candidate| candidate["motif"].clone())
            .collect()
    };
    let mut families = BTreeSet::new();
    for item in raw {
        let motif = item
            .as_str()
            .ok_or("target motifs must be DNA strings")?
            .to_ascii_uppercase();
        if motif.is_empty() || motif.len() > 64 || !is_dna(motif.as_bytes()) {
            return Err("target motifs require 1..64 A/C/G/T bases".into());
        }
        let unit = &motif.as_bytes()[..primitive_len(motif.as_bytes())];
        families.insert(canonical(unit).0);
        if families.len() > 256 {
            return Err("target limit exceeded: at most 256 unique primitive families".into());
        }
    }
    if families.is_empty() {
        return Err("at least one target motif is required".into());
    }
    Ok(families.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    const CLEAN: &[u8] = b"AGTCGATCGGACTACGATCA";

    fn tract(start: usize, end: usize, strand: bool) -> Tract {
        Tract {
            start,
            end,
            identity: 1.0,
            strand: Some(strand),
        }
    }

    fn rc(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|&base| read_complement(base))
            .collect()
    }

    #[test]
    fn optional_period_annotation_preserves_legacy_json_and_counters() {
        assert_eq!(
            std::mem::size_of::<Option<Box<FlankPeriodAnnotation>>>(),
            std::mem::size_of::<usize>()
        );
        let off_options = Options::default();
        assert!(serde_json::to_value(&off_options)
            .unwrap()
            .get("flank_period_annotation")
            .is_none());
        let old_options: Options =
            serde_json::from_value(serde_json::json!({"flank_bases":20,"end_slack":3})).unwrap();
        assert!(!old_options.flank_period_annotation);
        let on_options = Options {
            flank_period_annotation: true,
            ..Options::default()
        };
        let motif = b"AACCCCT";
        let mut off = JunctionEvidence::default();
        let mut on = JunctionEvidence::default();
        // The old 20 bp predicate cannot test period seven. The new annotation
        // requires all 21 adjacent bases; this intentionally injected tract
        // geometry tests evidence independently of tract extension behavior.
        assert!(usable_flank(&motif.repeat(3)[..20]));
        assert!(!usable_flank(&motif.repeat(3)));
        for flank in [
            motif.repeat(3)[..20].to_vec(),
            motif.repeat(3),
            CLEAN.to_vec(),
        ] {
            let sequence = [motif.repeat(10), flank].concat();
            let tracts = [tract(0, 70, true)];
            off.observe_with_motif(&sequence, &tracts, &off_options, motif);
            on.observe_with_motif(&sequence, &tracts, &on_options, motif);
        }
        let annotation = on.flank_period_annotation.as_ref().unwrap();
        assert_eq!(annotation.observed_events, 3);
        assert_eq!(annotation.classes.insufficient_window.events, 2);
        assert_eq!(annotation.classes.both.events, 1);
        assert!(annotation.measured_all_events);
        assert_eq!(on.junction_reads, 3);
        assert!(off.flank_period_annotation.is_none());
        let old_json = serde_json::to_value(&off).unwrap();
        let mut new_json = serde_json::to_value(&on).unwrap();
        new_json
            .as_object_mut()
            .unwrap()
            .remove("flank_period_annotation");
        assert_eq!(old_json, new_json);
        assert_eq!(
            serde_json::from_value::<JunctionEvidence>(serde_json::to_value(on).unwrap())
                .unwrap()
                .junction_reads,
            3
        );
    }

    #[test]
    fn annotation_eligibility_still_obeys_legacy_full_mask_and_both_boundary_events() {
        let motif = b"AACCCCT";
        let options = Options {
            flank_period_annotation: true,
            ..Options::default()
        };
        let sequence = [motif.repeat(100), CLEAN.to_vec()].concat();
        let mut full = JunctionEvidence::default();
        full.observe_with_motif(&sequence, &[tract(0, 700, true)], &options, motif);
        assert_eq!(
            full.flank_period_annotation
                .as_ref()
                .unwrap()
                .observed_events,
            0
        );
        let sequence = [motif.repeat(10), CLEAN.to_vec(), motif.repeat(10)].concat();
        let mut both = JunctionEvidence::default();
        both.observe_with_motif(
            &sequence,
            &[tract(0, 70, true), tract(90, 160, false)],
            &options,
            motif,
        );
        let annotation = both.flank_period_annotation.unwrap();
        assert_eq!(both.junction_reads, 1);
        assert_eq!(annotation.observed_events, 2);
        assert_eq!(annotation.classes.insufficient_window.events, 2);
        assert_eq!(
            annotation
                .classes
                .insufficient_window
                .quadrant_counts
                .left_forward,
            1
        );
        assert_eq!(
            annotation
                .classes
                .insufficient_window
                .quadrant_counts
                .right_reverse,
            1
        );
        let mut overlap = JunctionEvidence::default();
        overlap.observe_with_motif(
            &sequence,
            &[tract(0, 70, true), tract(80, 160, false)],
            &options,
            motif,
        );
        assert_eq!(overlap.flank_period_annotation.unwrap().observed_events, 0);
    }

    #[test]
    fn annotation_right_flank_keeps_physical_order_and_reverse_complement_quadrants() {
        let motif = b"AACCCCT";
        let options = Options {
            flank_period_annotation: true,
            ..Options::default()
        };
        let sequence = [motif.repeat(10), motif.repeat(3)].concat();
        let mut evidence = JunctionEvidence::default();
        evidence.observe_with_motif(&sequence, &[tract(0, 70, true)], &options, motif);
        evidence.observe_with_motif(
            &rc(&sequence),
            &[tract(21, 91, false)],
            &options,
            b"TTAGGGG",
        );
        let annotation = evidence.flank_period_annotation.unwrap();
        assert_eq!(annotation.classes.both.events, 2);
        assert_eq!(annotation.classes.both.quadrant_counts.left_forward, 1);
        assert_eq!(annotation.classes.both.quadrant_counts.right_reverse, 1);
    }

    #[test]
    fn initialized_zero_and_partially_measured_merges_are_distinct() {
        let mut zero = JunctionEvidence::default();
        zero.init_flank_period_annotation(b"TTAGGGG");
        assert!(
            zero.flank_period_annotation
                .as_ref()
                .unwrap()
                .measured_all_events
        );
        assert_eq!(
            zero.flank_period_annotation
                .as_ref()
                .unwrap()
                .observed_events,
            0
        );
        let mut absent = JunctionEvidence::default();
        absent.merge(zero.clone());
        assert_eq!(absent, zero);
        let mut measured_first = zero.clone();
        measured_first.merge(JunctionEvidence::default());
        assert_eq!(measured_first, zero);
        let sequence = [b"AACCCCT".repeat(10), CLEAN.to_vec()].concat();
        let mut unmeasured = JunctionEvidence::default();
        unmeasured.observe(&sequence, &[tract(0, 70, true)], &Options::default());
        assert_eq!(unmeasured.junction_reads, 1);
        let mut left = zero.clone();
        left.merge(unmeasured.clone());
        let mut right = unmeasured.clone();
        right.merge(zero);
        assert_eq!(left, right);
        assert!(
            !left
                .flank_period_annotation
                .as_ref()
                .unwrap()
                .measured_all_events
        );
        unmeasured.init_flank_period_annotation(b"TTAGGGG");
        assert!(
            !unmeasured
                .flank_period_annotation
                .as_ref()
                .unwrap()
                .measured_all_events
        );
        let mut missing_motif = JunctionEvidence::default();
        missing_motif.init_flank_period_annotation(b"TTAGGGG");
        missing_motif.observe(
            &sequence,
            &[tract(0, 70, true)],
            &Options {
                flank_period_annotation: true,
                ..Options::default()
            },
        );
        assert!(
            !missing_motif
                .flank_period_annotation
                .unwrap()
                .measured_all_events
        );
    }

    #[test]
    fn annotation_counts_eligible_event_even_when_old_context_retention_drops_it() {
        let sequence = [b"AACCCCT".repeat(10), vec![b'A'; 20], CLEAN.to_vec()].concat();
        let mut evidence = JunctionEvidence::default();
        evidence.observe_with_motif(
            &sequence,
            &[tract(0, 70, true)],
            &Options {
                flank_bases: 40,
                flank_period_annotation: true,
                ..Options::default()
            },
            b"AACCCCT",
        );
        assert_eq!(evidence.junction_reads, 1);
        assert_eq!(evidence.validation.flank_contexts.observed_events, 0);
        assert_eq!(evidence.flank_period_annotation.unwrap().observed_events, 1);
    }

    #[test]
    fn reverse_complement_reads_and_physical_flanks_share_signatures() {
        let sequence = [b"AACCCT".repeat(10), CLEAN.to_vec()].concat();
        let reverse = rc(&sequence);
        let mut evidence = JunctionEvidence::default();
        evidence.observe(&sequence, &[tract(0, 60, true)], &Options::default());
        evidence.observe(&reverse, &[tract(20, 80, false)], &Options::default());
        assert_eq!(evidence.junction_reads, 2);
        assert_eq!(evidence.left_forward_junction_reads, 1);
        assert_eq!(evidence.right_reverse_junction_reads, 1);
        let reads = &evidence.validation.junction_read_signatures;
        assert_eq!(reads.observed_records, 2);
        assert_eq!(reads.distinct_lower_bound, 1);
        assert_eq!(reads.entries[0].records, 2);
        assert_eq!(reads.entries[0].quadrant_counts.left_forward, 1);
        assert_eq!(reads.entries[0].quadrant_counts.right_reverse, 1);
        assert_eq!(reads.entries[0].normalized_quadrant_counts.left_forward, 2);
        let contexts = &evidence.validation.flank_contexts;
        assert_eq!(contexts.observed_events, 2);
        assert_eq!(contexts.distinct_lower_bound, 1);
        assert_eq!(contexts.entries[0].sequence.as_bytes(), CLEAN);
        assert_eq!(contexts.entries[0].events, 2);
        assert_eq!(
            contexts.entries[0].normalized_quadrant_counts.left_forward,
            2
        );
    }

    #[test]
    fn full_repeat_duplicates_are_never_junction_signatures() {
        let mut evidence = JunctionEvidence::default();
        for _ in 0..100 {
            evidence.observe(
                &b"AACCCT".repeat(20),
                &[tract(0, 120, true)],
                &Options::default(),
            );
        }
        assert_eq!(evidence.full_repeat_records, 100);
        assert_eq!(
            evidence
                .validation
                .junction_read_signatures
                .observed_records,
            0
        );
        assert!(evidence
            .validation
            .junction_read_signatures
            .entries
            .is_empty());
        assert_eq!(evidence.validation.flank_contexts.observed_events, 0);
        assert_eq!(evidence.validation.clean_two_flank_repeat_records, 0);
    }

    #[test]
    fn wider_legacy_flank_does_not_fabricate_clean_twenty_base_context() {
        let sequence = [b"AACCCT".repeat(10), vec![b'A'; 20], CLEAN.to_vec()].concat();
        let mut evidence = JunctionEvidence::default();
        evidence.observe(
            &sequence,
            &[tract(0, 60, true)],
            &Options {
                flank_bases: 40,
                ..Options::default()
            },
        );
        assert_eq!(evidence.junction_reads, 1);
        assert_eq!(
            evidence
                .validation
                .junction_read_signatures
                .observed_records,
            1
        );
        assert_eq!(evidence.validation.flank_contexts.observed_events, 0);
    }

    #[test]
    fn two_boundary_events_still_profile_one_junction_read() {
        let sequence = [b"AACCCT".repeat(10), CLEAN.repeat(2), b"AACCCT".repeat(10)].concat();
        let mut evidence = JunctionEvidence::default();
        evidence.observe(
            &sequence,
            &[tract(0, 60, true), tract(100, 160, true)],
            &Options::default(),
        );
        assert_eq!(evidence.junction_reads, 1);
        assert_eq!(
            evidence
                .validation
                .junction_read_signatures
                .observed_records,
            1
        );
        assert_eq!(
            evidence.validation.junction_read_signatures.entries[0].records,
            1
        );
        assert_eq!(evidence.validation.flank_contexts.observed_events, 2);
    }

    #[test]
    fn clean_two_flank_evidence_does_not_change_legacy_internal_count() {
        let repeat = b"AACCCT".repeat(10);
        let sequences = [
            [CLEAN.to_vec(), repeat.clone(), CLEAN.to_vec()].concat(),
            [b"AGTC".to_vec(), repeat.clone(), b"AGTC".to_vec()].concat(),
            [vec![b'A'; 20], repeat.clone(), CLEAN.to_vec()].concat(),
        ];
        let mut evidence = JunctionEvidence::default();
        for (sequence, start) in sequences.iter().zip([20, 4, 20]) {
            evidence.observe(
                sequence,
                &[tract(start, start + 60, true)],
                &Options::default(),
            );
        }
        assert_eq!(evidence.internal_repeat_reads, 3);
        assert_eq!(evidence.validation.clean_two_flank_repeat_records, 1);
        assert_eq!(evidence.junction_reads, 0);
        // Multiple qualifying tracts still count only one record.
        let sequence = [
            CLEAN.to_vec(),
            repeat.clone(),
            CLEAN.to_vec(),
            repeat,
            CLEAN.to_vec(),
        ]
        .concat();
        evidence.observe(
            &sequence,
            &[tract(20, 80, true), tract(100, 160, true)],
            &Options::default(),
        );
        assert_eq!(evidence.validation.clean_two_flank_repeat_records, 2);
        // A same-family tract overlaps the only otherwise-clean left window.
        let sequence = sequences[0].clone();
        let mut overlap = JunctionEvidence::default();
        overlap.observe(
            &sequence,
            &[tract(0, 10, true), tract(20, 80, true)],
            &Options::default(),
        );
        assert_eq!(overlap.validation.clean_two_flank_repeat_records, 0);
    }

    #[test]
    fn clean_two_flanks_remain_visible_when_one_tract_covers_ninety_percent() {
        let sequence = [CLEAN.to_vec(), b"AACCCT".repeat(100), CLEAN.to_vec()].concat();
        let mut evidence = JunctionEvidence::default();
        evidence.observe(&sequence, &[tract(20, 620, true)], &Options::default());
        assert_eq!(evidence.full_repeat_records, 1);
        assert_eq!(evidence.validation.clean_two_flank_repeat_records, 1);
        assert_eq!(
            evidence
                .validation
                .junction_read_signatures
                .observed_records,
            0
        );
    }

    #[test]
    fn read_identity_preserves_iupac_symbols() {
        let sequence = b"ACGTRYSWKMBDHVNACGTA";
        assert_eq!(rc(&rc(sequence)), sequence);
        let mut reads = ReadSignatures::default();
        let quadrants = QuadrantCounts {
            left_forward: 1,
            ..QuadrantCounts::default()
        };
        reads.observe(sequence, quadrants.clone());
        reads.observe(&rc(sequence), quadrants.reverse_complement());
        reads.observe(b"ACGTNYSWKMBDHVNACGTA", quadrants);
        assert_eq!(reads.observed_records, 3);
        assert_eq!(reads.distinct_lower_bound, 2);
        assert!(reads.entries.iter().any(|entry| entry.records == 2));
    }

    #[test]
    fn legacy_deserialization_is_unmeasured_and_new_empty_measurement_is_measured() {
        let legacy: JunctionEvidence =
            serde_json::from_value(serde_json::json!({"junction_reads": 4})).unwrap();
        assert!(!legacy.validation.measured);
        let empty = JunctionEvidence::default();
        assert!(empty.validation.measured);
        assert_eq!(
            serde_json::from_value::<JunctionEvidence>(serde_json::to_value(&empty).unwrap())
                .unwrap(),
            empty
        );
    }

    fn random_sequence(mut seed: u64, length: usize) -> Vec<u8> {
        (0..length)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed & 3) as usize]
            })
            .collect()
    }

    #[test]
    fn saturated_bottom_k_counts_are_exact_and_merge_order_independent() {
        let mut all = ReadSignatures::default();
        let mut groups = vec![ReadSignatures::default(); 7];
        let mut truth = std::collections::BTreeMap::<String, u64>::new();
        for index in 1..=450 {
            let sequence = random_sequence(index, 100);
            let normalized = sequence.clone().min(rc(&sequence));
            let signature = format!("{:032x}", stable_hash(normalized));
            for occurrence in 0..(index % 7 + 1) {
                let q = QuadrantCounts {
                    right_forward: 1,
                    ..QuadrantCounts::default()
                };
                all.observe(&sequence, q.clone());
                groups[((index + occurrence) % 7) as usize].observe(&sequence, q);
                *truth.entry(signature.clone()).or_default() += 1;
            }
        }
        assert!(all.saturated);
        assert_eq!(all.entries.len(), SIGNATURE_CAPACITY);
        assert_eq!(all.distinct_lower_bound, SIGNATURE_CAPACITY as u64 + 1);
        let expected: Vec<_> = truth
            .iter()
            .take(SIGNATURE_CAPACITY)
            .map(|(key, count)| (key.clone(), *count))
            .collect();
        assert_eq!(
            all.entries
                .iter()
                .map(|entry| (entry.signature.clone(), entry.records))
                .collect::<Vec<_>>(),
            expected
        );
        let mut forward = ReadSignatures::default();
        for group in groups.clone() {
            forward.merge(group);
        }
        let mut reverse = ReadSignatures::default();
        for group in groups.into_iter().rev() {
            reverse.merge(group);
        }
        assert_eq!(all, forward);
        assert_eq!(all, reverse);
    }

    #[test]
    fn saturated_flank_context_merge_keeps_exact_counts_and_representatives() {
        let mut all = FlankContexts::default();
        let mut groups = vec![FlankContexts::default(); 5];
        let mut truth = std::collections::BTreeMap::<(String, String), (u64, String)>::new();
        for index in 1..=400 {
            let flank = random_sequence(index, 80);
            if !usable_flank(&flank[..20]) {
                continue;
            }
            for occurrence in 0..3 {
                let n = [20, 40, 80][occurrence];
                let sequence = [b"AACCCT".repeat(10), flank[..n].to_vec()].concat();
                let tract = tract(0, 60, true);
                all.observe(&sequence, &tract, std::slice::from_ref(&tract), true, 0);
                groups[(index as usize + occurrence) % 5].observe(
                    &sequence,
                    &tract,
                    std::slice::from_ref(&tract),
                    true,
                    0,
                );
                let mut single = FlankContexts::default();
                single.observe(&sequence, &tract, std::slice::from_ref(&tract), true, 0);
                let entry = &single.entries[0];
                let target = truth
                    .entry((entry.rank.clone(), entry.signature.clone()))
                    .or_default();
                target.0 += 1;
                if entry.representative_flank.len() > target.1.len() {
                    target.1 = entry.representative_flank.clone();
                }
            }
        }
        assert!(all.saturated);
        assert_eq!(all.entries.len(), 128);
        for (entry, ((rank, signature), (count, representative))) in
            all.entries.iter().zip(truth.iter())
        {
            assert_eq!(
                (
                    &entry.rank,
                    &entry.signature,
                    entry.events,
                    &entry.representative_flank
                ),
                (rank, signature, *count, representative)
            );
        }
        let mut merged = FlankContexts::default();
        for group in groups.into_iter().rev() {
            merged.merge(group);
        }
        assert_eq!(merged, all);
    }

    #[test]
    fn representative_is_longest_then_lexical_and_stops_at_ambiguity_or_repeat() {
        let mut sketch = FlankContexts::default();
        for suffix in [b"TGCT".as_slice(), b"ACGT".as_slice()] {
            let sequence = [
                b"AACCCT".repeat(10),
                CLEAN.to_vec(),
                suffix.to_vec(),
                b"NCCCC".to_vec(),
            ]
            .concat();
            sketch.observe(
                &sequence,
                &tract(0, 60, true),
                &[tract(0, 60, true)],
                true,
                0,
            );
        }
        assert_eq!(
            sketch.entries[0].representative_flank.as_bytes(),
            [CLEAN, b"ACGT"].concat()
        );
        let sequence = [b"AACCCT".repeat(10), CLEAN.to_vec(), b"AACCCT".repeat(10)].concat();
        let mut bounded = FlankContexts::default();
        bounded.observe(
            &sequence,
            &tract(0, 60, true),
            &[tract(0, 60, true), tract(80, 140, true)],
            true,
            0,
        );
        assert_eq!(bounded.entries[0].representative_flank.as_bytes(), CLEAN);
    }

    #[test]
    fn nonrepeat_flanks_reject_ambiguity_and_short_periods() {
        assert!(usable_flank(b"AGTCGATCGGACTACGATCA"));
        assert!(!usable_flank(b"AGTCGATCNNNACTACGATCA"));
        assert!(!usable_flank(b"AAAAAAAAAAAAAAAAAAAA"));
        assert!(!usable_flank(b"TTAGGCTTAGGCTTAGGCTTA"));
    }

    #[test]
    fn histogram_is_bounded_and_counts_every_record() {
        let mut qc = QualitySummary::default();
        for length in [20, 101, 1000, 1001, 1023, 1024, 2048] {
            qc.observe(&vec![b'A'; length]);
        }
        qc.observe(b"ACGTNN");
        let value = qc.json();
        assert_eq!(value["records"], 8);
        assert_eq!(value["ambiguous_bases"], 2);
        assert_eq!(value["min_record_bases"], 6);
        let bins = value["length_histogram"].as_array().unwrap();
        assert_eq!(
            bins.iter()
                .map(|bin| bin["records"].as_u64().unwrap())
                .sum::<u64>(),
            8
        );
        assert!(bins
            .iter()
            .any(|bin| bin["lower"] == 1001 && bin["upper"] == 1023 && bin["records"] == 2));
        assert!(bins
            .iter()
            .any(|bin| bin["lower"] == 1024 && bin["upper"] == 2047));
    }
}
