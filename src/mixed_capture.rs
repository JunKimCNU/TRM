//! Optional, bounded capture of qualifying mixed records in source order.
//!
//! A buffer holds a prefix of eligible records and the identity/length of its
//! first excluded record. Keeping that boundary witness is essential: after a
//! base cap blocks a long record, a later short record must not fill the hole.
//! Worker prefixes can be merged in any order. If a worker omitted a record,
//! its earlier local prefix already exhausts a cap; inserting additional earlier
//! records from other workers can only move the global stopping point earlier.
//! Thus no globally needed record can have been discarded by a worker.
//!
//! The caps bound sequence bases and record counts, not total process RSS.
//! Quality payload is at most another max_bases; IDs are bounded separately.
//! Each worker retains its own bounded buffer; merging can briefly retain one
//! additional candidate. Hit metadata depends on definitions and record length,
//! and JSON construction/serialization has additional allocations. Ordinary
//! input parsing and the existing scanner still accept whole large records.
use crate::{
    input,
    mixed::{MixedDefinition, MixedTract},
};
use serde::Serialize;
use std::collections::BTreeMap;

pub(crate) const MAX_ID_BYTES: usize = 4096;
type Key = (usize, u64);

#[derive(Clone, Debug, Serialize)]
pub(crate) struct Options {
    pub enabled: bool,
    pub max_records: usize,
    pub max_bases: usize,
    pub max_record_bases: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            enabled: false,
            max_records: 1000,
            max_bases: 1_000_000,
            max_record_bases: 100_000,
        }
    }
}
impl Options {
    pub fn validate(&self) -> Result<(), String> {
        if !(1..=1_000_000).contains(&self.max_records)
            || !(1..=1_000_000_000).contains(&self.max_bases)
            || !(1..=10_000_000).contains(&self.max_record_bases)
            || self.max_record_bases > self.max_bases
        {
            return Err("mixed capture requires max-records 1..1000000, max-bases 1..1000000000, and max-record-bases 1..10000000 <= max-bases".into());
        }
        Ok(())
    }
}

/// Numeric group index avoids copying arbitrarily long definition IDs into
/// every captured hit. The index addresses mixed_arrays.definitions.
#[derive(Debug, Serialize)]
struct Hit {
    group_index: usize,
    start: usize,
    end: usize,
    core_start: usize,
    core_end: usize,
    strand: Option<bool>,
    unit_counts: Vec<usize>,
    unit_counts_max: Vec<usize>,
    segmentation_ambiguous: bool,
    core_bounds_ambiguous: bool,
    orientation_ambiguous: bool,
    exact_covered_bases: usize,
    partial_edge_bases: usize,
    distinct_units: usize,
}
#[derive(Debug, Serialize)]
struct CapturedRecord {
    input_index: usize,
    record_ordinal: u64,
    record_id: String,
    record_length: usize,
    sequence: String,
    quality: Option<String>,
    quality_status: &'static str,
    hits: Vec<Hit>,
}
impl CapturedRecord {
    fn key(&self) -> Key {
        (self.input_index, self.record_ordinal)
    }
}
#[derive(Clone, Copy, Debug, Serialize)]
struct Boundary {
    input_index: usize,
    record_ordinal: u64,
    record_length: usize,
}
impl Boundary {
    fn key(&self) -> Key {
        (self.input_index, self.record_ordinal)
    }
}
#[derive(Default, Debug, Serialize)]
struct Counts {
    positive_records_in_analyzed_input: u64,
    eligible_positive_records: u64,
    eligible_positive_bases: u64,
    omitted_oversize_record_records: u64,
    omitted_oversize_id_records: u64,
}

pub(crate) struct Buffer {
    options: Options,
    records: BTreeMap<Key, CapturedRecord>,
    bases: usize,
    boundary: Option<Boundary>,
    counts: Counts,
}
impl Buffer {
    pub fn new(options: Options) -> Self {
        Self {
            options,
            records: BTreeMap::new(),
            bases: 0,
            boundary: None,
            counts: Counts::default(),
        }
    }
    pub fn observe(
        &mut self,
        input_index: usize,
        ordinal: u64,
        record: input::Record,
        hits: Vec<MixedTract>,
        definitions: &[MixedDefinition],
    ) {
        if !hits.iter().any(|h| h.is_mixed) {
            return;
        }
        self.counts.positive_records_in_analyzed_input += 1;
        if record.seq.len() > self.options.max_record_bases {
            self.counts.omitted_oversize_record_records += 1;
            return;
        }
        if record.id.len() > MAX_ID_BYTES {
            self.counts.omitted_oversize_id_records += 1;
            return;
        }
        self.counts.eligible_positive_records += 1;
        self.counts.eligible_positive_bases += record.seq.len() as u64;
        let key = (input_index, ordinal);
        if self.boundary.is_some_and(|b| key >= b.key()) {
            return;
        }
        // Streaming workers receive ordered batches. Avoid allocating hit
        // metadata after saturation; retain only this small boundary witness.
        if self
            .records
            .last_key_value()
            .map_or(true, |(last, _)| key > *last)
            && (self.records.len() == self.options.max_records
                || self.bases + record.seq.len() > self.options.max_bases)
        {
            self.set_boundary(Boundary {
                input_index,
                record_ordinal: ordinal,
                record_length: record.seq.len(),
            });
            return;
        }
        let quality_status = if record.quality.is_some() {
            "fastq_characters_preserved"
        } else {
            "absent_in_fasta"
        };
        let captured_hits = hits
            .into_iter()
            .filter(|h| h.is_mixed)
            .map(|h| Hit {
                group_index: definitions
                    .iter()
                    .position(|d| d.id == h.definition_id)
                    .expect("scanner definition exists"),
                start: h.start,
                end: h.end,
                core_start: h.core_start,
                core_end: h.core_end,
                strand: h.strand,
                unit_counts: h.unit_counts,
                unit_counts_max: h.unit_counts_max,
                segmentation_ambiguous: h.segmentation_ambiguous,
                core_bounds_ambiguous: h.core_bounds_ambiguous,
                orientation_ambiguous: h.orientation_ambiguous,
                exact_covered_bases: h.exact_covered_bases,
                partial_edge_bases: h.partial_edge_bases,
                distinct_units: h.distinct_units,
            })
            .collect();
        self.retain(CapturedRecord {
            input_index,
            record_ordinal: ordinal,
            record_length: record.seq.len(),
            record_id: record.id,
            sequence: String::from_utf8(record.seq).expect("validated ASCII DNA"),
            quality: record
                .quality
                .map(|q| String::from_utf8(q).expect("validated ASCII quality")),
            quality_status,
            hits: captured_hits,
        });
    }
    fn set_boundary(&mut self, boundary: Boundary) {
        if self.boundary.map_or(true, |old| boundary.key() < old.key()) {
            self.boundary = Some(boundary);
            while self
                .records
                .last_key_value()
                .is_some_and(|(key, _)| *key >= boundary.key())
            {
                let (_, record) = self.records.pop_last().unwrap();
                self.bases -= record.record_length;
            }
        }
    }
    fn retain(&mut self, record: CapturedRecord) {
        if self.boundary.is_some_and(|b| record.key() >= b.key()) {
            return;
        }
        self.bases += record.record_length;
        assert!(
            self.records.insert(record.key(), record).is_none(),
            "source ordinal must be unique within an input"
        );
        while self.records.len() > self.options.max_records || self.bases > self.options.max_bases {
            let (_, last) = self.records.pop_last().unwrap();
            self.bases -= last.record_length;
            self.set_boundary(Boundary {
                input_index: last.input_index,
                record_ordinal: last.record_ordinal,
                record_length: last.record_length,
            });
        }
    }
    pub fn merge(&mut self, other: Self) {
        self.counts.positive_records_in_analyzed_input +=
            other.counts.positive_records_in_analyzed_input;
        self.counts.eligible_positive_records += other.counts.eligible_positive_records;
        self.counts.eligible_positive_bases += other.counts.eligible_positive_bases;
        self.counts.omitted_oversize_record_records += other.counts.omitted_oversize_record_records;
        self.counts.omitted_oversize_id_records += other.counts.omitted_oversize_id_records;
        if let Some(boundary) = other.boundary {
            self.set_boundary(boundary);
        }
        for (_, record) in other.records {
            self.retain(record);
        }
    }
    pub fn json(
        self,
        sampled_prefix: bool,
        workers: usize,
        batch_bases: usize,
    ) -> serde_json::Value {
        let stored_records = self.records.len();
        let omitted_caps = self.counts.eligible_positive_records - stored_records as u64;
        let missing_quality = self
            .records
            .values()
            .filter(|r| r.quality.is_none())
            .count();
        let cap_reason = self.boundary.map(|b| {
            match (
                stored_records == self.options.max_records,
                self.bases + b.record_length > self.options.max_bases,
            ) {
                (true, true) => "record_and_base_limits",
                (true, false) => "record_limit",
                (false, true) => "base_limit",
                (false, false) => unreachable!("first excluded eligible record must exceed a cap"),
            }
        });
        let complete_analyzed = omitted_caps == 0
            && self.counts.omitted_oversize_record_records == 0
            && self.counts.omitted_oversize_id_records == 0;
        let mut reasons = serde_json::Map::new();
        if let Some(reason) = cap_reason {
            reasons.insert(reason.into(), omitted_caps.into());
        }
        serde_json::json!({
            "schema_version": 1, "enabled": true,
            "settings": {"max_records": self.options.max_records, "max_bases": self.options.max_bases,
                "max_record_bases": self.options.max_record_bases, "max_record_id_bytes": MAX_ID_BYTES},
            "selection": "first eligible positive prefix by (input_index, 1-based record_ordinal); stop at first record/base cap, do not fill later smaller records",
            "input_scope": if sampled_prefix { "analyzed_prefix_only" } else { "complete_input_stream" },
            "source_eof_observed": !sampled_prefix,
            "complete_for_analyzed_positive_records": complete_analyzed,
            "complete_for_all_input_positive_records": complete_analyzed && !sampled_prefix,
            "status": if !complete_analyzed { "incomplete_omissions" } else if sampled_prefix { "complete_for_analyzed_prefix_only" } else { "complete_for_input_positive_records" },
            "counts": self.counts, "stored_records": stored_records, "stored_sequence_bases": self.bases,
            "omitted_after_prefix_cap_records": omitted_caps, "cap_reached": omitted_caps > 0,
            "cap_stop_reason": cap_reason, "omission_counts_by_first_cap_reason": reasons,
            "first_cap_excluded_record": self.boundary,
            "stored_fastq_quality_records": stored_records - missing_quality,
            "stored_fasta_without_quality_records": missing_quality,
            "quality_status": if missing_quality == 0 { "available_for_all_stored_records" } else { "fasta_quality_absent" },
            "record_identity": "input_index is zero-based; record_ordinal is one-based within each input, including nonpositive records; duplicate text IDs remain distinct",
            "sequence_representation": "parser-normalized uppercase IUPAC DNA; not original file bytes, header comments, whitespace or wrapping",
            "quality_representation": "validated FASTQ ASCII quality characters in original order with line endings/wrapping removed; encoding not calibrated; null for FASTA",
            "paired_context": "independent positive records only; no native pairing or automatic mate capture; uncaptured mate is unknown",
            "hits_scope": "qualifying mixed arrays only, group_index addresses mixed_arrays.definitions; pure/nonqualifying blockers require re-evaluation of the retained full sequence",
            "memory_scope": {"worker_count": workers, "batch_bases": batch_bases,
                "per_worker_sequence_base_cap": self.options.max_bases, "per_worker_quality_byte_cap": self.options.max_bases,
                "per_worker_record_cap": self.options.max_records, "queued_batch_capacity": workers * 2,
                "notes": "Each worker holds one capped prefix; merge temporarily adds at most one retained record. At most 3*workers+1 batches plus parser are in flight. Extra queued quality is bounded by (3*workers+1)*(batch_bases+max_record_bases), plus one parser record; oversized records retain no quality. IDs <=4096 bytes per captured record. Hit metadata scales with definitions and captured sequence length (at most 2*definition_count*captured_bases tracts); numeric group indices avoid repeated definition strings. Vec/map overhead, original parser/scanner allocations and JSON serialization are additional. These are payload/structure bounds, not an RSS guarantee; the existing memory policy is heuristic."},
            "records": self.records.into_values().collect::<Vec<_>>()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mixed::{MixedConfig, MixedScanner};
    fn sample(
        ordinal: u64,
        length: usize,
    ) -> (input::Record, Vec<MixedTract>, Vec<MixedDefinition>) {
        let scanner = MixedScanner::new(
            vec![MixedDefinition {
                id: "group".into(),
                units: vec!["TTAGGG".into(), "TTAGGGG".into()],
            }],
            MixedConfig::default(),
        )
        .unwrap();
        let mut seq = b"TTAGGGTTAGGGTTAGGGGTTAGGGG".to_vec();
        seq.resize(length, b'N');
        let hits = scanner.scan(&seq);
        (
            input::Record {
                id: "duplicate".into(),
                seq,
                quality: Some(vec![b'!' + (ordinal % 50) as u8; length]),
            },
            hits,
            scanner.definitions,
        )
    }
    fn partitioned(
        workers: usize,
        reverse: bool,
        max_records: usize,
        max_bases: usize,
    ) -> serde_json::Value {
        let options = Options {
            enabled: true,
            max_records,
            max_bases,
            max_record_bases: 160,
        };
        let mut local: Vec<_> = (0..workers).map(|_| Buffer::new(options.clone())).collect();
        for ordinal in 1..=55_u64 {
            let (record, hits, definitions) = sample(ordinal, 26 + (ordinal * 43 % 173) as usize);
            local[(ordinal as usize * 17) % workers].observe(
                0,
                ordinal,
                record,
                hits,
                &definitions,
            );
        }
        if reverse {
            local.reverse();
        }
        let mut combined = Buffer::new(options);
        for buffer in local {
            combined.merge(buffer);
        }
        let mut json = combined.json(false, workers, 1000);
        json.as_object_mut().unwrap().remove("memory_scope");
        json
    }
    #[test]
    fn worker_prefix_merge_matches_one_worker_with_long_blocker_and_short_followers() {
        for (records, bases) in [(7, 500), (3, 1000), (100, 161), (55, 100_000)] {
            let expected = partitioned(1, false, records, bases);
            for workers in [2, 3, 4, 7, 8] {
                for reverse in [false, true] {
                    assert_eq!(partitioned(workers, reverse, records, bases), expected);
                }
            }
        }
    }
    #[test]
    fn invalid_caps_are_finite_and_fail_before_retention() {
        let mut options = Options {
            max_records: 0,
            ..Options::default()
        };
        assert!(options.validate().is_err());
        options.max_records = 1;
        options.max_record_bases = options.max_bases + 1;
        assert!(options.validate().is_err());
        options.max_record_bases = 1;
        options.max_bases = usize::MAX;
        assert!(options.validate().is_err());
    }
}
