//! Explicit, empirical repeat-background screening after the sequence prior.
//!
//! Read-fraction Wilson bounds are an uncertainty heuristic: reads and mates
//! are dependent, so this module does not report p-values or an FDR. No motif
//! identity is prohibited; every family is evaluated with the same equations.

use crate::candidate::{prior_tier, same_evidence};
use crate::motif::{canonical, is_dna, primitive_len};
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};

const WARNING_PREFIX: &str = "Background support: ";
const DETECTOR_KEYS: [&str; 5] = ["min_k", "max_k", "seed_bases", "min_copies", "min_identity"];

#[derive(Clone, Debug, Serialize)]
pub struct SupportOptions {
    pub background_model: Option<String>,
    pub policy: String,
    pub min_enrichment: f64,
}

impl Default for SupportOptions {
    fn default() -> Self {
        Self {
            background_model: None,
            policy: "annotate".into(),
            min_enrichment: 2.0,
        }
    }
}

impl SupportOptions {
    pub fn validate(&self) -> Result<(), String> {
        if !["filter", "annotate"].contains(&self.policy.as_str()) {
            return Err("support policy must be filter or annotate".into());
        }
        if !self.min_enrichment.is_finite() || self.min_enrichment < 1.0 {
            return Err("support min-enrichment must be finite and at least 1".into());
        }
        if self
            .background_model
            .as_ref()
            .is_some_and(|path| path.trim().is_empty())
        {
            return Err("background model path must not be empty".into());
        }
        Ok(())
    }
}

fn unsigned(value: &Value, key: &str, context: &str) -> Result<u64, String> {
    value[key]
        .as_u64()
        .ok_or_else(|| format!("{context}: {key} must be an unsigned integer"))
}

fn nonempty(value: &Value, key: &str, context: &str) -> Result<String, String> {
    value[key]
        .as_str()
        .filter(|value| !value.trim().is_empty())
        .map(str::to_string)
        .ok_or_else(|| format!("{context}: {key} must be a nonempty string"))
}

fn validated_config(value: &Value, context: &str) -> Result<(), String> {
    if !value.is_object() {
        return Err(format!("{context}: config must be an object"));
    }
    let min_k = unsigned(value, "min_k", context)?;
    let max_k = unsigned(value, "max_k", context)?;
    if min_k == 0 || min_k > max_k || max_k > 64 {
        return Err(format!("{context}: require 1 <= min_k <= max_k <= 64"));
    }
    for key in ["seed_bases", "min_copies", "min_reads"] {
        if unsigned(value, key, context)? == 0 {
            return Err(format!("{context}: {key} must be positive"));
        }
    }
    if unsigned(value, "min_copies", context)? < 2
        || unsigned(value, "min_copies", context)? > 1_000_000
        || unsigned(value, "seed_bases", context)? > 100_000_000
    {
        return Err(format!(
            "{context}: invalid min_copies or seed_bases detector bounds"
        ));
    }
    if value.get("max_candidates").is_some() && unsigned(value, "max_candidates", context)? == 0 {
        return Err(format!("{context}: max_candidates must be positive"));
    }
    value["min_identity"]
        .as_f64()
        .filter(|identity| identity.is_finite() && (0.5..=1.0).contains(identity))
        .ok_or_else(|| format!("{context}: min_identity must be finite in [0.5, 1]"))?;
    if !matches!(value["mode"].as_str(), Some("reads" | "assembly")) {
        return Err(format!("{context}: mode must be reads or assembly"));
    }
    Ok(())
}

fn totals(value: &Value, context: &str) -> Result<(u64, u64, u64), String> {
    let records = unsigned(value, "records", context)?;
    let bases = unsigned(value, "bases", context)?;
    let max_record = unsigned(value, "max_record_bases", context)?;
    if records == 0 || bases == 0 || max_record == 0 || max_record > bases {
        return Err(format!(
            "{context}: records, bases and max_record_bases must be positive and consistent"
        ));
    }
    if bases as u128 > records as u128 * max_record as u128 {
        return Err(format!(
            "{context}: bases exceeds records * max_record_bases"
        ));
    }
    Ok((records, bases, max_record))
}

fn family(sequence: &str) -> Result<String, String> {
    if sequence.is_empty() || sequence.len() > 64 || !is_dna(sequence.as_bytes()) {
        return Err(format!("invalid background/candidate motif {sequence:?}"));
    }
    let primitive = &sequence.as_bytes()[..primitive_len(sequence.as_bytes())];
    Ok(canonical(primitive).0)
}

#[derive(Clone, Copy)]
struct Counts {
    reads: u64,
    repeat_bases: u64,
}

fn counts(
    value: &Value,
    records: u64,
    bases: u64,
    max_record: u64,
    context: &str,
) -> Result<Counts, String> {
    if !value.is_object() {
        return Err(format!("{context}: motif counts must be an object"));
    }
    let reads = unsigned(value, "reads", context)?;
    let repeat_bases = unsigned(value, "repeat_bases", context)?;
    if reads == 0 || reads > records || repeat_bases == 0 || repeat_bases > bases {
        return Err(format!(
            "{context}: motif counts must be positive and not exceed input totals"
        ));
    }
    if repeat_bases as u128 > reads as u128 * max_record as u128 {
        return Err(format!(
            "{context}: repeat_bases exceeds reads * max_record_bases"
        ));
    }
    Ok(Counts {
        reads,
        repeat_bases,
    })
}

struct Control {
    id: String,
    records: u64,
    bases: u64,
    max_record: u64,
    min_reads: u64,
    sampled_prefix: bool,
    config: Value,
    motifs: HashMap<String, Counts>,
}

fn validate_model(model: &Value) -> Result<Vec<Control>, String> {
    if model["schema_version"].as_u64() != Some(1)
        || model["model_type"] != "empirical_repeat_background"
    {
        return Err(
            "background model requires schema_version 1 and model_type empirical_repeat_background"
                .into(),
        );
    }
    nonempty(model, "name", "background model")?;
    nonempty(model, "created_utc", "background model")?;
    if !model["limitations"]
        .as_array()
        .is_some_and(|values| values.iter().all(Value::is_string))
    {
        return Err("background model limitations must be an array of strings".into());
    }
    let raw_controls = model["controls"]
        .as_array()
        .filter(|controls| !controls.is_empty())
        .ok_or("background model requires at least one control")?;
    let mut ids = HashSet::new();
    let mut controls = Vec::with_capacity(raw_controls.len());
    for raw in raw_controls {
        let id = nonempty(raw, "control_id", "background control")?;
        let context = format!("background control {id}");
        if !ids.insert(id.clone()) {
            return Err(format!("duplicate background control_id {id}"));
        }
        nonempty(raw, "source_result", &context)?;
        let hash = nonempty(raw, "source_result_sha256", &context)?;
        if hash.len() != 64 || !hash.bytes().all(|byte| byte.is_ascii_hexdigit()) {
            return Err(format!(
                "{context}: source_result_sha256 must contain 64 hex digits"
            ));
        }
        validated_config(&raw["config"], &context)?;
        if raw["config"]["mode"] != "reads" {
            return Err(format!(
                "{context}: background controls must use reads mode"
            ));
        }
        let (records, bases, max_record) = totals(raw, &context)?;
        let sampled_prefix = raw["sampled_prefix"]
            .as_bool()
            .ok_or_else(|| format!("{context}: sampled_prefix must be boolean"))?;
        if sampled_prefix {
            return Err(format!("{context}: prefix-sampled controls are not supported; use a complete or explicitly randomized control dataset"));
        }
        let raw_motifs = raw["motifs"]
            .as_object()
            .ok_or_else(|| format!("{context}: motifs must be an object"))?;
        let min_reads = unsigned(&raw["config"], "min_reads", &context)?;
        let min_k = unsigned(&raw["config"], "min_k", &context)? as usize;
        let max_k = unsigned(&raw["config"], "max_k", &context)? as usize;
        let mut motifs = HashMap::with_capacity(raw_motifs.len());
        for (motif, value) in raw_motifs {
            if family(motif)? != *motif || motif.len() < min_k || motif.len() > max_k {
                return Err(format!("{context}: motif {motif:?} must be a canonical primitive DNA family within min_k..max_k"));
            }
            let observed = counts(value, records, bases, max_record, &context)?;
            if observed.reads < min_reads {
                return Err(format!(
                    "{context}: observed motif {motif} is below min_reads"
                ));
            }
            motifs.insert(motif.clone(), observed);
        }
        controls.push(Control {
            id,
            records,
            bases,
            max_record,
            min_reads,
            sampled_prefix,
            config: raw["config"].clone(),
            motifs,
        });
    }
    Ok(controls)
}

fn read_model(path: &str) -> Result<(Value, Vec<Control>), String> {
    let raw = std::fs::read(path)
        .map_err(|error| format!("cannot read background model {path}: {error}"))?;
    let model: Value = serde_json::from_slice(&raw)
        .map_err(|error| format!("invalid background model {path}: {error}"))?;
    let controls = validate_model(&model)?;
    Ok((model, controls))
}

fn compare_config(config: &Value, controls: &[Control]) -> Result<(), String> {
    validated_config(config, "sample")?;
    if config["mode"] != "assembly" {
        for control in controls {
            for key in DETECTOR_KEYS {
                if config[key].as_f64() != control.config[key].as_f64() {
                    return Err(format!(
                        "background control {} detector mismatch for {key}: sample {}, control {}",
                        control.id, config[key], control.config[key]
                    ));
                }
            }
        }
    }
    Ok(())
}

/// Validate an explicit model and detector compatibility before scanning reads.
pub fn validate_config(config: &Value, opts: &SupportOptions) -> Result<(), String> {
    opts.validate()?;
    if let Some(path) = &opts.background_model {
        let (_, controls) = read_model(path)?;
        compare_config(config, &controls)?;
    }
    Ok(())
}

/// Conventional z=1.96 Wilson bounds, used only as a transparent heuristic.
fn wilson(successes: u64, trials: u64) -> (f64, f64) {
    let n = trials as f64;
    let p = successes as f64 / n;
    let z_squared = 1.96_f64.powi(2);
    let denominator = 1.0 + z_squared / n;
    let center = (p + z_squared / (2.0 * n)) / denominator;
    let radius = (z_squared * (p * (1.0 - p) / n + z_squared / (4.0 * n * n))).sqrt() / denominator;
    ((center - radius).max(0.0), (center + radius).min(1.0))
}

fn clear_annotations(report: &mut Value) {
    let Some(object) = report.as_object_mut() else {
        return;
    };
    object.remove("support_model");
    object.remove("support_options");
    for key in ["candidates", "rejected_candidates"] {
        if let Some(rows) = object.get_mut(key).and_then(Value::as_array_mut) {
            for row in rows {
                if let Some(row) = row.as_object_mut() {
                    row.remove("support");
                    row.remove("candidate_eligible");
                }
            }
        }
    }
    if let Some(summary) = object.get_mut("summary").and_then(Value::as_object_mut) {
        for key in [
            "sequence_prior_candidates",
            "background_supported_candidates",
            "joint_eligible_candidates",
        ] {
            summary.remove(key);
        }
    }
    if let Some(warnings) = object.get_mut("warnings").and_then(Value::as_array_mut) {
        warnings.retain(|warning| {
            !warning
                .as_str()
                .is_some_and(|warning| warning.starts_with(WARNING_PREFIX))
        });
    }
}

/// Apply an explicitly provided model without reading FASTQ or changing counts.
/// All validation is performed before replacing any fields in the input report.
pub fn apply(report: &mut Value, opts: &SupportOptions) -> Result<(), String> {
    opts.validate()?;
    let Some(path) = &opts.background_model else {
        clear_annotations(report);
        return Ok(());
    };
    let (model, controls) = read_model(path)?;
    if !report.is_object() || !report["summary"].is_object() {
        return Err("support requires a report object with a summary object".into());
    }
    compare_config(&report["config"], &controls)?;
    let (records, bases, max_record) = totals(&report["summary"], "sample")?;
    let assembly = report["config"]["mode"] == "assembly";
    let prior_policy = report["prior_options"]["policy"]
        .as_str()
        .filter(|policy| ["canonical", "telomere", "discovery"].contains(policy))
        .ok_or("support requires prior_options.policy canonical, telomere or discovery")?;
    let candidates = report["candidates"]
        .as_array()
        .ok_or("support requires a candidates array")?;
    let rejected: &[Value] = match report.get("rejected_candidates") {
        Some(rows) => rows
            .as_array()
            .ok_or("support rejected_candidates must be an array")?,
        None => &[],
    };
    let mut seen: HashMap<String, &Value> = HashMap::new();
    let mut rows = Vec::with_capacity(candidates.len() + rejected.len());
    for candidate in candidates.iter().chain(rejected) {
        if !candidate.is_object() {
            return Err("support candidate must be an object".into());
        }
        let motif = nonempty(candidate, "motif", "candidate")?;
        let family = family(&motif)?;
        if let Some(previous) = seen.get(&family) {
            if !same_evidence(previous, candidate) {
                return Err(format!("conflicting duplicate candidate family {family}"));
            }
            continue;
        }
        seen.insert(family.clone(), candidate);
        let prior_eligible = candidate["telomere_prior"]["eligible"]
            .as_bool()
            .ok_or("support requires a boolean telomere_prior.eligible for each candidate")?;
        if candidate["telomere_prior"]["discovery_rank"]
            .as_u64()
            .is_none()
        {
            return Err("support requires telomere_prior.discovery_rank for each candidate".into());
        }
        let observed = counts(candidate, records, bases, max_record, "sample candidate")?;
        let full_repeat_fraction = match candidate.get("full_repeat_reads") {
            Some(count) => {
                let count = count
                    .as_u64()
                    .filter(|count| *count <= observed.reads)
                    .ok_or("full_repeat_reads must be unsigned and at most reads")?;
                Some(count as f64 / observed.reads as f64)
            }
            None => None,
        };
        let mut annotation = json!({
            "status": "not_applicable",
            "eligible": true,
            "reads_per_million": observed.reads as f64 / records as f64 * 1e6,
            "repeat_bases_per_million": observed.repeat_bases as f64 / bases as f64 * 1e6,
            "sample_read_fraction_lower": null,
            "background_read_fraction_upper": null,
            "read_enrichment_lower": null,
            "smoothed_background_base_fraction": null,
            "repeat_base_enrichment": null,
            "background_motif_observed": null,
            "reasons": ["assembly_not_comparable_to_read_background"],
        });
        if let Some(fraction) = full_repeat_fraction {
            annotation["full_repeat_fraction"] = json!(fraction);
        }
        if !assembly {
            let sample_lower = wilson(observed.reads, records).0;
            let mut background_upper: f64 = 0.0;
            let mut background_base: f64 = 0.0;
            let mut censored_controls = Vec::new();
            let mut background_observed = false;
            let mut max_censored_reads = 0;
            for control in &controls {
                let background = match control.motifs.get(&family) {
                    Some(counts) => {
                        background_observed = true;
                        *counts
                    }
                    None => {
                        censored_controls.push(control.id.clone());
                        let reads = control.min_reads.saturating_sub(1).min(control.records);
                        max_censored_reads = max_censored_reads.max(reads);
                        Counts {
                            reads,
                            repeat_bases: (reads as u128 * control.max_record as u128)
                                .min(control.bases as u128)
                                as u64,
                        }
                    }
                };
                background_upper =
                    background_upper.max(wilson(background.reads, control.records).1);
                let smoothed_base = (background.repeat_bases as f64 + control.max_record as f64)
                    / (control.bases as f64 + control.max_record as f64);
                background_base = background_base.max(smoothed_base);
            }
            let read_enrichment = sample_lower / background_upper;
            let base_enrichment = observed.repeat_bases as f64 / bases as f64 / background_base;
            let mut reasons = Vec::new();
            if read_enrichment < opts.min_enrichment {
                reasons.push("read_fraction_enrichment_below_threshold");
            }
            if base_enrichment < opts.min_enrichment {
                reasons.push("repeat_base_enrichment_below_threshold");
            }
            let eligible = reasons.is_empty();
            let status = if eligible {
                "enriched"
            } else if !background_observed && read_enrichment < opts.min_enrichment {
                reasons.push("calibration_depth_limits_unseen_motif");
                "insufficient_background_resolution"
            } else if !background_observed {
                "insufficient_support"
            } else {
                "background_compatible"
            };
            annotation["status"] = json!(status);
            annotation["eligible"] = json!(eligible);
            annotation["sample_read_fraction_lower"] = json!(sample_lower);
            annotation["background_read_fraction_upper"] = json!(background_upper);
            annotation["read_enrichment_lower"] = json!(read_enrichment);
            annotation["smoothed_background_base_fraction"] = json!(background_base);
            annotation["repeat_base_enrichment"] = json!(base_enrichment);
            annotation["reasons"] = json!(reasons);
            annotation["censored_controls"] = json!(censored_controls);
            annotation["background_motif_observed"] = json!(background_observed);
            annotation["max_background_censored_reads"] = json!(max_censored_reads);
        }
        let mut row = candidate.clone();
        row["candidate_eligible"] = json!(
            prior_eligible
                && (annotation["eligible"] == true
                    || (prior_policy == "canonical" && opts.policy == "annotate"))
        );
        row["support"] = annotation;
        rows.push(row);
    }
    // Restore the sequence prior's stable ordering even when filter/annotate
    // policies have moved rows between the primary and rejected arrays.
    rows.sort_by_key(|row| row["telomere_prior"]["discovery_rank"].as_u64().unwrap());
    if prior_policy != "canonical" {
        rows.sort_by_key(prior_tier);
    }
    let raw_count = rows.len();
    let prior_count = rows
        .iter()
        .filter(|row| row["telomere_prior"]["eligible"] == true)
        .count();
    let supported_count = rows
        .iter()
        .filter(|row| row["support"]["status"] == "enriched")
        .count();
    let joint_count = rows
        .iter()
        .filter(|row| {
            row["telomere_prior"]["eligible"] == true && row["support"]["eligible"] == true
        })
        .count();
    let eligible_count = if prior_policy == "canonical" && opts.policy == "annotate" {
        prior_count
    } else {
        joint_count
    };
    let (rows, rejected_rows): (Vec<Value>, Vec<Value>) =
        if prior_policy != "discovery" && opts.policy == "filter" {
            rows.into_iter()
                .partition(|row| row["candidate_eligible"] == true)
        } else {
            let rejected_rows = rows
                .iter()
                .filter(|row| row["candidate_eligible"] != true)
                .cloned()
                .collect();
            if prior_policy != "discovery" {
                rows.retain(|row| row["telomere_prior"]["eligible"] == true);
            }
            (rows, rejected_rows)
        };
    let mut warnings: Vec<Value> = match report.get("warnings") {
        Some(value) => value
            .as_array()
            .ok_or("report warnings must be an array")?
            .iter()
            .filter(|warning| {
                !warning
                    .as_str()
                    .is_some_and(|warning| warning.starts_with(WARNING_PREFIX))
            })
            .cloned()
            .collect(),
        None => Vec::new(),
    };
    warnings.push(json!(format!("{WARNING_PREFIX}Wilson read-fraction bounds (z=1.96) are an uncertainty heuristic, not p-values or FDR: paired reads, PCR duplicates and repeat loci are dependent. Base-fraction smoothing is not a confidence interval.")));
    warnings.push(json!(format!("{WARNING_PREFIX}Normalized read/base support is not genomic coverage, telomere length or proof of chromosome-terminal location; comparison to these controls is not a universal absence/background model across the tree of life.")));
    if assembly {
        warnings.push(json!(format!("{WARNING_PREFIX}Assembly records are not comparable to WGS read fractions; background screening is not applicable and sequence-prior eligibility is retained.")));
    } else {
        let mean_read = bases as f64 / records as f64;
        if report["summary"]["sampled_prefix"] == true {
            warnings.push(json!(format!("{WARNING_PREFIX}Sample was scanned as a prefix, which is not random sampling and can bias normalized support.")));
        }
        for control in &controls {
            let mean_control = control.bases as f64 / control.records as f64;
            if (mean_read / mean_control - 1.0).abs() > 0.2 {
                warnings.push(json!(format!("{WARNING_PREFIX}Mean read length differs by more than 20% from control {} (sample {:.1} bp, control {:.1} bp); read-length-dependent detection limits applicability. Means, not medians, are available in these reports.", control.id, mean_read, mean_control)));
            }
            if control.sampled_prefix {
                warnings.push(json!(format!("{WARNING_PREFIX}Control {} was scanned as a prefix sample, which may not represent the full library.", control.id)));
            }
        }
    }
    let mut support_model = model;
    support_model["evaluation_method"] = json!({
        "read_fraction": "sample Wilson lower / maximum control Wilson upper; z=1.96; dependence-aware interpretation required; not a p-value or FDR",
        "unobserved_motif": "censored reads=min(min_reads-1, records); repeat bases=min(censored reads * max_record_bases, bases)",
        "base_fraction": "sample repeat_bases/bases divided by maximum control (repeat_bases+max_record_bases)/(bases+max_record_bases); smoothing, not a confidence interval",
        "decision": if prior_policy == "canonical" && opts.policy == "annotate" {
            "Background enrichment is annotation only; insufficient control resolution does not erase observed canonical candidates. Reported identity is not terminal proof."
        } else { "both enrichment ratios must reach min_enrichment; final eligibility also requires sequence prior" },
    });
    let support_options = serde_json::to_value(opts).map_err(|error| error.to_string())?;
    let reported_count = rows.len();
    // Commit only after every model, row and report validation has succeeded.
    report["schema_version"] = json!(3);
    report["candidates"] = Value::Array(rows);
    report["rejected_candidates"] = Value::Array(rejected_rows);
    report["support_model"] = support_model;
    report["support_options"] = support_options;
    report["warnings"] = Value::Array(warnings);
    let summary = report["summary"].as_object_mut().unwrap();
    summary.insert("raw_reported_candidates".into(), json!(raw_count));
    summary.insert("reported_candidates".into(), json!(reported_count));
    summary.insert("sequence_prior_candidates".into(), json!(prior_count));
    summary.insert(
        "background_supported_candidates".into(),
        json!(supported_count),
    );
    summary.insert("joint_eligible_candidates".into(), json!(joint_count));
    summary.insert("eligible_telomere_candidates".into(), json!(eligible_count));
    summary.insert(
        "rejected_repeat_candidates".into(),
        json!(raw_count - eligible_count),
    );
    summary.insert(
        "call".into(),
        json!(if eligible_count > 0 {
            "candidates_only"
        } else {
            "no_call"
        }),
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_annotation_retains_finite_control_resolution_signal() {
        let file = ModelFile::new(&model());
        let mut value = report("AAGCCT", 2, 200, true);
        value["prior_options"]["policy"] = json!("canonical");
        let mut options = file.options();
        options.policy = "annotate".into();
        apply(&mut value, &options).unwrap();
        assert_eq!(value["candidates"].as_array().unwrap().len(), 1);
        assert_eq!(
            value["candidates"][0]["support"]["status"],
            "insufficient_background_resolution"
        );
        assert_eq!(value["candidates"][0]["candidate_eligible"], true);
        assert_eq!(value["summary"]["joint_eligible_candidates"], 0);
        assert_eq!(value["summary"]["eligible_telomere_candidates"], 1);
        assert_eq!(value["summary"]["call"], "candidates_only");
    }
    use std::sync::atomic::{AtomicUsize, Ordering};

    static NEXT_ID: AtomicUsize = AtomicUsize::new(0);

    fn config() -> Value {
        json!({"min_k": 3, "max_k": 32, "seed_bases": 24, "min_copies": 3,
            "min_identity": 0.85, "min_reads": 2, "mode": "reads"})
    }

    fn model() -> Value {
        json!({"schema_version": 1, "model_type": "empirical_repeat_background",
            "name": "unit test", "created_utc": "2026-09-20T00:00:00Z", "limitations": [],
            "controls": [{"control_id": "negative", "source_result": "test.json",
                "source_result_sha256": "a".repeat(64), "config": config(),
                "records": 100000, "bases": 15000000, "max_record_bases": 150,
                "sampled_prefix": false, "motifs": {"AAT": {"reads": 1000, "repeat_bases": 100000}}}]})
    }

    fn report(motif: &str, reads: u64, repeat_bases: u64, eligible: bool) -> Value {
        json!({"schema_version": 2, "config": config(), "prior_options": {"policy": "telomere"},
            "summary": {"records": 100000, "bases": 15000000, "max_record_bases": 150,
                "raw_reported_candidates": 1, "reported_candidates": 1}, "warnings": [],
            "candidates": [{"motif": motif, "reads": reads, "repeat_bases": repeat_bases,
                "full_repeat_reads": reads / 2, "telomere_prior": {"eligible": eligible,
                    "confidence_tier": "known_motif", "discovery_rank": 0}}],
            "rejected_candidates": []})
    }

    struct ModelFile(std::path::PathBuf);

    impl ModelFile {
        fn new(model: &Value) -> Self {
            let path = std::env::temp_dir().join(format!(
                "trm-support-unit-{}-{}.json",
                std::process::id(),
                NEXT_ID.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::write(&path, serde_json::to_vec(model).unwrap()).unwrap();
            Self(path)
        }

        fn options(&self) -> SupportOptions {
            SupportOptions {
                background_model: Some(self.0.to_string_lossy().into()),
                policy: "filter".into(),
                ..SupportOptions::default()
            }
        }
    }

    impl Drop for ModelFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    #[test]
    fn wilson_boundary_and_known_values() {
        let (lower, upper) = wilson(50, 100);
        assert!((lower - 0.40383).abs() < 0.00001);
        assert!((upper - 0.59617).abs() < 0.00001);
        assert_eq!(wilson(0, 100).0, 0.0);
        assert!((wilson(0, 100).1 - 0.0369948).abs() < 1e-6);
        assert!((wilson(100, 100).1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn high_abundance_background_fails_without_motif_ban() {
        let file = ModelFile::new(&model());
        let mut value = report("AAT", 1100, 110000, true);
        apply(&mut value, &file.options()).unwrap();
        assert!(value["candidates"].as_array().unwrap().is_empty());
        assert_eq!(
            value["rejected_candidates"][0]["support"]["status"],
            "background_compatible"
        );
        assert_eq!(value["summary"]["call"], "no_call");
        let mut enriched = report("AAT", 5000, 500000, true);
        apply(&mut enriched, &file.options()).unwrap();
        assert_eq!(enriched["candidates"][0]["support"]["status"], "enriched");
        assert_eq!(enriched["candidates"][0]["candidate_eligible"], true);
    }

    #[test]
    fn unseen_motif_uses_censored_count_not_zero() {
        let file = ModelFile::new(&model());
        let mut value = report("AACCCT", 100, 10000, true);
        apply(&mut value, &file.options()).unwrap();
        let support = &value["candidates"][0]["support"];
        assert_eq!(
            support["background_read_fraction_upper"].as_f64().unwrap(),
            wilson(1, 100000).1
        );
        assert_eq!(support["censored_controls"], json!(["negative"]));
        assert!(
            (support["smoothed_background_base_fraction"]
                .as_f64()
                .unwrap()
                - 300.0 / 15000150.0)
                .abs()
                < 1e-15
        );
    }

    #[test]
    fn enrichment_never_bypasses_sequence_prior() {
        let file = ModelFile::new(&model());
        let mut value = report("AAC", 1000, 100000, false);
        apply(&mut value, &file.options()).unwrap();
        assert_eq!(value["summary"]["background_supported_candidates"], 1);
        assert_eq!(value["summary"]["eligible_telomere_candidates"], 0);
        assert_eq!(value["rejected_candidates"][0]["support"]["eligible"], true);
        assert_eq!(value["rejected_candidates"][0]["candidate_eligible"], false);
    }

    #[test]
    fn evidence_unchanged_and_annotate_then_filter_is_reversible() {
        let file = ModelFile::new(&model());
        let mut value = report("AAT", 1100, 110000, true);
        let original = value["candidates"][0].clone();
        let mut options = file.options();
        options.policy = "annotate".into();
        apply(&mut value, &options).unwrap();
        assert_eq!(value["candidates"].as_array().unwrap().len(), 1);
        assert_eq!(value["rejected_candidates"].as_array().unwrap().len(), 1);
        assert!(same_evidence(&original, &value["candidates"][0]));
        assert_eq!(
            value["candidates"][0]["support"]["full_repeat_fraction"],
            0.5
        );
        options.policy = "filter".into();
        apply(&mut value, &options).unwrap();
        assert_eq!(value["summary"]["raw_reported_candidates"], 1);
        assert_eq!(value["summary"]["reported_candidates"], 0);
        options.policy = "annotate".into();
        apply(&mut value, &options).unwrap();
        assert_eq!(value["summary"]["reported_candidates"], 1);
        assert!(same_evidence(&original, &value["candidates"][0]));
    }

    #[test]
    fn discovery_policy_retains_all_evidence() {
        let file = ModelFile::new(&model());
        let mut value = report("AAT", 1100, 110000, false);
        value["prior_options"]["policy"] = json!("discovery");
        apply(&mut value, &file.options()).unwrap();
        assert_eq!(value["summary"]["reported_candidates"], 1);
        assert_eq!(value["summary"]["eligible_telomere_candidates"], 0);
        assert_eq!(value["candidates"][0]["candidate_eligible"], false);
    }

    #[test]
    fn mismatch_and_invalid_counts_do_not_partially_mutate_report() {
        let file = ModelFile::new(&model());
        for mut value in [
            report("AAT", 10, 1000, true),
            report("AAT", 100001, 100000, true),
        ] {
            if value["candidates"][0]["reads"] == 10 {
                value["config"]["seed_bases"] = json!(30);
            }
            let original = value.clone();
            assert!(apply(&mut value, &file.options()).is_err());
            assert_eq!(value, original);
        }
    }

    #[test]
    fn assembly_is_not_compared_to_wgs_fractions() {
        let file = ModelFile::new(&model());
        let mut value = report("AAT", 10, 1000, true);
        value["config"]["mode"] = json!("assembly");
        value["config"]["seed_bases"] = json!(30);
        apply(&mut value, &file.options()).unwrap();
        assert_eq!(
            value["candidates"][0]["support"]["status"],
            "not_applicable"
        );
        assert_eq!(value["candidates"][0]["candidate_eligible"], true);
        assert_eq!(value["summary"]["background_supported_candidates"], 0);
        assert_eq!(value["summary"]["eligible_telomere_candidates"], 1);
    }

    #[test]
    fn no_model_clears_stale_annotations_without_repartitioning() {
        let file = ModelFile::new(&model());
        let mut value = report("AAT", 10, 1000, true);
        apply(&mut value, &file.options()).unwrap();
        let old_candidates = value["candidates"].as_array().unwrap().len();
        let old_rejected = value["rejected_candidates"].as_array().unwrap().len();
        let old_count = value["summary"]["eligible_telomere_candidates"].clone();
        apply(&mut value, &SupportOptions::default()).unwrap();
        assert!(value.get("support_model").is_none());
        assert!(value["rejected_candidates"][0].get("support").is_none());
        assert!(value["rejected_candidates"][0]
            .get("candidate_eligible")
            .is_none());
        assert_eq!(
            value["candidates"].as_array().unwrap().len(),
            old_candidates
        );
        assert_eq!(
            value["rejected_candidates"].as_array().unwrap().len(),
            old_rejected
        );
        assert_eq!(value["summary"]["eligible_telomere_candidates"], old_count);
    }

    #[test]
    fn malformed_models_rejected() {
        for modification in 0..6 {
            let mut invalid = model();
            match modification {
                0 => {
                    invalid["controls"][0]["motifs"] =
                        json!({"TTA": {"reads": 5, "repeat_bases": 300}})
                }
                1 => invalid["controls"][0]["motifs"]["AAT"]["reads"] = json!(100001),
                2 => invalid["controls"][0]["config"]["seed_bases"] = Value::Null,
                3 => invalid["controls"][0]["config"]["mode"] = json!("assembly"),
                4 => invalid["controls"] = json!([invalid["controls"][0], invalid["controls"][0]]),
                _ => invalid["controls"][0]["source_result_sha256"] = json!("bad"),
            }
            assert!(validate_model(&invalid).is_err(), "case {modification}");
        }
    }

    #[test]
    fn maximum_background_across_controls_is_used() {
        let mut control_model = model();
        let mut low_background = control_model["controls"][0].clone();
        low_background["control_id"] = json!("low");
        low_background["motifs"] = json!({});
        control_model["controls"]
            .as_array_mut()
            .unwrap()
            .push(low_background);
        let file = ModelFile::new(&control_model);
        let mut value = report("AAT", 1100, 110000, true);
        apply(&mut value, &file.options()).unwrap();
        let support = &value["rejected_candidates"][0]["support"];
        assert_eq!(
            support["background_read_fraction_upper"].as_f64().unwrap(),
            wilson(1000, 100000).1
        );
    }

    #[test]
    fn unseen_weak_signal_is_a_resolution_limit_not_observed_background() {
        let file = ModelFile::new(&model());
        let mut value = report("AACCCT", 2, 300, true);
        apply(&mut value, &file.options()).unwrap();
        let support = &value["rejected_candidates"][0]["support"];
        assert_eq!(support["status"], "insufficient_background_resolution");
        assert_eq!(support["background_motif_observed"], false);
        assert_eq!(support["max_background_censored_reads"], 1);
        assert!(support["reasons"]
            .as_array()
            .unwrap()
            .contains(&json!("calibration_depth_limits_unseen_motif")));
    }

    #[test]
    fn independent_read_and_base_enrichment_are_both_required() {
        let file = ModelFile::new(&model());
        let mut base_fails = report("AAT", 5000, 100000, true);
        apply(&mut base_fails, &file.options()).unwrap();
        let support = &base_fails["rejected_candidates"][0]["support"];
        assert!(support["read_enrichment_lower"].as_f64().unwrap() >= 2.0);
        assert!(support["repeat_base_enrichment"].as_f64().unwrap() < 2.0);
        assert_eq!(
            support["reasons"],
            json!(["repeat_base_enrichment_below_threshold"])
        );
        let mut sparse_model = model();
        sparse_model["controls"][0]["motifs"]["AAT"]["repeat_bases"] = json!(10000);
        let file = ModelFile::new(&sparse_model);
        let mut reads_fail = report("AAT", 1100, 110000, true);
        apply(&mut reads_fail, &file.options()).unwrap();
        let support = &reads_fail["rejected_candidates"][0]["support"];
        assert!(support["read_enrichment_lower"].as_f64().unwrap() < 2.0);
        assert!(support["repeat_base_enrichment"].as_f64().unwrap() >= 2.0);
        assert_eq!(
            support["reasons"],
            json!(["read_fraction_enrichment_below_threshold"])
        );
    }

    #[test]
    fn preserve_known_then_one_edit_then_two_edit_prior_ordering() {
        let file = ModelFile::new(&model());
        let mut value = report("AACCCT", 100, 10000, true);
        let template = value["candidates"][0].clone();
        let mut rows = Vec::new();
        for (rank, motif, class, edits) in [
            (0, "AACACT", "ttaggg_variant", 2),
            (1, "AACCCT", "exact_known", 0),
            (2, "AACCC", "ttaggg_variant", 1),
        ] {
            let mut row = template.clone();
            row["motif"] = json!(motif);
            row["telomere_prior"]["class"] = json!(class);
            row["telomere_prior"]["edit_distance"] = json!(edits);
            row["telomere_prior"]["discovery_rank"] = json!(rank);
            rows.push(row);
        }
        value["candidates"] = Value::Array(rows);
        apply(&mut value, &file.options()).unwrap();
        let motifs: Vec<&str> = value["candidates"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| row["motif"].as_str().unwrap())
            .collect();
        assert_eq!(motifs, ["AACCCT", "AACCC", "AACACT"]);
    }

    #[test]
    fn mean_read_length_and_sample_prefix_warnings_are_explicit() {
        let file = ModelFile::new(&model());
        let mut value = report("AACCCT", 100, 10000, true);
        value["summary"]["bases"] = json!(10100000);
        value["summary"]["max_record_bases"] = json!(101);
        value["summary"]["sampled_prefix"] = json!(true);
        apply(&mut value, &file.options()).unwrap();
        let warnings = value["warnings"].as_array().unwrap();
        assert!(warnings.iter().any(|warning| warning
            .as_str()
            .unwrap()
            .contains("Mean read length differs")));
        assert!(warnings
            .iter()
            .any(|warning| warning.as_str().unwrap().contains("not random sampling")));
        assert_eq!(value["candidates"][0]["support"]["status"], "enriched");
    }

    #[test]
    fn early_validation_fails_missing_model_or_mismatch() {
        let file = ModelFile::new(&model());
        validate_config(&config(), &file.options()).unwrap();
        let mut changed_config = config();
        changed_config["min_identity"] = json!(1.0);
        assert!(validate_config(&changed_config, &file.options())
            .unwrap_err()
            .contains("detector mismatch"));
        let options = SupportOptions {
            background_model: Some("/nonexistent/trm-background-model.json".into()),
            ..SupportOptions::default()
        };
        assert!(validate_config(&config(), &options)
            .unwrap_err()
            .contains("cannot read background model"));
    }
}
