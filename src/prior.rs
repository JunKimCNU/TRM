//! Conservative, auditable biological priors applied after repeat discovery.
//!
//! Repeat evidence is never upgraded to a confirmed telomere by this module.
//! The embedded catalog supplies exact exceptions; approximate matching is only
//! to TTAGGG or to explicitly supplied motifs from related taxa.

use crate::candidate::{prior_tier, same_evidence};
use crate::motif::{canonical, is_dna, primitive_len, reverse_complement};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};

const CATALOG_TSV: &str = include_str!("../data/telomere_motifs.tsv");
const MAX_MOTIF_BASES: usize = 64;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct PriorOptions {
    pub policy: String,
    pub max_edits: usize,
    pub max_distance: f64,
    pub related_motifs: Vec<String>,
    pub lineage: Option<String>,
}

impl Default for PriorOptions {
    fn default() -> Self {
        Self {
            policy: "canonical".into(),
            max_edits: 2,
            max_distance: 0.34,
            related_motifs: Vec::new(),
            lineage: None,
        }
    }
}

impl PriorOptions {
    pub fn validate(&self) -> Result<(), String> {
        if !["canonical", "telomere", "discovery"].contains(&self.policy.as_str()) {
            return Err("prior policy must be canonical, telomere or discovery".into());
        }
        if self.max_edits > MAX_MOTIF_BASES {
            return Err("prior max-edits must be 0..64".into());
        }
        if !self.max_distance.is_finite() || !(0.0..=1.0).contains(&self.max_distance) {
            return Err("prior max-distance must be finite and between 0 and 1".into());
        }
        if let Some(lineage) = &self.lineage {
            if lineage.trim().is_empty() {
                return Err("prior lineage label must not be empty".into());
            }
            if self.related_motifs.is_empty() {
                return Err(
                    "a lineage label requires explicit related motifs; no phylogeny is inferred"
                        .into(),
                );
            }
        }
        for motif in &self.related_motifs {
            normalized_motif(motif)?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize)]
struct CatalogEntry {
    motif: String,
    lineage: String,
    evidence: String,
    source: String,
    canonical_motif: String,
}

fn normalized_motif(value: &str) -> Result<String, String> {
    let motif = value.to_ascii_uppercase();
    if motif.is_empty() || motif.len() > MAX_MOTIF_BASES || !is_dna(motif.as_bytes()) {
        return Err(format!(
            "invalid repeat motif {value:?}: require 1..64 unambiguous A/C/G/T bases"
        ));
    }
    Ok(motif[..primitive_len(motif.as_bytes())].to_string())
}

fn catalog() -> Result<Vec<CatalogEntry>, String> {
    let mut entries = Vec::new();
    let mut columns = None;
    for (line_number, line) in CATALOG_TSV.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if columns.is_none() {
            let positions: Result<Vec<usize>, String> = ["motif", "lineage", "evidence", "source"]
                .iter()
                .map(|name| {
                    fields
                        .iter()
                        .position(|field| field == name)
                        .ok_or_else(|| format!("embedded motif catalog lacks {name:?} column"))
                })
                .collect();
            columns = Some(positions?);
            continue;
        }
        let cols = columns.as_ref().expect("catalog header parsed");
        if cols.iter().any(|&index| index >= fields.len()) {
            return Err(format!(
                "invalid embedded motif catalog row {}",
                line_number + 1
            ));
        }
        let motif = normalized_motif(fields[cols[0]])?;
        if cols[1..].iter().any(|&index| fields[index].is_empty()) {
            return Err(format!(
                "empty embedded motif catalog evidence at row {}",
                line_number + 1
            ));
        }
        entries.push(CatalogEntry {
            canonical_motif: canonical(motif.as_bytes()).0,
            motif,
            lineage: fields[cols[1]].to_string(),
            evidence: fields[cols[2]].to_string(),
            source: fields[cols[3]].to_string(),
        });
    }
    if entries.is_empty() {
        return Err("embedded telomere motif catalog is empty".into());
    }
    Ok(entries)
}

fn base_counts(sequence: &[u8]) -> [usize; 4] {
    let mut counts = [0; 4];
    for base in sequence {
        let i = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!("motifs validated before comparison"),
        };
        counts[i] += 1;
    }
    counts
}

/// Each substitution can repair at most one excess and one missing base.
fn composition_bound(a: &[usize; 4], b: &[usize; 4]) -> usize {
    let excess: usize = a.iter().zip(b).map(|(a, b)| a.saturating_sub(*b)).sum();
    let deficit: usize = a.iter().zip(b).map(|(a, b)| b.saturating_sub(*a)).sum();
    excess.max(deficit)
}

/// Banded Levenshtein DP using rotation offsets, without allocating rotated DNA.
fn bounded_linear_distance(
    a: &[u8],
    a_offset: usize,
    b: &[u8],
    b_offset: usize,
    limit: usize,
) -> Option<usize> {
    if a.len().abs_diff(b.len()) > limit {
        return None;
    }
    let sentinel = limit + 1;
    let mut previous = [sentinel; MAX_MOTIF_BASES + 1];
    let mut current = [sentinel; MAX_MOTIF_BASES + 1];
    for (j, cell) in previous.iter_mut().enumerate().take(b.len().min(limit) + 1) {
        *cell = j;
    }
    for i in 1..=a.len() {
        current.fill(sentinel);
        if i <= limit {
            current[0] = i;
        }
        let start = i.saturating_sub(limit).max(1);
        let end = (i + limit).min(b.len());
        let mut row_min = current[0];
        for j in start..=end {
            let substitution =
                usize::from(a[(a_offset + i - 1) % a.len()] != b[(b_offset + j - 1) % b.len()]);
            current[j] = (previous[j] + 1)
                .min(current[j - 1] + 1)
                .min(previous[j - 1] + substitution);
            row_min = row_min.min(current[j]);
        }
        if row_min > limit {
            return None;
        }
        std::mem::swap(&mut previous, &mut current);
    }
    (previous[b.len()] <= limit).then_some(previous[b.len()])
}

/// Minimum over rotations of BOTH primitive motifs and both relative strands.
/// Rotating both is important when an optimal alignment crosses an indel.
fn circular_distance(a: &[u8], b: &[u8], limit: usize) -> Option<usize> {
    if a.len().abs_diff(b.len()) > limit {
        return None;
    }
    if canonical(a).0 == canonical(b).0 {
        return Some(0);
    }
    let counts = base_counts(a);
    let reverse = reverse_complement(b);
    let mut best = limit + 1;
    for reference in [b, reverse.as_slice()] {
        let lower = composition_bound(&counts, &base_counts(reference));
        if lower >= best {
            continue;
        }
        for a_offset in 0..a.len() {
            for b_offset in 0..reference.len() {
                if let Some(distance) =
                    bounded_linear_distance(a, a_offset, reference, b_offset, best - 1)
                {
                    best = distance;
                    if best == a.len().abs_diff(b.len()).max(1) {
                        return Some(best);
                    }
                    if best == lower {
                        break;
                    }
                }
            }
            if best == lower {
                break;
            }
        }
    }
    (best <= limit).then_some(best)
}

#[derive(Clone)]
struct Match {
    reference: String,
    distance: usize,
    normalized_distance: f64,
    sources: Vec<String>,
    reference_lineages: Vec<String>,
    reference_evidence: Vec<String>,
}

fn qualifying_match(candidate: &str, reference: &str, opts: &PriorOptions) -> Option<Match> {
    let denominator = candidate.len().max(reference.len());
    let limit = opts
        .max_edits
        .min((denominator as f64 * opts.max_distance + 1e-12).floor() as usize);
    let distance = circular_distance(candidate.as_bytes(), reference.as_bytes(), limit)?;
    let normalized_distance = distance as f64 / denominator as f64;
    if normalized_distance > opts.max_distance {
        return None;
    }
    Some(Match {
        reference: reference.into(),
        distance,
        normalized_distance,
        sources: Vec::new(),
        reference_lineages: Vec::new(),
        reference_evidence: Vec::new(),
    })
}

fn add_catalog_evidence(matched: &mut Match, catalog: &[CatalogEntry]) {
    let family = canonical(matched.reference.as_bytes()).0;
    for reference in catalog
        .iter()
        .filter(|entry| entry.canonical_motif == family)
    {
        if !matched.sources.contains(&reference.source) {
            matched.sources.push(reference.source.clone());
        }
        if !matched.reference_lineages.contains(&reference.lineage) {
            matched.reference_lineages.push(reference.lineage.clone());
        }
        if !matched.reference_evidence.contains(&reference.evidence) {
            matched.reference_evidence.push(reference.evidence.clone());
        }
    }
}

fn global_match(
    candidate: &str,
    catalog: &[CatalogEntry],
    opts: &PriorOptions,
) -> Option<(bool, Match)> {
    if let Some(entry) = catalog
        .iter()
        .find(|entry| entry.canonical_motif == candidate)
    {
        let mut matched = Match {
            reference: entry.motif.clone(),
            distance: 0,
            normalized_distance: 0.0,
            sources: Vec::new(),
            reference_lineages: Vec::new(),
            reference_evidence: Vec::new(),
        };
        add_catalog_evidence(&mut matched, catalog);
        return Some((true, matched));
    }
    let mut matched = qualifying_match(candidate, "TTAGGG", opts)?;
    add_catalog_evidence(&mut matched, catalog);
    Some((false, matched))
}

fn established_evidence(evidence: &str) -> bool {
    matches!(
        evidence,
        "experimental"
            | "experimental_terminal"
            | "assembly_terminal"
            | "published_terminal_long_read"
            | "long_read_terminal"
            | "clade_consensus"
    )
}

fn canonical_annotation(candidate: &str, catalog: &[CatalogEntry], opts: &PriorOptions) -> Value {
    let known = catalog
        .iter()
        .find(|entry| entry.canonical_motif == candidate && established_evidence(&entry.evidence));
    let related = opts.related_motifs.iter().find(|motif| {
        canonical(normalized_motif(motif).expect("validated motif").as_bytes()).0 == candidate
    });
    let reference = related
        .map(String::as_str)
        .or_else(|| known.map(|entry| entry.motif.as_str()));
    let references: Vec<_> = catalog
        .iter()
        .filter(|entry| entry.canonical_motif == candidate && established_evidence(&entry.evidence))
        .collect();
    let eligible = reference.is_some();
    json!({
        "class": if related.is_some() { "related_exact" } else if known.is_some() { "exact_known" } else { "uncharacterized_repeat" },
        "eligible": eligible,
        "confidence_tier": if related.is_some() { "related_reference" } else if known.is_some() { "known_reference" } else { "uncharacterized" },
        "matched_reference": reference,
        "observed_motif_display": reference.unwrap_or(candidate),
        "canonical_family": candidate,
        "expected_family_match": if opts.related_motifs.is_empty() { None } else { Some(related.is_some()) },
        "edit_distance": reference.map(|_| 0),
        "normalized_distance": reference.map(|_| 0.0),
        "source_lineage": opts.lineage,
        "evidence_source": references.iter().map(|entry| entry.source.as_str()).collect::<Vec<_>>(),
        "reference_lineages": references.iter().map(|entry| entry.lineage.as_str()).collect::<Vec<_>>(),
        "reference_evidence": references.iter().map(|entry| entry.evidence.as_str()).collect::<Vec<_>>(),
        "explanation": if eligible { "Exact published or explicitly supplied family; a clade expectation does not exclude other known families. Identity is not proof of chromosome-terminal location." } else { "No exact established reference match. Raw repeat evidence is retained for independent novel-candidate assessment; sequence similarity alone does not establish a canonical TRM." },
    })
}

fn annotation(candidate: &str, catalog: &[CatalogEntry], opts: &PriorOptions) -> Value {
    if opts.policy == "canonical" {
        return canonical_annotation(candidate, catalog, opts);
    }
    let global = global_match(candidate, catalog, opts);
    let context = !opts.related_motifs.is_empty();
    let mut related: Option<Match> = None;
    for reference in &opts.related_motifs {
        // validate() has checked every reference before processing report rows.
        let reference = normalized_motif(reference).expect("validated related motif");
        if let Some(mut matched) = qualifying_match(candidate, &reference, opts) {
            let replace = related.as_ref().map_or(true, |best| {
                matched.distance < best.distance
                    || (matched.distance == best.distance
                        && matched.normalized_distance < best.normalized_distance)
            });
            if replace {
                add_catalog_evidence(&mut matched, catalog);
                related = Some(matched);
            }
        }
    }
    let (class, eligible, matched, explanation) = if let Some(matched) = related {
        let class = if matched.distance == 0 {
            "related_exact"
        } else {
            "related_variant"
        };
        (class, true, Some(matched),
            "Matches an explicitly supplied motif from related taxa; this is a sequence prior, not proof of telomeric location.")
    } else if context {
        match global {
            Some((_, matched)) => ("lineage_conflict", false, Some(matched),
                "Matches a global telomere reference but not the explicitly supplied related-taxon references under the configured thresholds."),
            None => ("unsupported_repeat", false, None,
                "Does not match the supplied related-taxon references or the global telomere prior; retained as repeat evidence."),
        }
    } else {
        match global {
            Some((true, matched)) => ("exact_known", true, Some(matched),
                "Exact rotation/strand-equivalent match to a documented telomere motif; species lineage and telomeric location are not inferred."),
            Some((false, matched)) => ("ttaggg_variant", true, Some(matched),
                "Within the configured circular edit-distance thresholds of TTAGGG; variants require independent validation, especially at two or more edits."),
            None => ("unsupported_repeat", false, None,
                "Neither an exact documented telomere motif nor sufficiently close to TTAGGG; retained as repeat evidence."),
        }
    };
    let confidence = match (class, matched.as_ref().map(|m| m.distance)) {
        ("exact_known", _) => "known_reference",
        ("related_exact", _) => "related_reference",
        ("lineage_conflict", _) => "conflicting_reference",
        (_, Some(1)) => "near_reference",
        (_, Some(_)) => "divergent_variant",
        _ => "unsupported",
    };
    json!({
        "class": class,
        "eligible": eligible,
        "confidence_tier": confidence,
        "matched_reference": matched.as_ref().map(|m| m.reference.as_str()),
        "edit_distance": matched.as_ref().map(|m| m.distance),
        "normalized_distance": matched.as_ref().map(|m| m.normalized_distance),
        "source_lineage": opts.lineage,
        "evidence_source": matched.as_ref().map(|m| m.sources.clone()).unwrap_or_default(),
        "reference_lineages": matched.as_ref().map(|m| m.reference_lineages.clone()).unwrap_or_default(),
        "reference_evidence": matched.as_ref().map(|m| m.reference_evidence.clone()).unwrap_or_default(),
        "explanation": explanation,
    })
}

/// Reclassify a report atomically, preserving all original repeat evidence.
pub fn apply(report: &mut Value, opts: &PriorOptions) -> Result<(), String> {
    opts.validate()?;
    let catalog = catalog()?;
    let object = report
        .as_object()
        .ok_or("repeat report must be a JSON object")?;
    let candidates = object
        .get("candidates")
        .and_then(Value::as_array)
        .ok_or("repeat report requires a candidates array")?;
    if !object.get("summary").is_some_and(Value::is_object) {
        return Err("repeat report requires a summary object".into());
    }
    let rejected: &[Value] = match object.get("rejected_candidates") {
        Some(value) => value
            .as_array()
            .ok_or("rejected_candidates must be an array")?
            .as_slice(),
        None => &[],
    };
    let mut rows = Vec::with_capacity(candidates.len() + rejected.len());
    let mut families: HashMap<String, &Value> = HashMap::new();
    let mut ranks = HashSet::new();
    for (index, candidate) in candidates.iter().chain(rejected).enumerate() {
        let row_object = candidate
            .as_object()
            .ok_or_else(|| format!("candidate {} is not an object", index + 1))?;
        let raw_motif = row_object
            .get("motif")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("candidate {} lacks a string motif", index + 1))?;
        let motif = normalized_motif(raw_motif)?;
        let family = canonical(motif.as_bytes()).0;
        if let Some(previous) = families.get(&family) {
            if !same_evidence(previous, candidate) {
                return Err(format!("conflicting duplicate candidate family {family}"));
            }
            continue;
        }
        families.insert(family.clone(), candidate);
        let rank = match candidate
            .get("telomere_prior")
            .and_then(|prior| prior.get("discovery_rank"))
        {
            Some(value) => value
                .as_u64()
                .ok_or("telomere_prior discovery_rank must be an unsigned integer")?,
            None => index as u64,
        };
        if !ranks.insert(rank) {
            return Err(format!("duplicate discovery rank {rank} in repeat report"));
        }
        let mut prior = annotation(&family, &catalog, opts);
        prior["discovery_rank"] = json!(rank);
        let mut evidence = candidate.clone();
        evidence["telomere_prior"] = prior;
        rows.push(evidence);
    }
    // Restore discovery ordering before stable tier ordering. In particular,
    // this retains the scanner's terminal-evidence ordering within each tier.
    rows.sort_by_key(|row| row["telomere_prior"]["discovery_rank"].as_u64().unwrap());
    if opts.policy != "canonical" {
        rows.sort_by_key(prior_tier);
    }
    let raw_count = rows.len();
    let eligible_count = rows
        .iter()
        .filter(|row| row["telomere_prior"]["eligible"] == true)
        .count();
    let (rows, rejected_rows): (Vec<Value>, Vec<Value>) = if opts.policy != "discovery" {
        rows.into_iter()
            .partition(|row| row["telomere_prior"]["eligible"] == true)
    } else {
        let rejected = rows
            .iter()
            .filter(|row| row["telomere_prior"]["eligible"] != true)
            .cloned()
            .collect();
        (rows, rejected)
    };
    let reported_count = rows.len();
    report["schema_version"] = json!(2);
    report["candidates"] = Value::Array(rows);
    report["rejected_candidates"] = Value::Array(rejected_rows);
    report["prior_options"] = serde_json::to_value(opts).map_err(|err| err.to_string())?;
    report["prior_catalog"] = json!({
        "name": "curated_exact_telomere_motifs",
        "format": "embedded TSV: motif, lineage, evidence, source",
        "exact_match_only": true,
        "approximate_global_reference": if opts.policy == "canonical" { None } else { Some("TTAGGG") },
        "distance": "minimum Levenshtein distance over rotations of primitive motifs and both relative strands",
        "normalization": "edit_distance / max(primitive candidate length, primitive reference length)",
        "lineage_policy": if opts.policy == "canonical" { "Explicit expectations annotate agreement; other established families remain eligible. No species identity is inferred." } else { "Legacy policy: explicit related motifs override global eligibility." },
        "confidence_policy": if opts.policy == "canonical" { "Exact established references only; observation ordering preserved. Novel candidates require independent evidence assessment." } else { "Legacy categorical reference support; edit-distance tiers are not calibrated probabilities." },
        "entries": catalog,
    });
    let summary = report["summary"]
        .as_object_mut()
        .expect("summary validated before edits");
    summary.insert("reported_candidates".into(), json!(reported_count));
    summary.insert("raw_reported_candidates".into(), json!(raw_count));
    summary.insert("eligible_telomere_candidates".into(), json!(eligible_count));
    summary.insert(
        "rejected_repeat_candidates".into(),
        json!(raw_count - eligible_count),
    );
    summary.insert(
        "call".into(),
        json!(if eligible_count == 0 {
            "no_call"
        } else {
            "candidates_only"
        }),
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_default_preserves_observation_order_and_other_clade_families() {
        let mut value = report(&["TTAGAC", "TTAGGG", "TTAGGC", "ACGACACC"]);
        let options = PriorOptions {
            related_motifs: vec!["TTAGGC".into()],
            lineage: Some("Nematoda clade expectation".into()),
            ..PriorOptions::default()
        };
        apply(&mut value, &options).unwrap();
        let rows = value["candidates"].as_array().unwrap();
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0]["motif"], "TTAGAC");
        assert_eq!(rows[0]["telomere_prior"]["matched_reference"], "TTAGAC");
        assert_eq!(rows[0]["telomere_prior"]["expected_family_match"], false);
        assert_eq!(rows[2]["telomere_prior"]["expected_family_match"], true);
        assert_eq!(value["rejected_candidates"][0]["motif"], "ACGACACC");
        assert_eq!(value["summary"]["raw_reported_candidates"], 4);
    }

    #[test]
    fn canonical_does_not_promote_edit_neighbours_or_catalog_candidates() {
        // Keep the weak-evidence fixture explicit when the curated catalog changes.
        // TTATTGGG is now an established assembly-terminal reference, tested below.
        let entries = catalog().unwrap();
        let weak = entries
            .iter()
            .find(|entry| entry.motif == "TTTATTAGGG")
            .unwrap();
        assert_eq!(weak.evidence, "comparative");
        assert!(!established_evidence(&weak.evidence));
        let mut value = report(&["TTAGAT", "TTTATTAGGG"]);
        apply(&mut value, &PriorOptions::default()).unwrap();
        assert!(value["candidates"].as_array().unwrap().is_empty());
        let rejected = value["rejected_candidates"].as_array().unwrap();
        assert_eq!(rejected.len(), 2);
        assert!(rejected
            .iter()
            .all(|row| row["telomere_prior"]["matched_reference"].is_null()));
        assert_eq!(value["summary"]["raw_reported_candidates"], 2);
    }

    #[test]
    fn canonical_accepts_curated_assembly_reference_without_confirming_terminal_location() {
        let mut value = report(&["TTATTGGG"]);
        apply(&mut value, &PriorOptions::default()).unwrap();
        let rows = value["candidates"].as_array().unwrap();
        assert_eq!(rows.len(), 1);
        assert!(value["rejected_candidates"].as_array().unwrap().is_empty());
        assert_eq!(rows[0]["telomere_prior"]["class"], "exact_known");
        assert_eq!(rows[0]["telomere_prior"]["matched_reference"], "TTATTGGG");
        assert_eq!(
            rows[0]["telomere_prior"]["reference_evidence"],
            json!(["assembly_terminal"])
        );
        assert_eq!(rows[0]["reads"], 10);
        assert_eq!(rows[0]["repeat_bases"], 1000);
        assert_eq!(rows[0]["max_tract"], 100);
        assert_eq!(rows[0]["status"], "short_read_candidate");
        assert_eq!(value["summary"]["call"], "candidates_only");
    }

    fn legacy_options() -> PriorOptions {
        PriorOptions {
            policy: "telomere".into(),
            ..PriorOptions::default()
        }
    }

    fn report(motifs: &[&str]) -> Value {
        json!({
            "schema_version": 1,
            "summary": {"records": 1000, "bases": 150000, "reported_candidates": motifs.len(), "call": "candidates_only"},
            "candidates": motifs.iter().enumerate().map(|(i, motif)| json!({
                "motif": motif, "reads": 10, "repeat_bases": 1000, "max_tract": 100,
                "status": "short_read_candidate", "rank_score": 1000 - i,
            })).collect::<Vec<_>>()
        })
    }

    #[test]
    fn rejects_user_negative_examples_by_default() {
        let mut value = report(&["AAAAAAAG", "AAT", "ACGCGGCG", "ACAT", "AAAT"]);
        apply(&mut value, &legacy_options()).unwrap();
        assert_eq!(value["candidates"].as_array().unwrap().len(), 0);
        assert_eq!(value["rejected_candidates"].as_array().unwrap().len(), 5);
        assert_eq!(value["summary"]["call"], "no_call");
        assert!(value["rejected_candidates"]
            .as_array()
            .unwrap()
            .iter()
            .all(|row| row["telomere_prior"]["class"] == "unsupported_repeat"));
    }

    #[test]
    fn accepts_expected_variants_and_documented_allium_exception() {
        let mut value = report(&[
            "TTAGGG",
            "TTTAGGG",
            "TTCAGGG",
            "TTAGGC",
            "TTAGAC",
            "TTAGG",
            "CTCGGTTATGGG",
        ]);
        apply(&mut value, &legacy_options()).unwrap();
        assert_eq!(value["summary"]["eligible_telomere_candidates"], 7);
        assert_eq!(value["summary"]["records"], 1000);
        assert_eq!(value["summary"]["bases"], 150000);
        let rows = value["candidates"].as_array().unwrap();
        let allium = rows
            .iter()
            .find(|row| row["motif"] == "CTCGGTTATGGG")
            .unwrap();
        assert_eq!(allium["telomere_prior"]["class"], "exact_known");
        assert!(!allium["telomere_prior"]["evidence_source"]
            .as_array()
            .unwrap()
            .is_empty());
        assert!(rows
            .iter()
            .all(|row| row["telomere_prior"]["source_lineage"].is_null()));
    }

    #[test]
    fn primitive_rotation_and_reverse_complement_invariance() {
        let expected = canonical(b"TTAGGG").0;
        for sequence in ["TTAGGG", "AGGGTT", "CCCTAA", "TTAGGGTTAGGG", "ccctaaccctaa"] {
            let normalized = normalized_motif(sequence).unwrap();
            assert_eq!(canonical(normalized.as_bytes()).0, expected);
            let mut value = report(&[sequence]);
            apply(&mut value, &legacy_options()).unwrap();
            assert_eq!(value["candidates"][0]["telomere_prior"]["edit_distance"], 0);
        }
    }

    #[test]
    fn two_edit_variants_have_lower_support_and_respect_both_cutoffs() {
        let mut value = report(&["TTAGAT", "TTCAGGG", "TTAGGG"]);
        apply(&mut value, &legacy_options()).unwrap();
        let rows = value["candidates"].as_array().unwrap();
        assert_eq!(rows[0]["motif"], "TTAGGG");
        assert_eq!(rows[1]["motif"], "TTCAGGG");
        assert_eq!(rows[2]["motif"], "TTAGAT");
        assert_eq!(rows[2]["telomere_prior"]["edit_distance"], 2);
        assert_eq!(
            rows[2]["telomere_prior"]["confidence_tier"],
            "divergent_variant"
        );
        assert_eq!(
            rows[0]["telomere_prior"]["reference_evidence"],
            json!(["experimental"])
        );
        for opts in [
            PriorOptions {
                max_edits: 1,
                ..legacy_options()
            },
            PriorOptions {
                max_distance: 0.30,
                ..legacy_options()
            },
        ] {
            let mut strict = report(&["TTAGAT"]);
            apply(&mut strict, &opts).unwrap();
            assert_eq!(strict["summary"]["eligible_telomere_candidates"], 0);
        }
    }

    #[test]
    fn global_divergent_exceptions_are_exact_only_unless_context_is_supplied() {
        let mut value = report(&["ATCGGTTATGGG"]);
        apply(&mut value, &legacy_options()).unwrap();
        assert_eq!(value["summary"]["eligible_telomere_candidates"], 0);
        let opts = PriorOptions {
            related_motifs: vec!["CTCGGTTATGGG".into()],
            lineage: Some("Allium".into()),
            ..legacy_options()
        };
        apply(&mut value, &opts).unwrap();
        assert_eq!(
            value["candidates"][0]["telomere_prior"]["class"],
            "related_variant"
        );
        assert_eq!(value["candidates"][0]["telomere_prior"]["edit_distance"], 1);
    }

    #[test]
    fn explicit_related_context_accepts_distant_reference_and_flags_conflict() {
        let opts = PriorOptions {
            related_motifs: vec!["ACGCGGCG".into()],
            lineage: Some("user-supplied sister taxon".into()),
            ..legacy_options()
        };
        let mut value = report(&["TTAGGG", "ACGCGGCG"]);
        apply(&mut value, &opts).unwrap();
        assert_eq!(value["candidates"][0]["motif"], "ACGCGGCG");
        assert_eq!(
            value["candidates"][0]["telomere_prior"]["class"],
            "related_exact"
        );
        assert_eq!(
            value["candidates"][0]["telomere_prior"]["source_lineage"],
            "user-supplied sister taxon"
        );
        assert_eq!(
            value["rejected_candidates"][0]["telomere_prior"]["class"],
            "lineage_conflict"
        );
    }

    #[test]
    fn discovery_retains_repeats_but_does_not_turn_them_into_a_telomere_call() {
        let mut value = report(&["AAT", "AAAAAAAG"]);
        let opts = PriorOptions {
            policy: "discovery".into(),
            ..legacy_options()
        };
        apply(&mut value, &opts).unwrap();
        assert_eq!(value["candidates"].as_array().unwrap().len(), 2);
        assert_eq!(value["summary"]["eligible_telomere_candidates"], 0);
        assert_eq!(value["summary"]["call"], "no_call");
        let before = value.clone();
        apply(&mut value, &opts).unwrap();
        assert_eq!(value, before);
    }

    #[test]
    fn policy_switching_is_lossless_and_idempotent() {
        let mut value = report(&["AAAAAAAG", "TTAGAC", "TTAGGG", "AAT"]);
        let original = value.clone();
        let opts = legacy_options();
        apply(&mut value, &opts).unwrap();
        assert_eq!(value["candidates"][0]["motif"], "TTAGAC");
        assert_eq!(value["summary"]["raw_reported_candidates"], 4);
        let classified = value.clone();
        apply(
            &mut value,
            &PriorOptions {
                policy: "discovery".into(),
                ..opts.clone()
            },
        )
        .unwrap();
        apply(&mut value, &opts).unwrap();
        assert_eq!(value, classified);
        for row in value["candidates"]
            .as_array()
            .unwrap()
            .iter()
            .chain(value["rejected_candidates"].as_array().unwrap())
        {
            let mut evidence = row.clone();
            evidence.as_object_mut().unwrap().remove("telomere_prior");
            assert!(original["candidates"]
                .as_array()
                .unwrap()
                .contains(&evidence));
        }
    }

    #[test]
    fn malformed_report_is_rejected_without_partial_mutation() {
        for motif in ["", "NNAGGG", "TT AGGG", "TTAGGG\n", "🍄"] {
            let mut value = report(&[motif]);
            let before = value.clone();
            assert!(apply(&mut value, &legacy_options()).is_err());
            assert_eq!(value, before);
        }
        let mut value = report(&["TTAGGG", "TTAGGG"]);
        assert!(apply(&mut value, &legacy_options()).is_err());
        let opts = PriorOptions {
            lineage: Some("plants".into()),
            ..legacy_options()
        };
        assert!(opts.validate().is_err());
    }

    fn full_distance(a: &[u8], b: &[u8]) -> usize {
        let mut row: Vec<usize> = (0..=b.len()).collect();
        for (i, left) in a.iter().enumerate() {
            let mut diagonal = row[0];
            row[0] = i + 1;
            for (j, right) in b.iter().enumerate() {
                let old = row[j + 1];
                row[j + 1] = (row[j] + 1)
                    .min(old + 1)
                    .min(diagonal + usize::from(left != right));
                diagonal = old;
            }
        }
        row[b.len()]
    }

    #[test]
    fn bounded_circular_distance_matches_exhaustive_indel_and_strand_reference() {
        let motifs = [
            "TTAGGG", "TTAGAC", "TTCAGGG", "ATG", "ACGT", "AAAT", "ACTGAC", "GGATACG",
        ];
        for left in motifs {
            for right in motifs {
                let mut expected = usize::MAX;
                for reference in [
                    right.as_bytes().to_vec(),
                    reverse_complement(right.as_bytes()),
                ] {
                    for i in 0..left.len() {
                        let a =
                            [left.as_bytes()[i..].to_vec(), left.as_bytes()[..i].to_vec()].concat();
                        for j in 0..reference.len() {
                            let b = [reference[j..].to_vec(), reference[..j].to_vec()].concat();
                            expected = expected.min(full_distance(&a, &b));
                        }
                    }
                }
                for limit in 0..=3 {
                    assert_eq!(
                        circular_distance(left.as_bytes(), right.as_bytes(), limit),
                        (expected <= limit).then_some(expected),
                        "{left} vs {right} limit {limit}"
                    );
                }
            }
        }
    }
}
