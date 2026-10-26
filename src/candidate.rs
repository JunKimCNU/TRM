//! Shared operations on saved candidate evidence and derived annotations.
use serde_json::Value;

/// Compare original observations exactly, ignoring only derived annotations.
/// Callers validate both candidate objects before comparing them.
pub(crate) fn same_evidence(left: &Value, right: &Value) -> bool {
    let ignored = ["telomere_prior", "support", "candidate_eligible"];
    left.as_object()
        .expect("candidate object already validated")
        .iter()
        .filter(|(key, _)| !ignored.contains(&key.as_str()))
        .eq(right
            .as_object()
            .expect("candidate object already validated")
            .iter()
            .filter(|(key, _)| !ignored.contains(&key.as_str())))
}

/// Legacy prior ordering, shared by prior and background reclassification.
/// Canonical policy preserves observation ordering instead of using this key.
pub(crate) fn prior_tier(row: &Value) -> u8 {
    match row["telomere_prior"]["class"].as_str() {
        Some("related_exact" | "exact_known") => 0,
        Some("related_variant" | "ttaggg_variant") => {
            if row["telomere_prior"]["edit_distance"] == 1 {
                1
            } else {
                2
            }
        }
        _ => 3,
    }
}
