use crate::motif::{canonical, is_dna, primitive_len, reverse_complement};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, serde::Serialize)]
pub struct Config {
    pub min_k: usize,
    pub max_k: usize,
    pub seed_bases: usize,
    pub min_copies: usize,
    pub min_tract: usize,
    pub min_identity: f64,
    pub end_window: usize,
    pub min_reads: u64,
    pub max_candidates: usize,
    pub mode: String,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            min_k: 3,
            max_k: 32,
            seed_bases: 24,
            min_copies: 3,
            min_tract: 500,
            min_identity: 0.85,
            end_window: 50,
            min_reads: 2,
            max_candidates: 100_000,
            mode: "reads".into(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Tract {
    pub start: usize,
    pub end: usize,
    pub identity: f64,
    pub strand: Option<bool>,
}

/// Extend a seed with a substitution-only X-drop alignment to its periodic motif.
/// Ambiguous bases are hard boundaries. Keeping the highest-scoring prefix avoids
/// adding trailing random sequence. This is not an indel-aware alignment.
fn extend(
    seq: &[u8],
    anchor: usize,
    motif: &[u8],
    start: usize,
    end: usize,
    identity: f64,
    family: &(String, Option<bool>),
) -> (usize, usize) {
    let k = motif.len();
    let penalty = if identity >= 1.0 {
        i64::MAX / 4
    } else {
        (identity / (1.0 - identity)).ceil() as i64
    };
    let drop = (2 * k).max(12) as i64;
    let expected = |p: usize| motif[(p + k - anchor % k) % k];
    let opposite_seed = |pair: &[u8]| {
        let unit = &pair[..k];
        if unit != &pair[k..] || !is_dna(unit) {
            return false;
        }
        let other = canonical(unit);
        other.0 == family.0 && other.1.is_some() && other.1 != family.1
    };
    let mut left = start;
    let (mut score, mut best) = (0i64, 0i64);
    for p in (0..start).rev() {
        if !is_dna(&seq[p..p + 1]) {
            break;
        }
        if seq[p] != expected(p)
            && family.1.is_some()
            && p + 1 >= 2 * k
            && opposite_seed(&seq[p + 1 - 2 * k..p + 1])
        {
            break;
        }
        score += if seq[p] == expected(p) { 1 } else { -penalty };
        if score > best {
            best = score;
            left = p;
        }
        if best - score > drop {
            break;
        }
    }
    let mut right = end;
    score = 0;
    best = 0;
    for (p, &base) in seq.iter().enumerate().skip(end) {
        if !matches!(base, b'A' | b'C' | b'G' | b'T') {
            break;
        }
        if base != expected(p)
            && family.1.is_some()
            && p + 2 * k <= seq.len()
            && opposite_seed(&seq[p..p + 2 * k])
        {
            break;
        }
        score += if base == expected(p) { 1 } else { -penalty };
        if score > best {
            best = score;
            right = p + 1;
        }
        if best - score > drop {
            break;
        }
    }
    (left, right)
}

/// Block comparisons cost O(n * sum(1/k)) on nonrepetitive reads; comparisons
/// short-circuit and canonicalization runs only after finding an exact seed.
/// Approximate extension is additional work; a bounded accepted-tract cache
/// avoids re-extending contained seeds without skipping other motif templates.
pub fn discover(seq: &[u8], cfg: &Config) -> HashMap<String, Vec<Tract>> {
    let mut hits: HashMap<String, Vec<Tract>> = HashMap::new();
    for k in cfg.min_k..=cfg.max_k.min(seq.len() / cfg.min_copies) {
        let seed_copies = cfg.min_copies.max(cfg.seed_bases.div_ceil(k));
        let required = seed_copies * k;
        let mut covered = Vec::new();
        if seed_copies == 2 {
            // Two-copy seeds can straddle the block grid without containing
            // two complete aligned blocks. Scan lag equality for this case.
            let mut p = k;
            while p < seq.len() {
                if seq[p] != seq[p - k] {
                    p += 1;
                    continue;
                }
                let start = p - k;
                while p < seq.len() && seq[p] == seq[p - k] {
                    p += 1;
                }
                let unit = &seq[start..start + k];
                if p - start >= required && is_dna(unit) && primitive_len(unit) == k {
                    record_discovery_seed(seq, cfg, start, unit, start..p, &mut hits, &mut covered);
                }
            }
            continue;
        }
        let mut pos = 0;
        // Every >=3-copy seed contains two complete blocks on this grid,
        // even when the seed's actual start is not a multiple of k.
        while pos + 2 * k <= seq.len() {
            let unit = &seq[pos..pos + k];
            if unit != &seq[pos + k..pos + 2 * k] {
                pos += k;
                continue;
            }
            let mut end = pos + 2 * k;
            while end + k <= seq.len() && unit == &seq[end..end + k] {
                end += k;
            }
            if !is_dna(unit) || primitive_len(unit) != k {
                pos = end - k;
                continue;
            }
            let mut start = pos;
            while start > 0 && seq[start - 1] == unit[(start - 1 + k - pos % k) % k] {
                start -= 1;
            }
            while end < seq.len() && seq[end] == unit[(end + k - pos % k) % k] {
                end += 1;
            }
            if end - start < required {
                // Keep the next probe on the block grid after a partial end.
                pos = (end / k).max(pos / k + 1) * k;
                continue;
            }
            record_discovery_seed(seq, cfg, pos, unit, start..end, &mut hits, &mut covered);
            // Different period-k units can share up to k-1 boundary bases.
            // Retain the partial final block so a minimum-length neighboring
            // seed is not skipped, while guaranteeing forward progress. Only
            // skip the exact run: an approximate extension can contain valid
            // seeds of a different family of this same period.
            pos = (end / k).max(pos / k + 1) * k;
        }
    }
    deduplicate(&mut hits);
    hits
}

// Local to one period in one record. A template records expected bases at
// absolute coordinates modulo k, so phase and strand are both distinguished.
// Eviction can only cause recomputation; it never hides a seed or joins tracts.
struct CoveredSeed {
    template: u128,
    start: usize,
    end: usize,
}

fn record_discovery_seed(
    seq: &[u8],
    cfg: &Config,
    anchor: usize,
    unit: &[u8],
    seed: std::ops::Range<usize>,
    hits: &mut HashMap<String, Vec<Tract>>,
    covered: &mut Vec<CoveredSeed>,
) {
    let k = unit.len();
    let mut template = 0u128;
    for position in 0..k {
        template =
            (template << 2) | u128::from(base_bits(unit[(position + k - anchor % k) % k]).unwrap());
    }
    if covered
        .iter()
        .any(|c| c.template == template && c.start <= seed.start && seed.end <= c.end)
    {
        return;
    }
    let seed_start = seed.start;
    let result = record_seed_result(seq, cfg, anchor, unit, seed, hits);
    if result.accepted {
        covered.retain(|c| c.end > seed_start && c.template != template);
        const MAX_CACHED_TEMPLATES: usize = 64;
        if covered.len() == MAX_CACHED_TEMPLATES {
            covered.remove(0);
        }
        covered.push(CoveredSeed {
            template,
            start: result.left,
            end: result.right,
        });
    }
}

fn deduplicate(hits: &mut HashMap<String, Vec<Tract>>) {
    // Remove duplicate / overlapping seeds, preserving a real observed span.
    // Never join two short tracts into apparent >=500 bp continuity.
    for tracts in hits.values_mut() {
        tracts.sort_by_key(|t| (t.start, std::cmp::Reverse(t.end)));
        let mut keep: Vec<Tract> = Vec::new();
        for t in tracts.drain(..) {
            if keep
                .last()
                .is_some_and(|last| t.end <= last.end && t.strand == last.strand)
            {
                continue;
            }
            keep.push(t);
        }
        *tracts = keep;
    }
}

#[derive(Clone, Debug)]
struct TargetSeed {
    family_index: usize,
    unit: Vec<u8>,
    required: usize,
}

/// Fixed-window two-bit seed lookup over the entire read. A small exact bitmap
/// rejects most nonmatching windows before consulting the sparse hash table.
/// Only explicitly requested motifs are searched; this is not de novo discovery.
#[derive(Clone, Debug)]
pub struct TargetIndex {
    pub families: Vec<String>,
    probe_bases: usize,
    mask: u64,
    quick: Vec<u64>,
    seeds: HashMap<u64, Vec<TargetSeed>>,
}

fn base_bits(base: u8) -> Option<u64> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl TargetIndex {
    pub fn new(families: Vec<String>, cfg: &Config) -> Self {
        let required = |length: usize| cfg.min_copies.max(cfg.seed_bases.div_ceil(length)) * length;
        let probe_bases = families
            .iter()
            .map(|family| required(family.len()))
            .min()
            .unwrap_or(12)
            .min(12);
        let mut index = Self {
            families,
            probe_bases,
            mask: (1u64 << (2 * probe_bases)) - 1,
            quick: vec![0; 1024],
            seeds: HashMap::new(),
        };
        for (family_index, family) in index.families.iter().enumerate() {
            let forward = family.as_bytes();
            let reverse = reverse_complement(forward);
            let mut seen = HashSet::new();
            for sequence in [forward, reverse.as_slice()] {
                for offset in 0..sequence.len() {
                    let unit: Vec<u8> = (0..sequence.len())
                        .map(|i| sequence[(i + offset) % sequence.len()])
                        .collect();
                    if !seen.insert(unit.clone()) {
                        continue;
                    }
                    let mut key = 0u64;
                    for position in 0..probe_bases {
                        key = (key << 2) | base_bits(unit[position % unit.len()]).unwrap();
                    }
                    let quick_key = (key & 0xffff) as usize;
                    index.quick[quick_key / 64] |= 1u64 << (quick_key % 64);
                    index.seeds.entry(key).or_default().push(TargetSeed {
                        family_index,
                        required: required(unit.len()),
                        unit,
                    });
                }
            }
        }
        index
    }

    pub fn discover(&self, seq: &[u8], cfg: &Config) -> HashMap<String, Vec<Tract>> {
        let mut hits = HashMap::new();
        // Most WGS reads have no target seed. Allocate family bookkeeping only
        // after verifying one, rather than zeroing 256 counters per input read.
        let mut covered_until = Vec::<usize>::new();
        let (mut key, mut valid) = (0u64, 0usize);
        for (position, &base) in seq.iter().enumerate() {
            let Some(bits) = base_bits(base) else {
                key = 0;
                valid = 0;
                continue;
            };
            key = ((key << 2) | bits) & self.mask;
            valid += 1;
            if valid < self.probe_bases {
                continue;
            }
            let quick_key = (key & 0xffff) as usize;
            if self.quick[quick_key / 64] & (1u64 << (quick_key % 64)) == 0 {
                continue;
            }
            let Some(seeds) = self.seeds.get(&key) else {
                continue;
            };
            let anchor = position + 1 - self.probe_bases;
            for seed in seeds {
                let k = seed.unit.len();
                if anchor
                    < covered_until
                        .get(seed.family_index)
                        .copied()
                        .unwrap_or(0)
                        .saturating_sub(k - 1)
                    || seq.len() - anchor < seed.required
                    || !(0..seed.required)
                        .all(|offset| seq[anchor + offset] == seed.unit[offset % k])
                {
                    continue;
                }
                if covered_until.is_empty() {
                    covered_until.resize(self.families.len(), 0);
                }
                let mut left = anchor;
                let mut right = anchor + seed.required;
                while left > 0 && seq[left - 1] == seed.unit[(left - 1 + k - anchor % k) % k] {
                    left -= 1;
                }
                while right < seq.len() && seq[right] == seed.unit[(right - anchor) % k] {
                    right += 1;
                }
                covered_until[seed.family_index] = covered_until[seed.family_index].max(
                    record_seed(seq, cfg, anchor, &seed.unit, left..right, &mut hits),
                );
            }
        }
        deduplicate(&mut hits);
        hits
    }
}

fn record_seed(
    seq: &[u8],
    cfg: &Config,
    anchor: usize,
    unit: &[u8],
    seed: std::ops::Range<usize>,
    hits: &mut HashMap<String, Vec<Tract>>,
) -> usize {
    record_seed_result(seq, cfg, anchor, unit, seed, hits).right
}

struct SeedResult {
    left: usize,
    right: usize,
    accepted: bool,
}

fn record_seed_result(
    seq: &[u8],
    cfg: &Config,
    anchor: usize,
    unit: &[u8],
    seed: std::ops::Range<usize>,
    hits: &mut HashMap<String, Vec<Tract>>,
) -> SeedResult {
    let k = unit.len();
    let family = canonical(unit);
    let (left, right) = extend(
        seq,
        anchor,
        unit,
        seed.start,
        seed.end,
        cfg.min_identity,
        &family,
    );
    let matched = (left..right)
        .filter(|&p| seq[p] == unit[(p + k - anchor % k) % k])
        .count();
    let identity = matched as f64 / (right - left) as f64;
    let accepted = identity + 1e-12 >= cfg.min_identity;
    if accepted {
        let (key, strand) = family;
        hits.entry(key).or_default().push(Tract {
            start: left,
            end: right,
            identity,
            strand,
        });
    }
    SeedResult {
        left,
        right,
        accepted,
    }
}

#[cfg(test)]
#[path = "detect_traversal_tests.rs"]
mod traversal_tests;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn targeted_seed_search_all_phases_and_strands() {
        let cfg = Config::default();
        for motif in ["TTAGGG", "TTAGGC", "CTCGGTTATGGG", "CCCAAGGTTTAAAAGG"] {
            let family = canonical(motif.as_bytes()).0;
            let index = TargetIndex::new(vec![family.clone()], &cfg);
            for unit in [
                motif.as_bytes().to_vec(),
                reverse_complement(motif.as_bytes()),
            ] {
                for shift in 0..unit.len() * 2 {
                    let copies = cfg.min_copies.max(cfg.seed_bases.div_ceil(unit.len()));
                    let mut sequence = vec![b'N'; 151 + shift];
                    sequence.extend(unit.repeat(copies));
                    sequence.extend([b'N'; 20]);
                    let targeted = index.discover(&sequence, &cfg);
                    let exhaustive = discover(&sequence, &cfg);
                    let tract = &targeted[&family][0];
                    assert_eq!(
                        (tract.start, tract.end),
                        (151 + shift, 151 + shift + unit.len() * copies)
                    );
                    assert_eq!(tract.strand, exhaustive[&family][0].strand);
                    assert_eq!(targeted.len(), 1);
                }
            }
        }
    }

    #[test]
    fn targeted_multiple_tracts_do_not_bridge_or_double_count() {
        let cfg = Config::default();
        let family = canonical(b"TTAGGG").0;
        let index = TargetIndex::new(vec![family.clone()], &cfg);
        let sequence = [b"TTAGGG".repeat(50), vec![b'N'; 23], b"CCCTAA".repeat(50)].concat();
        let hits = index.discover(&sequence, &cfg);
        assert_eq!(hits[&family].len(), 2);
        assert_eq!(hits[&family][0].end - hits[&family][0].start, 300);
        assert_eq!(hits[&family][1].end - hits[&family][1].start, 300);
        assert_ne!(hits[&family][0].strand, hits[&family][1].strand);
    }
    #[test]
    fn all_phases_and_primitive_motifs() {
        let cfg = Config::default();
        for shift in 0..12 {
            let mut s = vec![b'N'; shift];
            s.extend(b"CTCGGTTATGGG".repeat(50));
            s.extend([b'N'; 17]);
            let hits = discover(&s, &cfg);
            let key = canonical(b"CTCGGTTATGGG").0;
            let t = &hits[&key][0];
            assert_eq!((t.start, t.end), (shift, shift + 600));
            assert_eq!(hits.len(), 1);
        }
    }
    #[test]
    fn minimum_length_seeds_at_every_phase() {
        for k in 3..=64 {
            // A unique G guarantees that the unit has primitive period k.
            let mut unit: Vec<u8> = (0..k)
                .map(|i| if i % 2 == 0 { b'A' } else { b'C' })
                .collect();
            unit[k - 1] = b'G';
            for copies in [2, 3] {
                let cfg = Config {
                    min_k: k,
                    max_k: k,
                    min_copies: copies,
                    seed_bases: 1,
                    ..Config::default()
                };
                let key = canonical(&unit).0;
                for shift in 0..k {
                    let mut seq = vec![b'N'; shift];
                    seq.extend(unit.repeat(copies));
                    seq.extend([b'N'; 2]);
                    let hits = discover(&seq, &cfg);
                    let tract = &hits
                        .get(&key)
                        .unwrap_or_else(|| panic!("missed k={k}, copies={copies}, shift={shift}"))
                        [0];
                    assert_eq!(
                        (tract.start, tract.end),
                        (shift, shift + k * copies),
                        "k={k}, copies={copies}, shift={shift}"
                    );
                    assert_eq!(hits.len(), 1);
                }
            }
        }
    }
    #[test]
    fn default_minimum_seed_is_phase_independent() {
        for shift in 0..6 {
            let mut seq = vec![b'N'; shift];
            seq.extend(b"TTAGGG".repeat(4));
            seq.push(b'N');
            let hits = discover(&seq, &Config::default());
            assert_eq!(hits[&canonical(b"TTAGGG").0][0].end - shift, 24);
        }
    }
    #[test]
    fn overlapping_minimum_seeds_remain_visible() {
        let cfg = Config {
            min_k: 3,
            max_k: 3,
            seed_bases: 1,
            min_identity: 1.0,
            ..Config::default()
        };
        // The adjacent units share their boundary CG, but form two different
        // exact three-copy arrays. The first tract ends between grid probes.
        let hits = discover(b"NACGACGACGCCGCCGCN", &cfg);
        assert_eq!(hits.len(), 2);
        let first = &hits[&canonical(b"ACG").0][0];
        let second = &hits[&canonical(b"CGC").0][0];
        assert_eq!((first.start, first.end), (1, 10));
        assert_eq!((second.start, second.end), (8, 17));
    }
    #[test]
    fn no_homopolymer_harmonics_or_ambiguous_seeds() {
        for s in [b"A".repeat(600), b"AT".repeat(300), b"N".repeat(600)] {
            assert!(discover(&s, &Config::default()).is_empty());
        }
    }
    #[test]
    fn do_not_bridge_ambiguity_or_strand_switch() {
        let mut s = b"TTAGGG".repeat(50);
        s.extend([b'N'; 12]);
        s.extend(b"TTAGGG".repeat(50));
        let hits = discover(&s, &Config::default());
        assert!(hits.values().flatten().all(|t| t.end - t.start <= 300));
        let s = [b"TTAGGG".repeat(50), b"CCCTAA".repeat(50)].concat();
        let hits = discover(&s, &Config::default());
        assert!(hits.values().flatten().all(|t| t.end - t.start < 500));
    }
    #[test]
    fn substitution_extension() {
        let mut s = b"TTAGGG".repeat(200);
        for p in (55..s.len()).step_by(37) {
            s[p] = if s[p] == b'A' { b'C' } else { b'A' };
        }
        let hits = discover(&s, &Config::default());
        let key = canonical(b"TTAGGG").0;
        assert!(hits[&key]
            .iter()
            .any(|t| t.end - t.start > 1000 && t.identity > 0.96));
    }
    #[test]
    fn near_symmetric_reverse_seed_stops_extension() {
        let motif = b"CACCGTATATACGGTT";
        let reverse = crate::motif::reverse_complement(motif);
        let seq = [motif.repeat(20), reverse.repeat(20)].concat();
        let cfg = Config {
            min_k: motif.len(),
            max_k: motif.len(),
            ..Config::default()
        };
        let hits = discover(&seq, &cfg);
        let tracts = &hits[&canonical(motif).0];
        assert!(tracts.iter().all(|t| t.end - t.start < 500));
        assert!(tracts.iter().any(|t| t.strand == Some(true)));
        assert!(tracts.iter().any(|t| t.strand == Some(false)));
    }
    #[test]
    fn near_symmetric_reverse_seeds_at_all_periods() {
        let mut random = 13_u64;
        for k in 3..=64 {
            for _ in 0..4 {
                let half: Vec<u8> = (0..k / 2)
                    .map(|_| {
                        random ^= random << 13;
                        random ^= random >> 7;
                        random ^= random << 17;
                        b"ACGT"[(random & 3) as usize]
                    })
                    .collect();
                let mut unit = half.clone();
                if k % 2 == 1 {
                    // Odd units cannot equal their own reverse complement;
                    // a central base yields the closest possible pair.
                    unit.push(b'A');
                }
                unit.extend(crate::motif::reverse_complement(&half));
                if k % 2 == 0 {
                    unit[0] = if unit[0] == b'A' { b'C' } else { b'A' };
                }
                if primitive_len(&unit) != k || canonical(&unit).1.is_none() {
                    continue;
                }
                let reverse = crate::motif::reverse_complement(&unit);
                let copies = 320_usize.div_ceil(k).max(3);
                for first_reverse in [false, true] {
                    let seq = if first_reverse {
                        [reverse.repeat(copies), unit.repeat(copies)].concat()
                    } else {
                        [unit.repeat(copies), reverse.repeat(copies)].concat()
                    };
                    let cfg = Config {
                        min_k: k,
                        max_k: k,
                        ..Config::default()
                    };
                    let hits = discover(&seq, &cfg);
                    let tracts = &hits[&canonical(&unit).0];
                    assert!(
                        tracts.iter().all(|t| t.end - t.start < 500),
                        "bridged strand switch at k={k}"
                    );
                    assert!(
                        tracts.iter().any(|t| t.strand == Some(true)),
                        "missing + at k={k}"
                    );
                    assert!(
                        tracts.iter().any(|t| t.strand == Some(false)),
                        "missing - at k={k}"
                    );
                }
            }
        }
    }
}
