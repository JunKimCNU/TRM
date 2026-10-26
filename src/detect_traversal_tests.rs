//! Independent traversal regressions. The oracle visits maximal exact runs by
//! lag equality, then uses the unchanged production extension without its cache.
use super::*;
use std::collections::BTreeSet;

type Span = (String, Option<bool>, usize, usize);

fn exact_run_oracle(seq: &[u8], cfg: &Config) -> HashMap<String, Vec<Tract>> {
    let mut hits = HashMap::new();
    for k in cfg.min_k..=cfg.max_k.min(seq.len() / cfg.min_copies) {
        let required = cfg.min_copies.max(cfg.seed_bases.div_ceil(k)) * k;
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
                // Ignore the returned extended endpoint: every maximal exact
                // run is independently inspected, including covered ones.
                record_seed(seq, cfg, start, unit, start..p, &mut hits);
            }
        }
    }
    hits
}

fn noncontained_spans(hits: &HashMap<String, Vec<Tract>>) -> BTreeSet<Span> {
    let mut spans = BTreeSet::new();
    for (family, tracts) in hits {
        for t in tracts {
            if !tracts.iter().any(|other| {
                other.strand == t.strand
                    && other.start <= t.start
                    && other.end >= t.end
                    && (other.start < t.start || other.end > t.end)
            }) {
                spans.insert((family.clone(), t.strand, t.start, t.end));
            }
        }
    }
    spans
}

fn assert_matches_oracle(seq: &[u8], cfg: &Config, context: &str) {
    let actual = discover(seq, cfg);
    let expected = exact_run_oracle(seq, cfg);
    assert_eq!(
        noncontained_spans(&actual),
        noncontained_spans(&expected),
        "{context}; k={}, copies={}, seed_bases={}; sequence={}",
        cfg.min_k,
        cfg.min_copies,
        cfg.seed_bases,
        String::from_utf8_lossy(seq)
    );
}

#[test]
fn cross_family_extension_does_not_hide_exact_seeds_at_any_phase() {
    let first = b"GGTGTTGG";
    let second = b"GGTGCTGG";
    assert_ne!(canonical(first).0, canonical(second).0);
    for copies in [2, 3] {
        for shift in 0..first.len() {
            for reverse in [false, true] {
                let mut seq = vec![b'N'; shift];
                seq.extend(first.repeat(8));
                seq.extend(second.repeat(10));
                seq.push(b'N');
                if reverse {
                    seq = reverse_complement(&seq);
                }
                let cfg = Config {
                    min_k: first.len(),
                    max_k: first.len(),
                    min_copies: copies,
                    seed_bases: 1,
                    ..Config::default()
                };
                assert_matches_oracle(&seq, &cfg, "neighbor family under extension");
                let hits = discover(&seq, &cfg);
                assert!(hits.contains_key(&canonical(first).0));
                assert!(hits.contains_key(&canonical(second).0));
                for motif in [first, second] {
                    let family = canonical(motif).0;
                    let target = TargetIndex::new(vec![family.clone()], &cfg).discover(&seq, &cfg);
                    let actual = HashMap::from([(family.clone(), hits[&family].clone())]);
                    assert_eq!(noncontained_spans(&actual), noncontained_spans(&target));
                }
            }
        }
    }
}

#[test]
fn exact_run_oracle_covers_supported_periods_phases_and_two_copy_path() {
    for k in 3..=64 {
        let mut unit: Vec<u8> = (0..k)
            .map(|i| if i % 2 == 0 { b'A' } else { b'C' })
            .collect();
        unit[k - 1] = b'G';
        for copies in [2, 3] {
            for seed_bases in [1, 24] {
                for shift in [0, 1, k - 1] {
                    for reverse in [false, true] {
                        let cfg = Config {
                            min_k: k,
                            max_k: k,
                            min_copies: copies,
                            seed_bases,
                            ..Config::default()
                        };
                        let needed = copies.max(seed_bases.div_ceil(k));
                        let mut seq = vec![b'N'; shift];
                        seq.extend(unit.repeat(needed));
                        seq.push(b'N');
                        let mut changed = unit.clone();
                        changed[0] = b'T';
                        seq.extend(changed.repeat(needed + 3));
                        seq.extend(unit[1..].iter());
                        seq.extend(unit.repeat(needed + 2));
                        seq.push(b'N');
                        if reverse {
                            seq = reverse_complement(&seq);
                        }
                        assert_matches_oracle(&seq, &cfg, "supported period/phase/strand");
                    }
                }
            }
        }
    }
}

fn random(state: &mut u64) -> usize {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state as usize
}

#[test]
fn imperfect_mosaics_match_uncached_exact_run_oracle() {
    let mut rng = 20260928_u64;
    for case in 0..600 {
        let k = [3, 4, 6, 7, 8, 12, 16, 31, 32, 64][random(&mut rng) % 10];
        let unit: Vec<u8> = (0..k).map(|_| b"ACGT"[random(&mut rng) % 4]).collect();
        if primitive_len(&unit) != k {
            continue;
        }
        let mut seq = vec![b'N'; random(&mut rng) % k];
        for _ in 0..2 + random(&mut rng) % 6 {
            let mut rotated = unit.clone();
            rotated.rotate_left(random(&mut rng) % k);
            if random(&mut rng) % 2 == 0 {
                rotated = reverse_complement(&rotated);
            }
            if random(&mut rng) % 4 == 0 {
                rotated[random(&mut rng) % k] = b"ACGT"[random(&mut rng) % 4];
            }
            let mut chunk = rotated.repeat(3 + random(&mut rng) % 18);
            let step = 25 + random(&mut rng) % 45;
            for p in (20 + random(&mut rng) % 30..chunk.len()).step_by(step) {
                chunk[p] = b"ACGTN"[random(&mut rng) % 5];
            }
            seq.extend(chunk);
            if random(&mut rng) % 5 == 0 {
                seq.extend(vec![b'N'; 1 + random(&mut rng) % 7]);
            }
        }
        let cfg = Config {
            min_k: k,
            max_k: k,
            min_copies: 2 + random(&mut rng) % 2,
            seed_bases: [1, 24][random(&mut rng) % 2],
            min_identity: [0.70, 0.85, 0.95, 1.0][random(&mut rng) % 4],
            ..Config::default()
        };
        assert_matches_oracle(&seq, &cfg, &format!("random mosaic {case}"));
    }
}

#[test]
fn sparse_errors_and_phase_switches_match_uncached_extensions() {
    for unit in [
        b"TTAGGG".as_slice(),
        b"TTTAGGG",
        b"GGTGTTGG",
        b"CACCGTATATACGGTT",
    ] {
        let k = unit.len();
        for copies in [2, 3] {
            let mut seq = unit.repeat(150);
            for p in (55..seq.len()).step_by(37) {
                seq[p] = if seq[p] == b'A' { b'C' } else { b'A' };
            }
            seq.extend(&unit[1..]);
            seq.extend(unit.repeat(12));
            seq.extend(reverse_complement(unit).repeat(12));
            seq.push(b'N');
            seq.extend(unit.repeat(6));
            let cfg = Config {
                min_k: k,
                max_k: k,
                min_copies: copies,
                seed_bases: 1,
                ..Config::default()
            };
            assert_matches_oracle(&seq, &cfg, "sparse errors and changed phase/strand");
        }
    }
}
