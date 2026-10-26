//! Repeat-family representation under rotation and reverse complementation.
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

pub fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

fn min_rotation(s: &[u8]) -> Vec<u8> {
    // Booth's algorithm: linear time, including long repeat units.
    let n = s.len();
    let (mut i, mut j, mut k) = (0, 1, 0);
    while i < n && j < n && k < n {
        let a = s[(i + k) % n];
        let b = s[(j + k) % n];
        if a == b {
            k += 1;
            continue;
        }
        if a > b {
            i += k + 1;
            if i <= j {
                i = j + 1;
            }
        } else {
            j += k + 1;
            if j <= i {
                j = i + 1;
            }
        }
        k = 0;
    }
    let start = i.min(j);
    (0..n).map(|x| s[(start + x) % n]).collect()
}

/// Returns the canonical family and strand relative to that family.
/// None strand denotes reverse-complement symmetry (orientation uninformative).
pub fn canonical(s: &[u8]) -> (String, Option<bool>) {
    let f = min_rotation(s);
    let r = min_rotation(&reverse_complement(s));
    let strand = if f == r { None } else { Some(f < r) };
    let key = if f <= r { f } else { r };
    (String::from_utf8(key).expect("ACGT motif"), strand)
}

pub fn primitive_len(s: &[u8]) -> usize {
    for k in 1..=s.len() / 2 {
        if s.len() % k == 0 && s.iter().enumerate().all(|(i, &b)| b == s[i % k]) {
            return k;
        }
    }
    s.len()
}

pub fn is_dna(s: &[u8]) -> bool {
    s.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn canonical_invariant_and_strand() {
        let key = canonical(b"TTAGGG").0;
        assert_eq!(key, "AACCCT");
        for s in [b"TAGGGT", b"AGGGTT", b"CCCTAA", b"CCTAAC"] {
            assert_eq!(canonical(s).0, key);
        }
        assert_ne!(canonical(b"TTAGGG").1, canonical(b"CCCTAA").1);
        assert_eq!(canonical(b"AT").1, None);
        assert_eq!(primitive_len(b"TTAGGGTTAGGG"), 6);
        assert_eq!(primitive_len(b"AAAAAA"), 1);
    }
    #[test]
    fn booth_matches_exhaustive_rotation() {
        for len in 1..=7 {
            for mut bits in 0..4usize.pow(len as u32) {
                let mut s = vec![b'A'; len];
                for b in &mut s {
                    *b = b"ACGT"[bits % 4];
                    bits /= 4;
                }
                let expected = (0..len)
                    .map(|i| [s[i..].to_vec(), s[..i].to_vec()].concat())
                    .min()
                    .unwrap();
                assert_eq!(min_rotation(&s), expected);
            }
        }
    }
}
