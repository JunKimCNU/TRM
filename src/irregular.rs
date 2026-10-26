//! Optional, seed-conditional per-unit edit model; no biological classification.
//!
//! Each complete supplied unit consumes its literal sequence or one S/I/D edit
//! (edits require unit length >=5). There are no free between-unit gap edges.
//! Every edge is locally bounded, so endpoint feasibility has no path-history
//! constraint: earliest start dominates later starts for maximal query coverage.
//! Among paths at that start we minimize total edits, keeping one deterministic
//! representative and flagging equal-cost alternatives. Counts describe that
//! representative, not guaranteed component counts across all possible parses.
//! Adjacent edited units can create adjacent edits; this is not a global 1-bp
//! gap bound. Query span, exact-unit bases and alignment edits remain separate.

use crate::motif::{canonical, is_dna, primitive_len, reverse_complement};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Definition {
    pub id: String,
    pub units: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Config {
    pub min_span: usize,
    pub min_exact_anchors: usize,
    pub max_nodes: usize,
    pub max_record_bases: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            min_span: 100,
            min_exact_anchors: 4,
            max_nodes: 2_000_000,
            max_record_bases: 25_000_000,
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Edge {
    pub start: usize,
    pub end: usize,
    pub cost: usize,
    pub variant: usize,
    pub ambiguous: bool,
}

#[derive(Clone, Debug)]
pub struct Path {
    pub start: usize,
    pub end: usize,
    pub cost: usize,
    pub ambiguous: bool,
    pub edges: Vec<usize>,
}

#[derive(Clone)]
struct Node {
    start: usize,
    end: usize,
    cost: usize,
    ambiguous: bool,
    edge: usize,
    previous: Option<usize>,
}

#[derive(Clone, Copy)]
struct State {
    node: usize,
    extended: bool,
}

/// Numeric DAG kernel. Edges must be sorted by nondecreasing query end.
/// All edges are allowed independently of their predecessor/path history.
pub fn solve_tokens(
    n: usize,
    edges: &[Edge],
    minimum: usize,
    limit: usize,
) -> Result<Vec<Path>, String> {
    if edges.is_empty() {
        return Ok(Vec::new());
    }
    let mut last = 0;
    let mut max_len = 0;
    for edge in edges {
        if edge.start >= edge.end || edge.end > n || edge.end < last {
            return Err("invalid or unsorted token coordinates".into());
        }
        max_len = max_len.max(edge.end - edge.start);
        last = edge.end;
    }
    let mut ring: Vec<Option<State>> = vec![None; max_len + 1];
    let mut nodes: Vec<Node> = Vec::new();
    let mut completed = Vec::new();
    let mut next = 0;
    for end in 1..=n {
        let slot = end % ring.len();
        if let Some(old) = ring[slot].take() {
            if !old.extended {
                completed.push(old.node);
            }
        }
        let mut best: Option<Node> = None;
        while next < edges.len() && edges[next].end == end {
            let edge = edges[next];
            let previous_slot = edge.start % ring.len();
            let previous = ring[previous_slot]
                .as_mut()
                .filter(|state| nodes[state.node].end == edge.start);
            let mut candidate = Node {
                start: edge.start,
                end,
                cost: edge.cost,
                ambiguous: edge.ambiguous,
                edge: next,
                previous: None,
            };
            if let Some(state) = previous {
                state.extended = true;
                let node = &nodes[state.node];
                candidate.start = node.start;
                candidate.cost = candidate
                    .cost
                    .checked_add(node.cost)
                    .ok_or("edit count overflow")?;
                candidate.ambiguous |= node.ambiguous;
                candidate.previous = Some(state.node);
            }
            match &mut best {
                None => best = Some(candidate),
                Some(current) => {
                    let a = (candidate.start, candidate.cost);
                    let b = (current.start, current.cost);
                    match a.cmp(&b) {
                        std::cmp::Ordering::Less => best = Some(candidate),
                        std::cmp::Ordering::Equal => current.ambiguous = true,
                        std::cmp::Ordering::Greater => {}
                    }
                }
            }
            next += 1;
        }
        if let Some(best) = best {
            if nodes.len() >= limit {
                return Err("per-orientation token node limit exceeded; no complete result".into());
            }
            let node = nodes.len();
            nodes.push(best);
            ring[slot] = Some(State {
                node,
                extended: false,
            });
        }
    }
    for state in ring.into_iter().flatten() {
        if !state.extended {
            completed.push(state.node);
        }
    }
    completed.sort_by_key(|&i| {
        (
            nodes[i].start,
            std::cmp::Reverse(nodes[i].end),
            nodes[i].cost,
        )
    });
    let mut out = Vec::new();
    let mut furthest = 0;
    for i in completed {
        let node = &nodes[i];
        if node.end - node.start < minimum || node.end <= furthest {
            continue;
        }
        furthest = node.end;
        let mut trace = Vec::new();
        let mut current = Some(i);
        while let Some(index) = current {
            trace.push(nodes[index].edge);
            current = nodes[index].previous;
        }
        trace.reverse();
        out.push(Path {
            start: node.start,
            end: node.end,
            cost: node.cost,
            ambiguous: node.ambiguous,
            edges: trace,
        });
    }
    Ok(out)
}

#[derive(Clone, Debug)]
struct Variant {
    unit: usize,
    op: char,
    offset: usize,
    edits: usize,
    ambiguous: bool,
}

struct Grammar {
    units: Vec<Vec<u8>>,
    index: HashMap<(usize, u128), usize>,
    variants: Vec<Variant>,
    lengths: Vec<usize>,
    max_len: usize,
}

fn base(b: u8) -> Option<u128> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn code(bytes: &[u8]) -> u128 {
    bytes.iter().fold(0, |value, &b| {
        (value << 2) | base(b).expect("validated DNA")
    })
}

impl Grammar {
    fn new(units: Vec<Vec<u8>>) -> Self {
        let mut value = Self {
            units,
            index: HashMap::new(),
            variants: Vec::new(),
            lengths: Vec::new(),
            max_len: 0,
        };
        for i in 0..value.units.len() {
            let unit = value.units[i].clone();
            value.add(
                &unit,
                Variant {
                    unit: i,
                    op: '=',
                    offset: 0,
                    edits: 0,
                    ambiguous: false,
                },
            );
            if unit.len() < 5 {
                continue;
            }
            for pos in 0..unit.len() {
                for b in b"ACGT" {
                    if *b != unit[pos] {
                        let mut observed = unit.clone();
                        observed[pos] = *b;
                        value.add(
                            &observed,
                            Variant {
                                unit: i,
                                op: 'S',
                                offset: pos,
                                edits: 1,
                                ambiguous: false,
                            },
                        );
                    }
                }
                let mut observed = unit.clone();
                observed.remove(pos);
                value.add(
                    &observed,
                    Variant {
                        unit: i,
                        op: 'D',
                        offset: pos,
                        edits: 1,
                        ambiguous: false,
                    },
                );
            }
            for pos in 0..=unit.len() {
                for b in b"ACGT" {
                    let mut observed = unit.clone();
                    observed.insert(pos, *b);
                    value.add(
                        &observed,
                        Variant {
                            unit: i,
                            op: 'I',
                            offset: pos,
                            edits: 1,
                            ambiguous: false,
                        },
                    );
                }
            }
        }
        value.lengths = value
            .index
            .keys()
            .map(|x| x.0)
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        value.lengths.sort_unstable();
        value.max_len = *value.lengths.last().expect("nonempty grammar");
        value
    }

    fn add(&mut self, observed: &[u8], variant: Variant) {
        let key = (observed.len(), code(observed));
        if let Some(&i) = self.index.get(&key) {
            let old = &mut self.variants[i];
            match variant.edits.cmp(&old.edits) {
                std::cmp::Ordering::Less => *old = variant,
                std::cmp::Ordering::Equal => old.ambiguous = true,
                std::cmp::Ordering::Greater => {}
            }
        } else {
            self.index.insert(key, self.variants.len());
            self.variants.push(variant);
        }
    }

    /// Earliest-finish interval scheduling counts disjoint literal anchors.
    fn anchor_support(&self, seq: &[u8], required: usize) -> usize {
        let mut last_end = 0;
        let mut count = 0;
        for end in 1..=seq.len() {
            if self.units.iter().any(|unit| {
                end >= unit.len()
                    && end - unit.len() >= last_end
                    && seq[end - unit.len()..end] == unit[..]
            }) {
                count += 1;
                last_end = end;
                if count >= required {
                    break;
                }
            }
        }
        count
    }

    fn edges(&self, seq: &[u8], limit: usize) -> Result<Vec<Edge>, String> {
        let mut rolling = 0_u128;
        let mut valid = 0;
        let mask = (1_u128 << (self.max_len * 2)) - 1;
        let mut out = Vec::new();
        for (pos, &b) in seq.iter().enumerate() {
            let Some(value) = base(b) else {
                valid = 0;
                rolling = 0;
                continue;
            };
            rolling = ((rolling << 2) | value) & mask;
            valid += 1;
            for &len in &self.lengths {
                if valid < len {
                    continue;
                }
                let key = (len, rolling & ((1_u128 << (2 * len)) - 1));
                if let Some(&index) = self.index.get(&key) {
                    if out.len() >= limit {
                        return Err(
                            "per-orientation token edge limit exceeded; no complete result".into(),
                        );
                    }
                    let v = &self.variants[index];
                    out.push(Edge {
                        start: pos + 1 - len,
                        end: pos + 1,
                        cost: v.edits,
                        variant: index,
                        ambiguous: v.ambiguous,
                    });
                }
            }
        }
        Ok(out)
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct TokenTrace {
    pub start0: usize,
    pub end0: usize,
    pub unit_index: usize,
    pub edit: char,
    pub edit_offset: usize,
    pub ambiguous: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct ModelEvidence {
    pub matches: usize,
    pub substitutions: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub model_bases: usize,
    pub exact_unit_bases: usize,
    pub identity: f64,
    pub literal_gate_hits: usize,
    pub representative_unit_counts: Vec<usize>,
    pub optimal_parse_ambiguous: bool,
    pub max_consecutive_edited_units: usize,
    pub max_unmatched_query_run: usize,
    pub trace: Vec<TokenTrace>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Array {
    pub group_id: String,
    pub start0: usize,
    pub end0: usize,
    pub span_bp: usize,
    pub strand: String,
    pub alignment_strand: String,
    pub model: ModelEvidence,
}

#[derive(Default, Serialize)]
pub struct ScanResult {
    pub arrays: Vec<Array>,
    pub orientations_considered: usize,
    pub orientations_passing_literal_gate: usize,
}

pub struct Scanner {
    pub definitions: Vec<Definition>,
    pub config: Config,
    grammars: Vec<[Grammar; 2]>,
}

fn orientation_ambiguities(intervals: &[(usize, usize, bool)]) -> Vec<bool> {
    let mut indices = [Vec::new(), Vec::new()];
    for &(start, end, strand) in intervals {
        indices[usize::from(strand)].push((start, end));
    }
    for rows in &mut indices {
        rows.sort_unstable();
        let mut maximum = 0;
        for (_, end) in rows {
            maximum = maximum.max(*end);
            *end = maximum;
        }
    }
    intervals
        .iter()
        .map(|&(start, end, strand)| {
            let opposite = &indices[usize::from(!strand)];
            let last = opposite.partition_point(|&(a, _)| a < end);
            last > 0 && opposite[last - 1].1 > start
        })
        .collect()
}

impl Scanner {
    pub fn new(mut definitions: Vec<Definition>, config: Config) -> Result<Self, String> {
        if definitions.is_empty() || definitions.len() > 8 {
            return Err("require 1..8 definition groups".into());
        }
        if config.min_span == 0
            || config.min_exact_anchors == 0
            || config.max_nodes == 0
            || config.max_record_bases == 0
        {
            return Err("span, anchor and computational limits must be positive".into());
        }
        let mut ids = HashSet::new();
        let mut grammars = Vec::new();
        for definition in &mut definitions {
            if definition.id.trim().is_empty() || !ids.insert(definition.id.clone()) {
                return Err("nonempty unique group IDs required".into());
            }
            if definition.units.is_empty() || definition.units.len() > 8 {
                return Err("require 1..8 units per group".into());
            }
            definition
                .units
                .iter_mut()
                .for_each(|s| *s = s.to_ascii_uppercase());
            definition.units.sort();
            let mut families = HashSet::new();
            for unit in &definition.units {
                if !(3..=32).contains(&unit.len())
                    || !is_dna(unit.as_bytes())
                    || primitive_len(unit.as_bytes()) != unit.len()
                {
                    return Err("units must be primitive 3..32-base ACGT sequences".into());
                }
                if !families.insert(canonical(unit.as_bytes()).0) {
                    return Err("rotation/RC-equivalent duplicate units".into());
                }
            }
            let forward: Vec<Vec<u8>> = definition
                .units
                .iter()
                .map(|x| x.as_bytes().to_vec())
                .collect();
            let reverse = forward.iter().map(|x| reverse_complement(x)).collect();
            grammars.push([Grammar::new(forward), Grammar::new(reverse)]);
        }
        Ok(Self {
            definitions,
            config,
            grammars,
        })
    }

    pub fn scan(&self, seq: &[u8]) -> Result<ScanResult, String> {
        if seq.len() > self.config.max_record_bases {
            return Err(
                "record length exceeds explicit computational limit; no complete result".into(),
            );
        }
        let mut result = ScanResult::default();
        for (definition, grammars) in self.definitions.iter().zip(&self.grammars) {
            let begin = result.arrays.len();
            for (orientation, grammar) in grammars.iter().enumerate() {
                result.orientations_considered += 1;
                if grammar.anchor_support(seq, self.config.min_exact_anchors)
                    < self.config.min_exact_anchors
                {
                    continue;
                }
                result.orientations_passing_literal_gate += 1;
                let edges = grammar.edges(seq, self.config.max_nodes)?;
                let paths = solve_tokens(
                    seq.len(),
                    &edges,
                    self.config.min_span,
                    self.config.max_nodes,
                )?;
                for path in paths {
                    let hits = grammar.anchor_support(&seq[path.start..path.end], usize::MAX);
                    if hits < self.config.min_exact_anchors {
                        continue;
                    }
                    let mut model = ModelEvidence {
                        matches: 0,
                        substitutions: 0,
                        insertions: 0,
                        deletions: 0,
                        model_bases: 0,
                        exact_unit_bases: 0,
                        identity: 0.0,
                        literal_gate_hits: hits,
                        representative_unit_counts: vec![0; grammar.units.len()],
                        optimal_parse_ambiguous: path.ambiguous,
                        max_consecutive_edited_units: 0,
                        max_unmatched_query_run: 0,
                        trace: Vec::new(),
                    };
                    let (mut edited, mut unmatched) = (0, 0);
                    for edge_id in path.edges {
                        let edge = edges[edge_id];
                        let v = &grammar.variants[edge.variant];
                        let k = grammar.units[v.unit].len();
                        model.model_bases += k;
                        model.representative_unit_counts[v.unit] += 1;
                        model.matches += k - usize::from(v.op == 'S' || v.op == 'D');
                        model.substitutions += usize::from(v.op == 'S');
                        model.insertions += usize::from(v.op == 'I');
                        model.deletions += usize::from(v.op == 'D');
                        if v.op == '=' {
                            model.exact_unit_bases += k;
                            edited = 0;
                        } else {
                            edited += 1;
                            model.max_consecutive_edited_units =
                                model.max_consecutive_edited_units.max(edited);
                        }
                        // Query-consuming operations, including consecutive edits across a boundary.
                        for q in 0..edge.end - edge.start {
                            let bad = (v.op == 'I' || v.op == 'S') && q == v.offset;
                            unmatched = if bad { unmatched + 1 } else { 0 };
                            model.max_unmatched_query_run =
                                model.max_unmatched_query_run.max(unmatched);
                        }
                        model.trace.push(TokenTrace {
                            start0: edge.start,
                            end0: edge.end,
                            unit_index: v.unit,
                            edit: v.op,
                            edit_offset: v.offset,
                            ambiguous: v.ambiguous,
                        });
                    }
                    debug_assert_eq!(
                        model.substitutions + model.insertions + model.deletions,
                        path.cost
                    );
                    debug_assert_eq!(
                        model.matches + model.substitutions + model.insertions,
                        path.end - path.start
                    );
                    model.identity = model.matches as f64
                        / (model.matches + model.substitutions + model.insertions + model.deletions)
                            as f64;
                    let strand = if orientation == 0 { "+" } else { "-" }.to_string();
                    result.arrays.push(Array {
                        group_id: definition.id.clone(),
                        start0: path.start,
                        end0: path.end,
                        span_bp: path.end - path.start,
                        strand: strand.clone(),
                        alignment_strand: strand,
                        model,
                    });
                }
            }
            let coordinates: Vec<_> = result.arrays[begin..]
                .iter()
                .map(|x| (x.start0, x.end0, x.alignment_strand == "+"))
                .collect();
            for (array, ambiguous) in result.arrays[begin..]
                .iter_mut()
                .zip(orientation_ambiguities(&coordinates))
            {
                if ambiguous {
                    array.strand = "ambiguous".into();
                }
            }
        }
        result.arrays.sort_by(|a, b| {
            (&a.group_id, a.start0, a.end0, &a.alignment_strand).cmp(&(
                &b.group_id,
                b.start0,
                b.end0,
                &b.alignment_strand,
            ))
        });
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Pure integer graph oracle: no constructed nucleotide sequences.
    fn enumerate(end: usize, edges: &[Edge]) -> Vec<(usize, usize)> {
        let mut paths = Vec::new();
        for e in edges.iter().filter(|e| e.end == end) {
            paths.push((e.start, e.cost));
            for (start, cost) in enumerate(e.start, edges) {
                paths.push((start, cost + e.cost));
            }
        }
        paths
    }

    #[test]
    fn integer_dags_match_exhaustive_path_oracle() {
        let edges: Vec<_> = (1..=5)
            .flat_map(|end| {
                (0..end).map(move |start| Edge {
                    start,
                    end,
                    cost: (start + end) % 2,
                    variant: 0,
                    ambiguous: false,
                })
            })
            .collect();
        for mask in 0u32..(1 << edges.len()) {
            let chosen: Vec<_> = edges
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
                .map(|(_, e)| *e)
                .collect();
            let mut oracle: Vec<_> = (1..=5)
                .filter_map(|end| {
                    enumerate(end, &chosen)
                        .into_iter()
                        .min()
                        .map(|(s, c)| (s, end, c))
                })
                .collect();
            oracle = oracle
                .iter()
                .copied()
                .filter(|&(s, e, _)| {
                    !oracle
                        .iter()
                        .any(|&(a, b, _)| a <= s && b >= e && (a < s || b > e))
                })
                .collect();
            oracle.sort_unstable();
            let actual = solve_tokens(5, &chosen, 1, 100).unwrap();
            let mut got: Vec<_> = actual.iter().map(|p| (p.start, p.end, p.cost)).collect();
            got.sort_unstable();
            assert_eq!(got, oracle, "edge mask {mask}");
            for p in actual {
                let mut pos = p.start;
                let mut cost = 0;
                for id in p.edges {
                    let e = chosen[id];
                    assert_eq!(e.start, pos);
                    pos = e.end;
                    cost += e.cost;
                }
                assert_eq!((pos, cost), (p.end, p.cost));
            }
        }
    }

    #[test]
    fn integer_limits_and_unconnected_intervals() {
        let e = [
            Edge {
                start: 0,
                end: 5,
                cost: 0,
                variant: 0,
                ambiguous: false,
            },
            Edge {
                start: 8,
                end: 13,
                cost: 1,
                variant: 1,
                ambiguous: false,
            },
        ];
        assert_eq!(solve_tokens(13, &e, 5, 100).unwrap().len(), 2);
        assert!(solve_tokens(13, &e, 10, 100).unwrap().is_empty());
        assert!(solve_tokens(13, &e, 1, 1).is_err());
        assert!(solve_tokens(12, &e, 1, 100).is_err());
        assert!(solve_tokens(13, &[e[1], e[0]], 1, 100).is_err());
    }

    #[test]
    fn integer_tied_paths_are_flagged() {
        let e = [
            Edge {
                start: 0,
                end: 3,
                cost: 0,
                variant: 0,
                ambiguous: false,
            },
            Edge {
                start: 0,
                end: 6,
                cost: 1,
                variant: 1,
                ambiguous: false,
            },
            Edge {
                start: 3,
                end: 6,
                cost: 1,
                variant: 2,
                ambiguous: false,
            },
        ];
        let out = solve_tokens(6, &e, 1, 100).unwrap();
        assert_eq!(out.len(), 1);
        assert!(out[0].ambiguous);
        assert_eq!((out[0].start, out[0].end, out[0].cost), (0, 6, 1));
    }

    #[test]
    fn integer_orientation_sweep_matches_pairwise_overlap() {
        let rows: Vec<_> = (0..12)
            .flat_map(|a| {
                ((a + 1)..=12)
                    .flat_map(move |b| [false, true].into_iter().map(move |strand| (a, b, strand)))
            })
            .collect();
        for step in 1..20 {
            let selection: Vec<_> = rows.iter().copied().step_by(step).collect();
            let expected: Vec<_> = selection
                .iter()
                .map(|&(a, b, s)| selection.iter().any(|&(x, y, t)| s != t && a < y && x < b))
                .collect();
            assert_eq!(orientation_ambiguities(&selection), expected);
        }
    }
}
