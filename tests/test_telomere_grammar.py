"""Exact grammar path contracts checked by exhaustive numeric path enumeration."""
from collections import defaultdict
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import telomere_grammar as grammar


def enumerate_paths(edges):
    outgoing = defaultdict(list)
    for edge in edges:
        outgoing[edge[0]].append(edge)
    paths = defaultdict(list)

    def visit(origin, end, chain):
        paths[origin, end].append(chain)
        for edge in outgoing[end]:
            visit(origin, edge[1], chain + [edge])

    for edge in edges:
        visit(edge[0], edge[1], [edge])
    return paths


class ExactGrammarPaths(unittest.TestCase):
    def check_graph(self, n, edges):
        enumerated = enumerate_paths(edges)
        actual, endpoints = grammar.solve_paths(n, edges)
        expected = []
        for end in sorted({b for a, b in enumerated}):
            start = min(a for a, b in enumerated if b == end)
            paths = enumerated[start, end]
            expected.append({'start0': start, 'end0': end, 'span_bp': end-start,
                'minimum_block_count': min(map(len, paths)), 'maximum_block_count': max(map(len, paths)),
                'path_count_capped2': min(2, len(paths))})
        self.assertEqual(endpoints, expected)
        maximal = [r for r in expected if not any(
            q['start0'] <= r['start0'] and q['end0'] >= r['end0'] and q != r for q in expected)]
        maximal.sort(key=lambda r: (r['start0'], -r['end0']))
        self.assertEqual([{k: v for k, v in row.items() if k != 'representative_tokens'} for row in actual], maximal)
        for row in actual:
            trace = [tuple(edge) for edge in row['representative_tokens']]
            self.assertIn(trace, enumerated[row['start0'], row['end0']])
            self.assertEqual(len(trace), row['minimum_block_count'])
        self.assertEqual(grammar.solve_paths(n, list(reversed(edges))), (actual, endpoints))

    def test_all_1024_small_graphs(self):
        universe = [(a, b, (a+b) % 3) for a in range(4) for b in range(a+1, 5)]
        for mask in range(1 << len(universe)):
            self.check_graph(4, [edge for i, edge in enumerate(universe) if mask & (1 << i)])

    def test_distinct_labels_and_block_count_ambiguity(self):
        self.check_graph(6, [(0, 2, 0), (0, 2, 1), (2, 4, 2), (4, 6, 3), (0, 6, 4)])

    def test_overlap_is_not_a_gapless_path(self):
        actual, _ = grammar.solve_paths(8, [(0, 3, 0), (2, 5, 1), (5, 8, 2)])
        self.assertEqual([(r['start0'], r['end0']) for r in actual], [(0, 3), (2, 8)])
        self.check_graph(8, [(0, 3, 0), (2, 5, 1), (5, 8, 2)])

    def test_invalid_graphs_fail(self):
        for n, edges in [(True, []), (-1, []), (3, [(0, 0, 0)]), (3, [(0, 4, 0)]),
                         (3, [(0, 1, -1)]), (3, [(0, 1, 0), (0, 1, 0)]), (3, [(0, 1)])]:
            with self.subTest(n=n, edges=edges), self.assertRaises(ValueError):
                grammar.solve_paths(n, edges)


if __name__ == '__main__':
    unittest.main()
