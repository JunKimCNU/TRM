"""Packed flank clustering must preserve the original character-level rule."""
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from trm_assessment import flank_clusters


def character_oracle(sequences, representatives=None):
    """Independent graph oracle for the frozen <=2 / >=30 / +/-16 rule."""
    seqs = sorted(set(sequences))
    reps = representatives or {s: s for s in seqs}
    adjacent = {s: set() for s in seqs}
    for i, left in enumerate(seqs):
        for right in seqs[:i]:
            linked = sum(a != b for a, b in zip(left, right)) <= 2
            a, b = reps[left], reps[right]
            if not linked:
                for offset in range(-16, 17):
                    x, y = max(offset, 0), max(-offset, 0)
                    length = min(len(a)-x, len(b)-y)
                    if length >= 30 and sum(a[x+k] != b[y+k] for k in range(length)) <= 2:
                        linked = True
                        break
            if linked:
                adjacent[left].add(right)
                adjacent[right].add(left)
    groups, seen = [], set()
    for sequence in seqs:
        if sequence in seen:
            continue
        pending, group = [sequence], []
        while pending:
            item = pending.pop()
            if item in seen:
                continue
            seen.add(item)
            group.append(item)
            pending.extend(adjacent[item] - seen)
        groups.append(sorted(group))
    return groups


class PackedFlankClusterTests(unittest.TestCase):
    def compare(self, values):
        representatives = {s[:20]: s for s in values}
        anchors = list(representatives)
        expected = character_oracle(anchors, representatives)
        self.assertEqual(flank_clusters(anchors, representatives), expected)
        self.assertEqual(flank_clusters(anchors[::-1] + anchors, representatives), expected)
        return expected

    def test_empty_singleton_and_duplicate_order(self):
        self.assertEqual(flank_clusters([]), [])
        self.assertEqual(flank_clusters(['A'*20]*3), [['A'*20]])
        self.compare(['T'*80, 'A'*80, 'C'*20])

    def test_all_substitution_codes_count_bases_not_bits(self):
        for first in 'ACGT':
            for second in 'ACGT':
                if first == second:
                    continue
                for mismatch in [1, 2, 3]:
                    with self.subTest(first=first, second=second, mismatch=mismatch):
                        groups = self.compare([first*20, second*mismatch+first*(20-mismatch)])
                        self.assertEqual(len(groups), 1 if mismatch <= 2 else 2)

    def test_transitive_components_preserve_sorted_representatives(self):
        values = ['A'*20, 'CC'+'A'*18, 'CCCC'+'A'*16]
        self.assertEqual(self.compare(values), [sorted(values)])

    def test_positive_and_negative_offset_limit(self):
        rng = random.Random(20261004)
        source = ''.join(rng.choices('ACGT', k=110))
        for offset in [1, 7, 15, 16, 17]:
            for values in [[source[:64], source[offset:offset+64]],
                           [source[offset:offset+64], source[:64]]]:
                with self.subTest(offset=offset, reversed=values[0] != source[:64]):
                    self.assertEqual(len(self.compare(values)), 1 if offset <= 16 else 2)

    def test_overlap_minimum_and_short_anchor_fallback(self):
        source = ''.join(random.Random(887).choices('ACGT', k=80))
        self.assertEqual(len(self.compare([source[:46], source[16:46]])), 1)
        self.assertEqual(len(self.compare([source[:45], source[16:45]])), 2)
        self.assertEqual(len(self.compare(['A'*20, 'CC'+'A'*18])), 1)

    def test_overlap_mismatch_threshold_and_masked_outside_bases(self):
        source = ''.join(random.Random(320).choices('ACGT', k=96))
        for mismatches in [0, 1, 2, 3]:
            # Changes at the far end are inside the 73-base aligned overlap,
            # while the unaligned leading 7 bases of the first rep are ignored.
            second = list(source[7:80])
            for position in range(mismatches):
                second[-position-1] = {'A':'T','T':'A','G':'C','C':'G'}[second[-position-1]]
            groups = self.compare([source[:80], ''.join(second)])
            self.assertEqual(len(groups), 1 if mismatches <= 2 else 2)

    def test_random_lengths_and_related_representatives(self):
        rng = random.Random(802604)
        for case in range(180):
            values = []
            source = ''.join(rng.choices('ACGT', k=120))
            for _ in range(rng.randrange(2, 18)):
                length = rng.randrange(20, 81)
                if rng.random() < .7:
                    offset = rng.randrange(0, 24)
                    value = list(source[offset:offset+length])
                    for pos in rng.sample(range(length), k=rng.randrange(0, 5)):
                        value[pos] = rng.choice('ACGT')
                    values.append(''.join(value))
                else:
                    values.append(''.join(rng.choices('ACGT', k=length)))
            with self.subTest(case=case):
                self.compare(values)

    def test_sketch_capacity_128_random_and_shifted_sets(self):
        rng = random.Random(128)
        self.compare([''.join(rng.choices('ACGT', k=80)) for _ in range(128)])
        contexts = [''.join(rng.choices('ACGT', k=120)) for _ in range(8)]
        self.compare([context[offset:offset+80] for context in contexts for offset in range(16)])

    def test_invalid_anchor_and_representative_validation_preserved(self):
        for value in ['A'*19, 'A'*21, 'A'*19+'N', 'a'*20]:
            with self.subTest(anchor=value), self.assertRaisesRegex(ValueError, 'exactly 20 ACGT'):
                flank_clusters([value])
        for rep in ['', 'A'*19, 'A'*81, 'C'+'A'*20, 'A'*20+'N']:
            with self.subTest(representative=rep), self.assertRaisesRegex(ValueError, '20..80 ACGT'):
                flank_clusters(['A'*20], {'A'*20:rep})
        with self.assertRaisesRegex(ValueError, '20..80 ACGT'):
            flank_clusters(['A'*20], {'C'*20:'C'*20})
        self.assertEqual(flank_clusters(['A'*20], {}), [['A'*20]])


if __name__ == '__main__':
    unittest.main()
