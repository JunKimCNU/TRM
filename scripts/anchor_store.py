"""Disk-backed bookkeeping for irregular-anchor validation; no biological model."""
import gzip
import json
import resource
import sqlite3


def apply_memory_limit(mib):
    """Limit this process's virtual address space, retaining stricter OS limits."""
    if mib is not None and (type(mib) is not int or mib <= 0):
        raise ValueError('memory limit must be a positive integer MiB')
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    if mib is not None:
        limits = [mib * 1024 * 1024]
        limits.extend(x for x in (soft, hard) if x != resource.RLIM_INFINITY)
        soft = min(limits)
        resource.setrlimit(resource.RLIMIT_AS, (soft, hard))
    return {'requested_memory_limit_mib': mib,
            'effective_address_space_limit_bytes': None if soft == resource.RLIM_INFINITY else soft,
            'memory_limit_kind': 'process_virtual_address_space_RLIMIT_AS',
            'hard_limit_unchanged': True}


class AlignmentGroup(list):
    """Selection representatives, with the original number of reported rows.

    select_alignment needs primary multiplicity (zero, one or >=two) and the
    greatest nonprimary AS. Rows are never deduplicated. With one primary, the
    exact primary object survives both iterations used by select_alignment.
    """
    def __init__(self, rows, reported_count):
        super().__init__(rows)
        self.reported_count = reported_count

    def __len__(self):
        return self.reported_count


class AnchorLengths:
    def __init__(self, store):
        self.store = store

    def get(self, name, default=None):
        value = self.store.length(name)
        return default if value is None else value

    def __len__(self):
        return self.store.db.execute('SELECT count(*) FROM anchors').fetchone()[0]


class AnchorStore:
    def __init__(self, path, cache_mib=8):
        if type(cache_mib) is not int or cache_mib <= 0:
            raise ValueError('index cache must be a positive integer MiB')
        self.db = sqlite3.connect(path)
        self.db.execute('PRAGMA cache_size = ' + str(-cache_mib * 1024))
        self.db.execute('PRAGMA mmap_size = 0')
        self.db.execute('PRAGMA temp_store = FILE')
        # Scratch can be reconstructed; no rollback copy is needed.
        self.db.execute('PRAGMA journal_mode = OFF')
        self.db.executescript('''
          CREATE TABLE anchors (id TEXT PRIMARY KEY, length TEXT NOT NULL, used INTEGER DEFAULT 0);
          CREATE TABLE arrays (id TEXT PRIMARY KEY);
          CREATE TABLE records (source TEXT, ordinal TEXT, name TEXT, length TEXT,
            joint INTEGER DEFAULT 0, long_joint INTEGER DEFAULT 0, nuclear INTEGER DEFAULT 0,
            PRIMARY KEY(source, ordinal));
          CREATE TABLE alignments (id TEXT PRIMARY KEY, count INTEGER, first_primary TEXT,
            second_primary TEXT, best_other TEXT);
        ''')
        self.lengths = AnchorLengths(self)

    def close(self):
        self.db.close()

    def add_anchor(self, name, length):
        if not name or length == 0:
            raise ValueError('empty FASTA record')
        try:
            self.db.execute('INSERT INTO anchors(id,length) VALUES (?,?)', (name, str(length)))
        except sqlite3.IntegrityError as exc:
            raise ValueError('duplicate FASTA identifier: ' + name) from exc

    def load_fasta(self, path):
        opener = gzip.open if str(path).endswith('.gz') else open
        name, length = None, 0
        with opener(path, 'rt', encoding='ascii') as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if name is not None:
                        self.add_anchor(name, length)
                    words = line[1:].split()
                    if not words:
                        raise ValueError('empty FASTA identifier')
                    name, length = words[0], 0
                else:
                    if name is None:
                        raise ValueError('sequence before FASTA identifier')
                    length += len(line)
        if name is not None:
            self.add_anchor(name, length)

    def length(self, name):
        row = self.db.execute('SELECT length FROM anchors WHERE id=?', (name,)).fetchone()
        return int(row[0]) if row else None

    def mark_anchor(self, name):
        changed = self.db.execute('UPDATE anchors SET used=1 WHERE id=? AND used=0', (name,)).rowcount
        if changed != 1:
            raise ValueError('duplicate or missing exported anchor ID')

    def add_array(self, row):
        try:
            self.db.execute('INSERT INTO arrays VALUES (?)', (row['array_id'],))
        except sqlite3.IntegrityError as exc:
            raise ValueError('duplicate array/anchor ID') from exc
        key = (str(row['source_index']), str(row['read_ordinal']))
        identity = (row['read_id'], str(row['read_length']))
        old = self.db.execute('SELECT name,length FROM records WHERE source=? AND ordinal=?', key).fetchone()
        if old is not None and old != identity:
            raise ValueError('inconsistent name/length for the same source record')
        if old is None:
            self.db.execute('INSERT INTO records(source,ordinal,name,length) VALUES (?,?,?,?)', key + identity)

    def mark_outcome(self, row):
        self.db.execute('''UPDATE records SET joint=max(joint,?), long_joint=max(long_joint,?),
            nuclear=max(nuclear,?) WHERE source=? AND ordinal=?''',
            (int(row['joint_terminal_context']), int(row['joint_long_model_span_context']),
             int(row['joint_nuclear_terminal_context']), str(row['source_index']), str(row['read_ordinal'])))

    def record_counts(self):
        return self.db.execute('SELECT count(*),coalesce(sum(joint),0),coalesce(sum(long_joint),0),coalesce(sum(nuclear),0) FROM records').fetchone()

    def finish_export(self, expectedarrays, expectedanchors):
        if self.db.execute('SELECT count(*) FROM anchors WHERE used=0').fetchone()[0]:
            raise ValueError('FASTA and exported anchor IDs differ')
        if (self.db.execute('SELECT count(*) FROM arrays').fetchone()[0] != expectedarrays or
                len(self.lengths) != expectedanchors):
            raise ValueError('array/anchor summary counts disagree with actual export')

    def add_alignment(self, row):
        key = row['query']
        previous = self.db.execute('SELECT count,first_primary,second_primary,best_other FROM alignments WHERE id=?', (key,)).fetchone()
        count, first, second, other = previous or (0, None, None, None)
        encoded = json.dumps(row, ensure_ascii=False)
        if row['alignment_type'] == 'P':
            if first is None:
                first = encoded
            elif second is None:
                second = encoded
        elif other is None or row['alignment_score'] > json.loads(other)['alignment_score']:
            other = encoded
        self.db.execute('INSERT OR REPLACE INTO alignments VALUES (?,?,?,?,?)',
                        (key, count + 1, first, second, other))

    def get(self, name, default=None):
        row = self.db.execute('SELECT count,first_primary,second_primary,best_other FROM alignments WHERE id=?', (name,)).fetchone()
        if row is None:
            return default
        count, first, second, other = row
        return AlignmentGroup([json.loads(x) for x in (first, second, other) if x is not None], count)
