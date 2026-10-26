"""Shared bounded I/O for saved reference and cross-period evidence.

Only file stability, strict JSON and canonical content seals belong here.
Each caller retains its error type, error prefix, size policy, schema and
biological validation. No data is read or validated on import.
"""
import hashlib
import json
import os
from pathlib import Path
import stat


class EvidenceIO:
    """Keep identical I/O behavior while preserving caller-specific errors."""

    def __init__(self, error_type, prefix):
        self.error_type = error_type
        self.prefix = prefix

    def fail(self, message):
        raise self.error_type(self.prefix + message)

    def unique_object(self, pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                self.fail('duplicate JSON key ' + repr(key))
            result[key] = value
        return result

    def reject_constant(self, value):
        self.fail('non-finite JSON number ' + value)

    def payload_sha256(self, document):
        """Seal canonical UTF-8 JSON excluding only the top-level seal."""
        if not isinstance(document, dict):
            self.fail('document must be an object')
        try:
            raw = json.dumps(
                {key: value for key, value in document.items() if key != 'payload_sha256'},
                sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False,
            ).encode('utf-8')
        except (TypeError, ValueError, UnicodeError) as error:
            raise self.error_type(self.prefix + 'non-JSON payload') from error
        return hashlib.sha256(raw).hexdigest()

    def read_bytes(self, path, *, max_bytes):
        """Bound reads and compare path/open-descriptor identities throughout."""
        path = Path(path).absolute()
        identity = lambda info: (
            info.st_dev, info.st_ino, info.st_size,
            info.st_mtime_ns, info.st_ctime_ns, info.st_mode,
        )
        try:
            resolved = path.resolve(strict=True)
            before = path.stat()
            if not stat.S_ISREG(before.st_mode) or before.st_size > max_bytes:
                self.fail('source must be a small regular file: ' + str(path))
            with path.open('rb') as handle:
                opened = os.fstat(handle.fileno())
                data = handle.read(max_bytes + 1)
                closed = os.fstat(handle.fileno())
            after = path.stat()
            if (len(data) > max_bytes or identity(before) != identity(opened)
                    or identity(opened) != identity(closed) or identity(closed) != identity(after)
                    or resolved != path.resolve(strict=True)):
                self.fail('source changed while reading: ' + str(path))
        except (OSError, RuntimeError) as error:
            raise self.error_type(self.prefix + 'cannot read ' + str(path)) from error
        return data

    def decode_json(self, data, path, *, object_pairs_hook, parse_constant):
        """Decode with the caller's hooks and retain the original byte hash."""
        try:
            document = json.loads(
                data, object_pairs_hook=object_pairs_hook, parse_constant=parse_constant,
            )
        except (ValueError, UnicodeError) as error:
            raise self.error_type(self.prefix + 'malformed JSON in ' + str(path)) from error
        if not isinstance(document, dict):
            self.fail('document must be an object')
        return document, hashlib.sha256(data).hexdigest()
