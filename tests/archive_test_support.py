"""Explicit opt-in for tests tied to preserved local study archives.

The release suite keeps these tests visible as skips. Set
TRM_RUN_ARCHIVE_TESTS=1 only in a workspace retaining the named historical
artifacts. This flag never downloads data or changes the test assertions.
"""
import os
import unittest

ARCHIVE_TESTS_ENABLED = os.environ.get('TRM_RUN_ARCHIVE_TESTS') == '1'
ARCHIVE_TEST_REASON = (
    'Archive-only integration: set TRM_RUN_ARCHIVE_TESTS=1 and retain the '
    'named historical study/binary artifacts; not part of the portable release suite'
)


def archive_only(test):
    return unittest.skipUnless(ARCHIVE_TESTS_ENABLED, ARCHIVE_TEST_REASON)(test)
