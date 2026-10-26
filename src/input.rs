//! Streaming FASTA/FASTQ input, including gzip streams and concatenated members.
//!
//! Sequence bases are normalized to uppercase. IUPAC ambiguity codes are
//! retained so callers can decide how to treat them. IDs contain the first
//! whitespace-delimited token of the header. After an error, discard the reader.

use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};

#[path = "profile.rs"]
pub mod profile;

const BUFFER_SIZE: usize = 256 * 1024;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub id: String,
    pub seq: Vec<u8>,
    /// Only populated by an explicitly quality-retaining reader. FASTA and
    /// records above that reader's retention limit have no retained quality.
    pub quality: Option<Vec<u8>>,
}

/// Open a path, or standard input when `path` is `-`. Compression is detected
/// from bytes rather than the filename. Decoder errors propagate while reading.
pub fn open(path: &str) -> Result<Box<dyn BufRead + Send>, String> {
    let source: Box<dyn Read + Send> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).map_err(|e| format!("cannot open {path:?}: {e}"))?)
    };
    autodetect(source).map_err(|e| format!("cannot read {path:?}: {e}"))
}

/// Profiling is explicitly opt-in; the ordinary path uses no timing wrappers.
pub fn open_profiled(
    path: &str,
    metrics: profile::SharedInput,
) -> Result<Box<dyn BufRead + Send>, String> {
    let source: Box<dyn Read + Send> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).map_err(|e| format!("cannot open {path:?}: {e}"))?)
    };
    let source = profile::timed_source(source, metrics.clone());
    autodetect_profiled(source, Some(metrics)).map_err(|e| format!("cannot read {path:?}: {e}"))
}

fn autodetect(source: Box<dyn Read + Send>) -> Result<Box<dyn BufRead + Send>, String> {
    autodetect_profiled(source, None)
}

fn autodetect_profiled(
    mut source: Box<dyn Read + Send>,
    metrics: Option<profile::SharedInput>,
) -> Result<Box<dyn BufRead + Send>, String> {
    // A pipe may supply only one byte per read. Do not assume fill_buf() can
    // provide both magic bytes in a single call.
    let mut magic = [0_u8; 2];
    let mut used = 0;
    while used < magic.len() {
        match source.read(&mut magic[used..]) {
            Ok(0) => break,
            Ok(n) => used += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.to_string()),
        }
    }
    let is_gzip = used == 2 && magic == [0x1f, 0x8b];
    let source = Cursor::new(magic[..used].to_vec()).chain(source);
    if let Some(metrics) = metrics {
        metrics.lock().map_err(|_| "profile lock poisoned")?.gzip = is_gzip;
        let decoded: Box<dyn Read + Send> = if is_gzip {
            Box::new(MultiGzDecoder::new(BufReader::with_capacity(
                BUFFER_SIZE,
                source,
            )))
        } else {
            Box::new(source)
        };
        return Ok(Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            profile::timed_decoded(decoded, metrics),
        )));
    }
    if is_gzip {
        let compressed = BufReader::with_capacity(BUFFER_SIZE, source);
        Ok(Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            MultiGzDecoder::new(compressed),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, source)))
    }
}

#[derive(Clone, Copy)]
enum Format {
    Fasta,
    Fastq,
}

pub struct Reader {
    source: Box<dyn BufRead + Send>,
    format: Option<Format>,
    line: Vec<u8>,
    line_number: u64,
    pending_header: Option<Vec<u8>>,
    finished: bool,
    pair_header: Option<String>,
    quality_max_bases: Option<usize>,
}

impl Reader {
    pub fn new(source: Box<dyn BufRead + Send>) -> Self {
        Self {
            source,
            format: None,
            line: Vec::new(),
            line_number: 0,
            pending_header: None,
            finished: false,
            pair_header: None,
            quality_max_bases: None,
        }
    }

    /// Keep only the most recently parsed header for streaming mate validation.
    pub fn paired(source: Box<dyn BufRead + Send>) -> Self {
        let mut reader = Self::new(source);
        reader.pair_header = Some(String::new());
        reader
    }

    /// Retain complete FASTQ quality only when the entire sequence fits this
    /// finite limit. Larger records are still parsed and validated completely.
    pub fn with_quality(source: Box<dyn BufRead + Send>, max_bases: usize) -> Self {
        let mut reader = Self::new(source);
        reader.quality_max_bases = Some(max_bases);
        reader
    }

    /// Combine bounded quality retention with optional streaming mate headers.
    pub fn retain_quality(&mut self, max_bases: usize) {
        self.quality_max_bases = Some(max_bases);
    }

    pub fn pair_header(&self) -> &str {
        self.pair_header.as_deref().unwrap_or("")
    }

    pub fn is_fastq(&self) -> bool {
        matches!(self.format, Some(Format::Fastq))
    }

    pub fn next_record(&mut self) -> Result<Option<Record>, String> {
        if self.finished {
            return Ok(None);
        }
        if let Some(header) = self.pending_header.take() {
            self.line = header;
        } else {
            loop {
                if !self.read_line()? {
                    self.finished = true;
                    return Ok(None);
                }
                // Tolerate leading blank lines, but reject blank lines between
                // FASTQ records because those often indicate broken wrapping.
                if self.format.is_none() && self.line.is_empty() {
                    continue;
                }
                break;
            }
        }

        let format = match self.format {
            Some(format) => format,
            None => {
                let format = match self.line.first() {
                    Some(b'>') => Format::Fasta,
                    Some(b'@') => Format::Fastq,
                    _ => return Err(self.error("expected a FASTA '>' or FASTQ '@' header")),
                };
                self.format = Some(format);
                format
            }
        };
        let marker = match format {
            Format::Fasta => b'>',
            Format::Fastq => b'@',
        };
        if self.line.first() != Some(&marker) {
            return Err(self.error(&format!("expected '{}' record header", marker as char)));
        }
        let id = self.header_id(&self.line[1..])?;
        if let Some(header) = &mut self.pair_header {
            header.clear();
            header.push_str(
                std::str::from_utf8(&self.line[1..])
                    .map_err(|_| "paired FASTQ header must be UTF-8")?,
            );
        }
        let (seq, quality) = match format {
            Format::Fasta => (self.read_fasta(&id)?, None),
            Format::Fastq => self.read_fastq(&id)?,
        };
        Ok(Some(Record { id, seq, quality }))
    }

    fn read_line(&mut self) -> Result<bool, String> {
        self.line.clear();
        let count = self
            .source
            .read_until(b'\n', &mut self.line)
            .map_err(|e| format!("input error near line {}: {e}", self.line_number + 1))?;
        if count == 0 {
            return Ok(false);
        }
        self.line_number += 1;
        if self.line.last() == Some(&b'\n') {
            self.line.pop();
            if self.line.last() == Some(&b'\r') {
                self.line.pop();
            }
        }
        Ok(true)
    }

    fn read_fasta(&mut self, id: &str) -> Result<Vec<u8>, String> {
        let mut seq = Vec::new();
        while self.read_line()? {
            if self.line.first() == Some(&b'>') {
                self.pending_header = Some(std::mem::take(&mut self.line));
                break;
            }
            self.append_sequence(&mut seq, id, true)?;
        }
        if seq.is_empty() {
            return Err(self.error(&format!("record {id:?} has an empty sequence")));
        }
        Ok(seq)
    }

    fn read_fastq(&mut self, id: &str) -> Result<(Vec<u8>, Option<Vec<u8>>), String> {
        let mut seq = Vec::new();
        loop {
            if !self.read_line()? {
                return Err(self.error(&format!(
                    "truncated FASTQ record {id:?}: missing '+' separator"
                )));
            }
            if self.line.first() == Some(&b'+') {
                let repeated = &self.line[1..];
                if !repeated.is_empty() {
                    let repeated_id = self.header_id(repeated)?;
                    if repeated_id != id {
                        return Err(self.error(&format!(
                            "FASTQ '+' identifier {repeated_id:?} does not match {id:?}"
                        )));
                    }
                }
                break;
            }
            if self.line.is_empty() {
                return Err(self.error(&format!("empty FASTQ sequence line in record {id:?}")));
            }
            self.append_sequence(&mut seq, id, false)?;
        }
        if seq.is_empty() {
            return Err(self.error(&format!("record {id:?} has an empty sequence")));
        }
        let mut quality_len = 0_usize;
        let mut quality = self
            .quality_max_bases
            .filter(|&limit| seq.len() <= limit)
            .map(|_| Vec::with_capacity(seq.len()));
        while quality_len < seq.len() {
            if !self.read_line()? {
                return Err(self.error(&format!(
                    "truncated FASTQ quality for {id:?}: expected {} characters, found {quality_len}",
                    seq.len()
                )));
            }
            if self.line.is_empty() {
                return Err(self.error(&format!("empty FASTQ quality line in record {id:?}")));
            }
            // A reduction over the whole quality line lets the compiler use
            // vector instructions; all characters still receive validation.
            if self.line.iter().fold(0_u8, |invalid, &b| {
                invalid | u8::from(b.wrapping_sub(33) > 93)
            }) != 0
            {
                return Err(
                    self.error(&format!("invalid FASTQ quality character in record {id:?}"))
                );
            }
            quality_len = quality_len.checked_add(self.line.len()).ok_or_else(|| {
                self.error(&format!("FASTQ quality length overflow in record {id:?}"))
            })?;
            if quality_len > seq.len() {
                return Err(self.error(&format!(
                    "FASTQ quality longer than sequence for {id:?}: {} bases, {quality_len} quality characters",
                    seq.len()
                )));
            }
            if let Some(retained) = &mut quality {
                retained.extend_from_slice(&self.line);
            }
        }
        Ok((seq, quality))
    }

    fn append_sequence(
        &self,
        seq: &mut Vec<u8>,
        id: &str,
        allow_space: bool,
    ) -> Result<(), String> {
        // Most WGS reads contain only uppercase A/C/G/T/N. Validate the
        // whole slice, then copy it in bulk instead of pushing each byte.
        // Keep the original path for lowercase, IUPAC, FASTA spacing and
        // precise invalid-byte diagnostics. The reduction is vectorizable.
        let invalid = self.line.iter().fold(0_u8, |invalid, &base| {
            invalid | u8::from(!matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
        });
        if invalid == 0 {
            seq.extend_from_slice(&self.line);
            return Ok(());
        }
        seq.reserve(self.line.len());
        for (offset, &raw) in self.line.iter().enumerate() {
            if allow_space && matches!(raw, b' ' | b'\t') {
                continue;
            }
            let base = raw.to_ascii_uppercase();
            if !matches!(
                base,
                b'A' | b'C'
                    | b'G'
                    | b'T'
                    | b'R'
                    | b'Y'
                    | b'S'
                    | b'W'
                    | b'K'
                    | b'M'
                    | b'B'
                    | b'D'
                    | b'H'
                    | b'V'
                    | b'N'
            ) {
                return Err(self.error(&format!(
                    "invalid DNA byte 0x{raw:02x} in record {id:?}, column {}",
                    offset + 1
                )));
            }
            seq.push(base);
        }
        Ok(())
    }

    fn header_id(&self, header: &[u8]) -> Result<String, String> {
        if header.iter().any(|&b| b.is_ascii_control() && b != b'\t') {
            return Err(self.error("control character in record header"));
        }
        let text = std::str::from_utf8(header)
            .map_err(|_| self.error("record header is not valid UTF-8"))?;
        text.split_whitespace()
            .next()
            .map(str::to_owned)
            .ok_or_else(|| self.error("record header has no identifier"))
    }

    fn error(&self, message: &str) -> String {
        format!("line {}: {message}", self.line_number)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::io::Write;

    fn parse(data: &[u8]) -> Result<Vec<Record>, String> {
        let mut reader = Reader::new(autodetect(Box::new(Cursor::new(data.to_vec())))?);
        let mut records = Vec::new();
        while let Some(record) = reader.next_record()? {
            records.push(record);
        }
        assert_eq!(reader.next_record()?, None);
        Ok(records)
    }

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn multiline_fasta_crlf_and_iupac() {
        let records =
            parse(b"\r\n>chr1 description\r\nac gt\r\nryswkmbdhvn\r\n\r\n>chr2\r\nT\tTAGGG")
                .unwrap();
        assert_eq!(
            records,
            vec![
                Record {
                    id: "chr1".into(),
                    seq: b"ACGTRYSWKMBDHVN".to_vec(),
                    quality: None,
                },
                Record {
                    id: "chr2".into(),
                    seq: b"TTAGGG".to_vec(),
                    quality: None,
                },
            ]
        );
    }

    #[test]
    fn wrapped_fastq_allows_header_markers_in_quality() {
        let records = parse(b"@read1 description\r\nacg\r\ntnn\r\n+read1 description\r\n@@@\r\n+++\r\n@read2\r\nT\r\n+\r\n!").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].seq, b"ACGTNN");
        assert_eq!(records[1].id, "read2");
    }

    #[test]
    fn rejects_empty_or_invalid_sequences_and_headers() {
        for input in [
            b">empty\n>next\nAC\n".as_slice(),
            b"@empty\n+\n\n",
            b"> \nAC\n",
            b"@\nAC\n+\n!!\n",
            b">a\nACU\n",
            b">a\nAC-\n",
            b">a\nAC\xff\n",
            b"@a\nAC GT\n+\n!!!!!\n",
            b">a\nAC\n@b\nTT\n",
            b"ACGT\n",
        ] {
            assert!(
                parse(input).is_err(),
                "accepted {:?}",
                String::from_utf8_lossy(input)
            );
        }
        assert!(parse(b"").unwrap().is_empty());
    }

    #[test]
    fn rejects_truncated_and_mismatched_fastq() {
        for input in [
            b"@a\nACGT\n".as_slice(),
            b"@a\nACGT\n+\n!!!\n",
            b"@a\nACGT\n+\n!!!!!\n",
            b"@a\nACGT\n+b\n!!!!\n",
            b"@a\nACGT\n+\n!! !!\n",
            b"@a\nACGT\n+\n!!!!\nnot-a-header\n",
            b"@a\nACGT\n+\n!!!!\n@b\n",
            b"@a\nACGT\n+\n!!!!\n\n",
        ] {
            assert!(
                parse(input).is_err(),
                "accepted {:?}",
                String::from_utf8_lossy(input)
            );
        }
    }

    #[test]
    fn gzip_concatenated_members_and_integrity() {
        let mut compressed = gzip(b">a\nAC\n");
        compressed.extend(gzip(b">b\nGT\n"));
        let records = parse(&compressed).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[1].seq, b"GT");

        assert!(parse(&compressed[..compressed.len() - 3]).is_err());
        let mut corrupt = gzip(b">a\nAC\n");
        let checksum = corrupt.len() - 8;
        corrupt[checksum] ^= 1;
        assert!(parse(&corrupt).is_err());

        compressed.extend_from_slice(b"trailing garbage");
        assert!(parse(&compressed).is_err());
    }

    #[test]
    fn sequence_fast_path_preserves_all_byte_classification_and_columns() {
        // Cover each byte on both sides of common vector-width boundaries.
        // The scalar expectation also checks the partial sequence on error.
        for allow_space in [false, true] {
            for offset in [0, 1, 15, 16, 31, 32, 63, 64, 127, 149] {
                for raw in 0..=u8::MAX {
                    let mut reader = Reader::new(Box::new(Cursor::new(Vec::<u8>::new())));
                    reader.line = vec![b'A'; 150];
                    reader.line[offset] = raw;
                    reader.line_number = 7;
                    let mut observed = b"TT".to_vec();
                    let result = reader.append_sequence(&mut observed, "read", allow_space);
                    let base = raw.to_ascii_uppercase();
                    let skipped = allow_space && matches!(raw, b' ' | b'\t');
                    let valid = b"ACGTRYSWKMBDHVN".contains(&base);
                    if skipped || valid {
                        let mut expected = b"TT".to_vec();
                        expected.extend_from_slice(&vec![b'A'; offset]);
                        if !skipped {
                            expected.push(base);
                        }
                        expected.extend_from_slice(&vec![b'A'; 149 - offset]);
                        assert!(result.is_ok(), "byte {raw}, offset {offset}");
                        assert_eq!(observed, expected);
                    } else {
                        assert_eq!(
                            result.unwrap_err(),
                            format!(
                                "line 7: invalid DNA byte 0x{raw:02x} in record \"read\", column {}",
                                offset + 1
                            )
                        );
                        assert_eq!(observed.len(), 2 + offset);
                    }
                }
            }
        }
    }

    #[test]
    fn quality_validation_covers_all_bytes_and_vector_boundaries() {
        for offset in [0, 1, 15, 16, 31, 32, 63, 64, 127, 149] {
            for raw in 0..=u8::MAX {
                let mut input = b"@read\n".to_vec();
                input.extend_from_slice(&[b'A'; 150]);
                input.extend_from_slice(b"\n+\n");
                let mut quality = [b'I'; 150];
                quality[offset] = raw;
                input.extend_from_slice(&quality);
                input.push(b'\n');
                assert_eq!(
                    parse(&input).is_ok(),
                    (33..=126).contains(&raw),
                    "byte {raw}, offset {offset}"
                );
            }
        }
    }

    #[test]
    fn fast_path_and_iupac_fallback_can_alternate_across_wrapped_lines() {
        let mut state = 52912_u64;
        let alphabet = b"acgtnryswkmbdhvACGTN";
        let mut expected = Vec::new();
        let mut input = b"@read 1:N:0:42\r\n".to_vec();
        for line in 0..25 {
            for _ in 0..(1 + line * 7) {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let raw = if line % 2 == 0 {
                    b"ACGTN"[(state >> 32) as usize % 5]
                } else {
                    alphabet[(state >> 32) as usize % alphabet.len()]
                };
                input.push(raw);
                expected.push(raw.to_ascii_uppercase());
            }
            input.extend_from_slice(b"\r\n");
        }
        input.extend_from_slice(b"+read description\r\n");
        for chunk in expected.chunks(47) {
            input.extend_from_slice(&vec![b'@'; chunk.len()]);
            input.extend_from_slice(b"\r\n");
        }
        let mut reader = Reader::paired(Box::new(Cursor::new(input)));
        assert_eq!(reader.next_record().unwrap().unwrap().seq, expected);
        assert_eq!(reader.pair_header(), "read 1:N:0:42");
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn magic_detection_handles_short_reads() {
        struct OneByte(Cursor<Vec<u8>>);
        impl Read for OneByte {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                let len = buf.len().min(1);
                self.0.read(&mut buf[..len])
            }
        }
        let source = OneByte(Cursor::new(gzip(b"@a\nAC\n+\n!!\n")));
        let mut reader = Reader::new(autodetect(Box::new(source)).unwrap());
        assert_eq!(reader.next_record().unwrap().unwrap().seq, b"AC");
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn quality_retention_is_opt_in_complete_and_wrapping_independent() {
        let data = b"@same comment\r\nacg\r\ntnr\r\n+same\r\n!@+\r\nI~?\r\n";
        let mut off = Reader::new(Box::new(Cursor::new(data.to_vec())));
        assert_eq!(off.next_record().unwrap().unwrap().quality, None);
        let mut on = Reader::with_quality(Box::new(Cursor::new(data.to_vec())), 6);
        let record = on.next_record().unwrap().unwrap();
        assert_eq!(record.id, "same");
        assert_eq!(record.seq, b"ACGTNR");
        assert_eq!(record.quality.unwrap(), b"!@+I~?");
        assert!(on.next_record().unwrap().is_none());
    }

    #[test]
    fn oversized_quality_is_discarded_but_still_validated() {
        let mut valid =
            Reader::with_quality(Box::new(Cursor::new(b"@x\nACGT\n+\nIIII\n".to_vec())), 3);
        let record = valid.next_record().unwrap().unwrap();
        assert_eq!(record.seq, b"ACGT");
        assert_eq!(record.quality, None);
        for quality in ["III", "II I", "IIIII"] {
            let mut bad = Reader::with_quality(
                Box::new(Cursor::new(
                    format!("@x\nACGT\n+\n{quality}\n").into_bytes(),
                )),
                3,
            );
            assert!(bad.next_record().is_err());
        }
    }

    #[test]
    fn fasta_never_invents_quality_when_capture_is_enabled() {
        let mut reader = Reader::with_quality(Box::new(Cursor::new(b">x\nac gt\n".to_vec())), 100);
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.seq, b"ACGT");
        assert_eq!(record.quality, None);
    }
}
