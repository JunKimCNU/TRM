//! Optional forensic export of every native targeted-family tract and its flanks.
//! This example does not change detection, candidate classification, or adequacy.
#[allow(dead_code)]
#[path = "../src/detect.rs"]
mod detect;
#[allow(dead_code)]
#[path = "../src/input.rs"]
mod input;
#[allow(dead_code)]
#[path = "../src/motif.rs"]
mod motif;

use detect::{Config, TargetIndex, Tract};
use serde_json::json;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Instant, UNIX_EPOCH};

const HELP: &str = "Usage: export_target_tracts INPUT --motif OUTWARD_MOTIF --output-dir NEW_DIR\n\
    [--min-anchor-bp 200] [--max-anchor-bp 1000]\n\
Single FASTA/FASTQ input, gzip detected from bytes. New output directory required.\n\
All native targeted tracts are exported (no 500 bp length or min-read filter).\n\
Coordinates are 0-based half-open; read and tract ordinals are 1-based.\n\
Strand is relative to the supplied motif; anchors keep the original read strand.\n\
Anchor lengths must satisfy 200 <= min <= max <= 1000.\n\
Requires sha256sum for bounded-memory hashes of completed output files.\n\
summary.json is written only after validated EOF and flushed output files.\n";

struct Args {
    input: PathBuf,
    motif: String,
    output: PathBuf,
    min_anchor: usize,
    max_anchor: usize,
}

fn arguments() -> Result<Option<Args>, String> {
    let mut input = None;
    let mut motif = None;
    let mut output = None;
    let (mut min_anchor, mut max_anchor) = (200, 1000);
    let mut min_seen = false;
    let mut max_seen = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" {
            print!("{HELP}");
            return Ok(None);
        }
        match arg.as_str() {
            "--motif" if motif.is_none() => {
                motif = Some(
                    args.next()
                        .ok_or("missing --motif value")?
                        .to_ascii_uppercase(),
                )
            }
            "--output-dir" if output.is_none() => {
                output = Some(PathBuf::from(
                    args.next().ok_or("missing --output-dir value")?,
                ))
            }
            "--min-anchor-bp" if !min_seen => {
                min_anchor = args
                    .next()
                    .ok_or("missing --min-anchor-bp value")?
                    .parse()
                    .map_err(|_| "invalid --min-anchor-bp value")?;
                min_seen = true;
            }
            "--max-anchor-bp" if !max_seen => {
                max_anchor = args
                    .next()
                    .ok_or("missing --max-anchor-bp value")?
                    .parse()
                    .map_err(|_| "invalid --max-anchor-bp value")?;
                max_seen = true;
            }
            _ if !arg.starts_with('-') && input.is_none() => input = Some(PathBuf::from(arg)),
            _ => return Err(format!("unknown, repeated, or unexpected argument {arg:?}")),
        }
    }
    let input = input.ok_or("one input FASTA/FASTQ path is required")?;
    let motif = motif.ok_or("--motif is required")?;
    let output = output.ok_or("--output-dir is required")?;
    if !(3..=32).contains(&motif.len()) || !motif::is_dna(motif.as_bytes()) {
        return Err("motif must contain 3..32 A/C/G/T bases".into());
    }
    if motif::primitive_len(motif.as_bytes()) != motif.len() {
        return Err("motif must be a primitive repeat unit, not copies of a shorter unit".into());
    }
    if min_anchor < 200 || min_anchor > max_anchor || max_anchor > 1000 {
        return Err("anchor lengths require 200 <= min <= max <= 1000".into());
    }
    Ok(Some(Args {
        input,
        motif,
        output,
        min_anchor,
        max_anchor,
    }))
}

fn create(path: &Path) -> Result<BufWriter<File>, String> {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map(BufWriter::new)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))
}

fn write_error(e: std::io::Error) -> String {
    format!("output error: {e}")
}

fn output_fingerprint(path: &Path) -> Result<serde_json::Value, String> {
    let source = File::open(path).map_err(write_error)?;
    let size = source.metadata().map_err(write_error)?.len();
    // Hash standard input: no shell, filename parsing, or in-memory file copy.
    let result = Command::new("sha256sum")
        .stdin(Stdio::from(source))
        .output()
        .map_err(|e| format!("sha256sum is required for output verification: {e}"))?;
    if !result.status.success() {
        return Err(format!(
            "sha256sum failed for {}: {}",
            path.display(),
            String::from_utf8_lossy(&result.stderr)
        ));
    }
    let stdout = String::from_utf8(result.stdout).map_err(|e| format!("sha256sum output: {e}"))?;
    let hash = stdout
        .split_whitespace()
        .next()
        .ok_or("empty sha256sum output")?;
    if hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err("invalid sha256sum digest".into());
    }
    Ok(json!({"sha256": hash.to_ascii_lowercase(), "size_bytes": size}))
}

fn safe_read_id(id: &str) -> String {
    let mut bytes = Vec::new();
    for b in id.bytes() {
        if b == b'%' || b.is_ascii_control() {
            bytes.extend_from_slice(format!("%{b:02X}").as_bytes());
        } else {
            bytes.push(b);
        }
    }
    String::from_utf8(bytes).expect("escaped valid UTF-8")
}

fn intervals(tracts: &[Tract]) -> Vec<(usize, usize)> {
    let mut union = Vec::<(usize, usize)>::new();
    for t in tracts {
        match union.last_mut() {
            Some(last) if t.start <= last.1 => last.1 = last.1.max(t.end),
            _ => union.push((t.start, t.end)),
        }
    }
    union
}

struct Anchor {
    start: usize,
    end: usize,
    stop: &'static str,
}

fn anchor(
    seq: &[u8],
    t: &Tract,
    union: &[(usize, usize)],
    index: usize,
    left: bool,
    maximum: usize,
) -> Anchor {
    let (merged_start, merged_end) = union[index];
    if left {
        if merged_start < t.start {
            return Anchor {
                start: t.start,
                end: t.start,
                stop: "same_family_tract",
            };
        }
        let hard = if index == 0 { 0 } else { union[index - 1].1 };
        let mut start = t.start.saturating_sub(maximum).max(hard);
        let mut stop = if start == hard {
            if index == 0 {
                "read_boundary"
            } else {
                "same_family_tract"
            }
        } else {
            "maximum_length"
        };
        if let Some(p) = (start..t.start)
            .rev()
            .find(|&p| !motif::is_dna(&seq[p..p + 1]))
        {
            start = p + 1;
            stop = "ambiguous_base";
        }
        Anchor {
            start,
            end: t.start,
            stop,
        }
    } else {
        if merged_end > t.end {
            return Anchor {
                start: t.end,
                end: t.end,
                stop: "same_family_tract",
            };
        }
        let hard = union.get(index + 1).map_or(seq.len(), |x| x.0);
        let mut end = t.end.saturating_add(maximum).min(hard);
        let mut stop = if end == hard {
            if index + 1 == union.len() {
                "read_boundary"
            } else {
                "same_family_tract"
            }
        } else {
            "maximum_length"
        };
        if let Some(p) = (t.end..end).find(|&p| !motif::is_dna(&seq[p..p + 1])) {
            end = p;
            stop = "ambiguous_base";
        }
        Anchor {
            start: t.end,
            end,
            stop,
        }
    }
}

fn export(args: Args) -> Result<(), String> {
    let started = Instant::now();
    let input_path = fs::canonicalize(&args.input).map_err(|e| format!("input path: {e}"))?;
    let before = fs::metadata(&input_path).map_err(|e| format!("input metadata: {e}"))?;
    if !before.is_file() {
        return Err("input must be a regular file".into());
    }
    let input_string = input_path.to_str().ok_or("input path must be UTF-8")?;
    let mut reader = input::Reader::new(input::open(input_string)?);
    if let Some(parent) = args.output.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| format!("output parent: {e}"))?;
    }
    fs::create_dir(&args.output).map_err(|e| format!("new output directory required: {e}"))?;
    let mut tsv = create(&args.output.join("tracts.tsv.partial"))?;
    let mut fasta = create(&args.output.join("anchors.fasta.partial"))?;
    writeln!(tsv, "tract_id\tread_ordinal\tread_id\tread_length\ttract_index\tmotif\tfamily_key\tstart0\tend0\ttract_length\tidentity\tstrand\tleft_end_distance\tright_end_distance\tleft_anchor_id\tleft_anchor_start0\tleft_anchor_end0\tleft_anchor_length\tleft_anchor_status\tleft_anchor_stop_reason\tright_anchor_id\tright_anchor_start0\tright_anchor_end0\tright_anchor_length\tright_anchor_status\tright_anchor_stop_reason").map_err(write_error)?;
    let cfg = Config::default();
    let (family, supplied_strand) = motif::canonical(args.motif.as_bytes());
    let targets = TargetIndex::new(vec![family.clone()], &cfg);
    let (mut records, mut bases, mut target_reads, mut tract_count, mut anchor_count) =
        (0u64, 0u64, 0u64, 0u64, 0u64);
    let (mut max_record, mut max_tract) = (0usize, 0usize);
    let (mut tracts_with_anchor, mut tracts_without_anchor) = (0u64, 0u64);
    while let Some(record) = reader.next_record()? {
        records += 1;
        bases = bases
            .checked_add(record.seq.len() as u64)
            .ok_or("base counter overflow")?;
        max_record = max_record.max(record.seq.len());
        let mut hits = targets.discover(&record.seq, &cfg);
        let Some(mut tracts) = hits.remove(&family) else {
            continue;
        };
        if tracts.is_empty() {
            continue;
        }
        target_reads += 1;
        tracts.sort_by_key(|t| (t.start, std::cmp::Reverse(t.end), t.strand));
        let union = intervals(&tracts);
        let mut union_index = 0;
        for (i, t) in tracts.iter().enumerate() {
            while union[union_index].1 <= t.start {
                union_index += 1;
            }
            let length = t.end - t.start;
            max_tract = max_tract.max(length);
            tract_count += 1;
            let tract_id = format!("r{records:012}_t{:06}", i + 1);
            let strand = match (t.strand, supplied_strand) {
                (Some(a), Some(b)) if a == b => "+",
                (Some(_), Some(_)) => "-",
                _ => "ambiguous",
            };
            write!(tsv, "{tract_id}\t{records}\t{}\t{}\t{}\t{}\t{family}\t{}\t{}\t{length}\t{:.12}\t{strand}\t{}\t{}",
                   safe_read_id(&record.id), record.seq.len(), i + 1, args.motif, t.start, t.end,
                   t.identity, t.start, record.seq.len() - t.end).map_err(write_error)?;
            let mut exported = false;
            for (left, suffix) in [(true, "L"), (false, "R")] {
                let a = anchor(&record.seq, t, &union, union_index, left, args.max_anchor);
                let length = a.end - a.start;
                let usable = length >= args.min_anchor;
                let query_id = if usable {
                    format!("{tract_id}_{suffix}")
                } else {
                    String::new()
                };
                let status = if usable { "exported" } else { "below_minimum" };
                write!(
                    tsv,
                    "\t{query_id}\t{}\t{}\t{length}\t{status}\t{}",
                    a.start, a.end, a.stop
                )
                .map_err(write_error)?;
                if usable {
                    writeln!(fasta, ">{query_id}").map_err(write_error)?;
                    fasta
                        .write_all(&record.seq[a.start..a.end])
                        .map_err(write_error)?;
                    fasta.write_all(b"\n").map_err(write_error)?;
                    anchor_count += 1;
                    exported = true;
                }
            }
            writeln!(tsv).map_err(write_error)?;
            if exported {
                tracts_with_anchor += 1;
            } else {
                tracts_without_anchor += 1;
            }
        }
    }
    if records == 0 {
        return Err("input contains no sequence records".into());
    }
    let after = fs::metadata(&input_path).map_err(|e| format!("input metadata after EOF: {e}"))?;
    if before.len() != after.len() || before.modified().ok() != after.modified().ok() {
        return Err("input changed while exporting".into());
    }
    tsv.flush().map_err(write_error)?;
    fasta.flush().map_err(write_error)?;
    drop(tsv);
    drop(fasta);
    for name in ["tracts.tsv", "anchors.fasta"] {
        fs::rename(
            args.output.join(format!("{name}.partial")),
            args.output.join(name),
        )
        .map_err(write_error)?;
    }
    let modified = before
        .modified()
        .ok()
        .and_then(|x| x.duration_since(UNIX_EPOCH).ok());
    let tsv_fingerprint = output_fingerprint(&args.output.join("tracts.tsv"))?;
    let fasta_fingerprint = output_fingerprint(&args.output.join("anchors.fasta"))?;
    let summary = json!({
        "schema_version": 1,
        "program": "export_target_tracts", "detector_version": env!("CARGO_PKG_VERSION"),
        "status": "complete", "complete_input": true, "validated_eof": true, "sampled_prefix": false,
        "input": {"path": input_path, "size_bytes": before.len(),
                  "modified_unix_seconds": modified.map(|x| x.as_secs()),
                  "modified_subsecond_nanoseconds": modified.map(|x| x.subsec_nanos()),
                  "format": if reader.is_fastq() { "fastq" } else { "fasta" }},
        "motif": args.motif, "family_key": family, "config": cfg,
        "length_filter_applied": false, "minimum_read_support_filter_applied": false,
        "targeted_strands": "both", "strand_reference": "supplied motif; reverse-complement symmetric families are ambiguous",
        "records": records, "bases": bases, "target_reads": target_reads, "tracts": tract_count,
        "max_record_bases": max_record, "max_tract_bases": max_tract,
        "anchors": anchor_count, "tracts_with_anchor": tracts_with_anchor, "tracts_without_anchor": tracts_without_anchor,
        "outputs": {"tracts.tsv": tsv_fingerprint, "anchors.fasta": fasta_fingerprint},
        "min_anchor_bases": args.min_anchor, "max_anchor_bases": args.max_anchor,
        "coordinates": "0-based half-open; read_ordinal and tract_index are 1-based",
        "read_id_encoding": "percent-escaped percent sign and ASCII control bytes; IDs do not determine tract/query uniqueness",
        "anchor_orientation": "original input strand; sequences are not reverse complemented",
        "anchor_boundaries": "stop at read edge, configured maximum length, non-ACGT base, or another detected same-family tract; short unseeded repeat fragments and other-family repeats are not masked",
        "memory_scope": "one input record and its native targeted-family tracts; no collection across records",
        "elapsed_seconds": started.elapsed().as_secs_f64(),
        "limitations": ["Native seed/extension and substitution-only alignment limitations still apply.",
                        "An exported flank is an adjacent sequence, not proof of unique subtelomeric identity or physical chromosome end.",
                        "Only detected tracts of the selected family bound anchors; other repeats and missed interrupted tracts can remain in an anchor."]
    });
    let mut out = create(&args.output.join("summary.json.partial"))?;
    serde_json::to_writer_pretty(&mut out, &summary).map_err(|e| format!("summary output: {e}"))?;
    writeln!(out).map_err(write_error)?;
    out.flush().map_err(write_error)?;
    drop(out);
    fs::rename(
        args.output.join("summary.json.partial"),
        args.output.join("summary.json"),
    )
    .map_err(write_error)?;
    Ok(())
}

fn main() {
    let result = arguments().and_then(|args| args.map_or(Ok(()), export));
    if let Err(error) = result {
        eprintln!("export_target_tracts: {error}");
        std::process::exit(1);
    }
}
