//! Optional, explicitly supplied irregular-array model. This exporter never
//! changes native TRM calls and never labels a modeled span an exact repeat.
#[allow(dead_code)]
#[path = "../src/input.rs"]
mod input;
#[allow(dead_code)]
#[path = "../src/irregular.rs"]
mod irregular;
#[allow(dead_code)]
#[path = "../src/memory.rs"]
mod memory;
#[allow(dead_code)]
#[path = "../src/motif.rs"]
mod motif;

use irregular::{Array, Config, Definition, ScanResult, Scanner};
use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs::{self, File, Metadata, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const HELP: &str = "Usage: scan_irregular_arrays INPUT [INPUT...] --definitions JSON --output-dir NEW_DIR\n\
    [--threads 2] [--memory-limit-mib N] [--batch-bases 2000000]\n\
    [--min-span 100] [--min-exact-anchors 4] [--max-record-bases 25000000]\n\
    [--max-nodes 2000000] [--export-anchors] [--max-records 0]\n\
Definitions: JSON array of {\"id\":\"group\",\"units\":[\"OUTWARD_UNIT\",...]}.\n\
Regular FASTA/FASTQ files only; gzip detected from bytes. New output directory required.\n\
Every accepted model array is emitted; computational limits abort without a summary.\n\
Modeled span, exact unit bases, edits and parse/orientation ambiguity remain distinct.\n\
Optional adjacent anchors are 200..1000 bp, in the original read strand.\n\
They stop at ambiguity, a read edge, or any reported array from any group/strand.\n\
Coordinates are zero-based half-open; source/read/array ordinals are one-based.\n\
--max-records > 0 explicitly permits an incomplete prefix; it is not a full-input result.\n\
--memory-limit-mib caps virtual address space, not RSS; a full record can exceed a batch.\n\
Requires sha256sum. Input sizes, modification times and hashes are checked before/after.\n\
summary.json is written last after validated requested input and output fingerprints.\n";

struct Args {
    inputs: Vec<PathBuf>, definitions: PathBuf, output: PathBuf,
    threads: usize, memory_limit_mib: Option<u64>, batch_bases: usize,
    config: Config, export_anchors: bool, max_records: u64,
}

fn parse_arguments(argv: impl IntoIterator<Item = String>) -> Result<Option<Args>, String> {
    let mut args = argv.into_iter();
    let (mut inputs, mut definitions, mut output) = (Vec::new(), None, None);
    let (mut threads, mut memory_limit_mib, mut batch_bases) = (2, None, 2_000_000);
    let (mut export_anchors, mut max_records) = (false, 0);
    let mut config = Config::default();
    let mut seen = HashSet::new();
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" { print!("{HELP}"); return Ok(None); }
        if arg == "--version" {
            println!("scan_irregular_arrays model 0.1.0 (trm {})", env!("CARGO_PKG_VERSION"));
            return Ok(None);
        }
        if !arg.starts_with('-') { inputs.push(PathBuf::from(arg)); continue; }
        if !seen.insert(arg.clone()) { return Err(format!("repeated argument {arg}")); }
        if arg == "--export-anchors" { export_anchors = true; continue; }
        let value = args.next().ok_or_else(|| format!("missing {arg} value"))?;
        let number = |name: &str| value.parse::<usize>().map_err(|_| format!("invalid {name} value"));
        match arg.as_str() {
            "--definitions" => definitions = Some(PathBuf::from(value)),
            "--output-dir" => output = Some(PathBuf::from(value)),
            "--threads" => threads = number(&arg)?,
            "--memory-limit-mib" => memory_limit_mib = Some(value.parse::<u64>().map_err(|_| "invalid --memory-limit-mib value")?),
            "--batch-bases" => batch_bases = number(&arg)?,
            "--min-span" => config.min_span = number(&arg)?,
            "--min-exact-anchors" => config.min_exact_anchors = number(&arg)?,
            "--max-record-bases" => config.max_record_bases = number(&arg)?,
            "--max-nodes" => config.max_nodes = number(&arg)?,
            "--max-records" => max_records = value.parse::<u64>().map_err(|_| "invalid --max-records value")?,
            _ => return Err(format!("unknown argument {arg}")),
        }
    }
    if inputs.is_empty() { return Err("at least one input FASTA/FASTQ path is required".into()); }
    if threads == 0 || batch_bases == 0 { return Err("threads and batch bases must be positive".into()); }
    if memory_limit_mib == Some(0) { return Err("memory limit must be positive".into()); }
    Ok(Some(Args { inputs, definitions: definitions.ok_or("--definitions is required")?,
        output: output.ok_or("--output-dir is required")?, threads, memory_limit_mib,
        batch_bases, config, export_anchors, max_records }))
}

fn io_error(error: std::io::Error) -> String { error.to_string() }

fn create(path: &Path) -> Result<BufWriter<File>, String> {
    OpenOptions::new().write(true).create_new(true).open(path).map(BufWriter::new)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct Fingerprint { sha256: String, size_bytes: u64 }

fn fingerprint(path: &Path) -> Result<Fingerprint, String> {
    let file = File::open(path).map_err(io_error)?;
    let size_bytes = file.metadata().map_err(io_error)?.len();
    let result = Command::new("sha256sum").stdin(Stdio::from(file)).output()
        .map_err(|e| format!("sha256sum is required: {e}"))?;
    if !result.status.success() {
        return Err(format!("sha256sum failed for {}: {}", path.display(), String::from_utf8_lossy(&result.stderr)));
    }
    let output = String::from_utf8(result.stdout).map_err(|e| format!("sha256sum output: {e}"))?;
    let hash = output.split_whitespace().next().ok_or("empty sha256sum output")?;
    if hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit()) { return Err("invalid sha256sum digest".into()); }
    Ok(Fingerprint { sha256: hash.to_ascii_lowercase(), size_bytes })
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct FileStamp { size_bytes: u64, modified_seconds: u64, modified_nanoseconds: u32 }

fn stamp(metadata: &Metadata) -> Result<FileStamp, String> {
    if !metadata.is_file() { return Err("input must be a regular file".into()); }
    let modified = metadata.modified().map_err(io_error)?.duration_since(UNIX_EPOCH)
        .map_err(|e| format!("file modification time: {e}"))?;
    Ok(FileStamp { size_bytes: metadata.len(), modified_seconds: modified.as_secs(), modified_nanoseconds: modified.subsec_nanos() })
}

#[derive(Serialize)]
struct Source {
    path: PathBuf,
    #[serde(flatten)]
    fingerprint: Fingerprint,
    metadata_before: FileStamp,
    metadata_after: Option<FileStamp>,
    sha256_after: Option<String>,
    records: u64, bases: u64, validated_eof: bool,
}

fn verified_source(path: &Path) -> Result<Source, String> {
    let path = fs::canonicalize(path).map_err(io_error)?;
    if path.to_str().is_none() { return Err("input path must be valid UTF-8".into()); }
    let before = stamp(&fs::metadata(&path).map_err(io_error)?)?;
    let fp = fingerprint(&path)?;
    let after = stamp(&fs::metadata(&path).map_err(io_error)?)?;
    if before != after || before.size_bytes != fp.size_bytes { return Err(format!("input changed while hashing: {}", path.display())); }
    Ok(Source { path, fingerprint: fp, metadata_before: before, metadata_after: None,
        sha256_after: None, records: 0, bases: 0, validated_eof: false })
}

fn verify_source_after(source: &mut Source) -> Result<(), String> {
    let before = stamp(&fs::metadata(&source.path).map_err(io_error)?)?;
    let fp = fingerprint(&source.path)?;
    let after = stamp(&fs::metadata(&source.path).map_err(io_error)?)?;
    if before != after || before != source.metadata_before || fp != source.fingerprint {
        return Err(format!("input changed during analysis: {}", source.path.display()));
    }
    source.metadata_after = Some(after);
    source.sha256_after = Some(fp.sha256);
    Ok(())
}

fn interval_union(arrays: &[Array]) -> Vec<(usize, usize)> {
    let mut rows: Vec<_> = arrays.iter().map(|x| (x.start0, x.end0)).collect();
    rows.sort_unstable();
    let mut union: Vec<(usize, usize)> = Vec::new();
    for (start, end) in rows {
        match union.last_mut() {
            Some(last) if start <= last.1 => last.1 = last.1.max(end),
            _ => union.push((start, end)),
        }
    }
    union
}

#[derive(Debug, Serialize)]
struct Anchor {
    id: String, start0: usize, end0: usize, length: usize,
    status: &'static str, stop_reason: &'static str,
}

fn flank(seq: &[u8], start: usize, end: usize, union: &[(usize, usize)], left: bool, enabled: bool, id: &str) -> Result<Anchor, String> {
    let edge = if left { start } else { end };
    if !enabled {
        return Ok(Anchor { id: String::new(), start0: edge, end0: edge, length: 0,
            status: "disabled", stop_reason: "export_not_requested" });
    }
    let index = union.partition_point(|&(a, _)| a <= start).checked_sub(1)
        .ok_or("array absent from interval union")?;
    let (merged_start, merged_end) = union[index];
    if merged_end < end { return Err("array not covered by interval union".into()); }
    let (mut a, mut b, mut stop) = if left {
        if merged_start < start { (start, start, "reported_array") }
        else {
            let hard = if index == 0 { 0 } else { union[index - 1].1 };
            let a = start.saturating_sub(1000).max(hard);
            (a, start, if a == hard { if index == 0 { "read_boundary" } else { "reported_array" } } else { "maximum_length" })
        }
    } else if merged_end > end { (end, end, "reported_array") }
    else {
        let hard = union.get(index + 1).map_or(seq.len(), |x| x.0);
        let b = end.saturating_add(1000).min(hard);
        (end, b, if b == hard { if index + 1 == union.len() { "read_boundary" } else { "reported_array" } } else { "maximum_length" })
    };
    if left {
        if let Some(p) = (a..b).rev().find(|&p| !motif::is_dna(&seq[p..p + 1])) { a = p + 1; stop = "ambiguous_base"; }
    } else if let Some(p) = (a..b).find(|&p| !motif::is_dna(&seq[p..p + 1])) { b = p; stop = "ambiguous_base"; }
    let length = b - a;
    Ok(Anchor { id: if length >= 200 { id.to_string() } else { String::new() },
        start0: a, end0: b, length, status: if length >= 200 { "exported" } else { "below_minimum" }, stop_reason: stop })
}

#[derive(Default, Serialize)]
struct GroupCounts {
    records_with_arrays: u64, arrays: u64, modeled_span_bp_sum: u64,
    exact_unit_bases_sum: u64, max_model_span_bp: usize,
    arrays_model_span_ge_500: u64, arrays_covering_ge_90_percent_of_record: u64,
    arrays_model_span_ge_500_and_covering_ge_90_percent_of_record: u64,
    orientation_ambiguous_arrays: u64, optimal_parse_ambiguous_arrays: u64,
}

#[derive(Default, Serialize)]
struct Totals {
    records: u64, bases: u64, arrays: u64, anchors: u64, records_with_arrays: u64,
    orientations_considered: u64, orientations_passing_literal_gate: u64,
    orientations_unscored_literal_gate: u64,
    max_model_span_bp: usize, arrays_model_span_ge_500: u64,
    arrays_covering_ge_90_percent_of_record: u64,
    arrays_model_span_ge_500_and_covering_ge_90_percent_of_record: u64,
    groups: BTreeMap<String, GroupCounts>,
}

fn full_record(span: usize, read_length: usize) -> bool {
    (span as u128) * 10 >= (read_length as u128) * 9
}

fn array_id(source_index: usize, ordinal: u64, array_index: usize) -> String {
    format!("s{source_index:03}_r{ordinal:012}_a{array_index:06}")
}

struct Outputs { arrays: BufWriter<File>, anchors: BufWriter<File> }

fn emit_record(record: &input::Record, result: &ScanResult, source_index: usize, ordinal: u64,
               enabled: bool, out: &mut Outputs, totals: &mut Totals) -> Result<(), String> {
    totals.records += 1;
    totals.bases += record.seq.len() as u64;
    totals.orientations_considered += result.orientations_considered as u64;
    totals.orientations_passing_literal_gate += result.orientations_passing_literal_gate as u64;
    totals.orientations_unscored_literal_gate += (result.orientations_considered - result.orientations_passing_literal_gate) as u64;
    if !result.arrays.is_empty() { totals.records_with_arrays += 1; }
    let union = interval_union(&result.arrays);
    let mut groups_seen = BTreeSet::new();
    for (index, array) in result.arrays.iter().enumerate() {
        let id = array_id(source_index, ordinal, index + 1);
        let left = flank(&record.seq, array.start0, array.end0, &union, true, enabled, &format!("{id}_L"))?;
        let right = flank(&record.seq, array.start0, array.end0, &union, false, enabled, &format!("{id}_R"))?;
        for anchor in [&left, &right] {
            if anchor.status == "exported" {
                writeln!(out.anchors, ">{}", anchor.id).map_err(io_error)?;
                out.anchors.write_all(&record.seq[anchor.start0..anchor.end0]).map_err(io_error)?;
                out.anchors.write_all(b"\n").map_err(io_error)?;
                totals.anchors += 1;
            }
        }
        let row = json!({"array_id":id, "source_index":source_index, "read_ordinal":ordinal,
            "read_id":record.id, "read_length":record.seq.len(), "group_id":array.group_id,
            "start0":array.start0, "end0":array.end0, "span_bp":array.span_bp,
            "strand":array.strand, "alignment_strand":array.alignment_strand,
            "left_end_distance":array.start0, "right_end_distance":record.seq.len() - array.end0,
            "model":array.model, "left_anchor":left, "right_anchor":right});
        serde_json::to_writer(&mut out.arrays, &row).map_err(|e| e.to_string())?;
        out.arrays.write_all(b"\n").map_err(io_error)?;
        let long = array.span_bp >= 500;
        let full = full_record(array.span_bp, record.seq.len());
        totals.arrays += 1;
        totals.max_model_span_bp = totals.max_model_span_bp.max(array.span_bp);
        totals.arrays_model_span_ge_500 += u64::from(long);
        totals.arrays_covering_ge_90_percent_of_record += u64::from(full);
        totals.arrays_model_span_ge_500_and_covering_ge_90_percent_of_record += u64::from(long && full);
        let group = totals.groups.get_mut(&array.group_id).ok_or("unknown emitted group")?;
        group.arrays += 1;
        group.modeled_span_bp_sum += array.span_bp as u64;
        group.exact_unit_bases_sum += array.model.exact_unit_bases as u64;
        group.max_model_span_bp = group.max_model_span_bp.max(array.span_bp);
        group.arrays_model_span_ge_500 += u64::from(long);
        group.arrays_covering_ge_90_percent_of_record += u64::from(full);
        group.arrays_model_span_ge_500_and_covering_ge_90_percent_of_record += u64::from(long && full);
        group.orientation_ambiguous_arrays += u64::from(array.strand == "ambiguous");
        group.optimal_parse_ambiguous_arrays += u64::from(array.model.optimal_parse_ambiguous);
        groups_seen.insert(&array.group_id);
    }
    for id in groups_seen { totals.groups.get_mut(id).ok_or("unknown group")?.records_with_arrays += 1; }
    Ok(())
}

fn scan_batch(scanner: &Scanner, records: &[input::Record], threads: usize) -> Result<Vec<ScanResult>, String> {
    if records.is_empty() { return Ok(Vec::new()); }
    let chunk_size = records.len().div_ceil(threads.min(records.len()));
    std::thread::scope(|scope| {
        let handles: Vec<_> = records.chunks(chunk_size).map(|chunk| {
            scope.spawn(move || chunk.iter().map(|record| scanner.scan(&record.seq)).collect::<Result<Vec<_>, _>>())
        }).collect();
        let mut results = Vec::with_capacity(records.len());
        for handle in handles {
            results.extend(handle.join().map_err(|_| "scanner worker panicked; no complete result")??);
        }
        Ok(results)
    })
}

fn run(args: Args) -> Result<(), String> {
    let started = Instant::now();
    let started_unix_seconds = SystemTime::now().duration_since(UNIX_EPOCH).map_err(|e| e.to_string())?.as_secs_f64();
    let memory_policy = memory::prepare(args.memory_limit_mib, args.threads, args.batch_bases)?;
    let mut definition_source = verified_source(&args.definitions)?;
    let definitions: Vec<Definition> = serde_json::from_reader(File::open(&definition_source.path).map_err(io_error)?)
        .map_err(|e| format!("definitions JSON: {e}"))?;
    verify_source_after(&mut definition_source)?;
    let scanner = Scanner::new(definitions, args.config)?;
    let mut sources = args.inputs.iter().map(|p| verified_source(p)).collect::<Result<Vec<_>, _>>()?;
    let mut unique = HashSet::new();
    if sources.iter().any(|s| !unique.insert(s.path.clone())) { return Err("duplicate canonical input path".into()); }
    let executable = fs::canonicalize(std::env::current_exe().map_err(io_error)?).map_err(io_error)?;
    let binary_fingerprint = fingerprint(&executable)?;
    let preflight_seconds = started.elapsed().as_secs_f64();
    fs::create_dir(&args.output).map_err(|e| format!("cannot create new output directory {}: {e}", args.output.display()))?;
    let arrays_partial = args.output.join("arrays.jsonl.partial");
    let anchors_partial = args.output.join("anchors.fasta.partial");
    let mut out = Outputs { arrays: create(&arrays_partial)?, anchors: create(&anchors_partial)? };
    let mut totals = Totals::default();
    for definition in &scanner.definitions { totals.groups.insert(definition.id.clone(), GroupCounts::default()); }
    let analysis_started = Instant::now();
    let mut sampled_prefix = false;
    for (source_zero, source) in sources.iter_mut().enumerate() {
        if args.max_records > 0 && totals.records >= args.max_records { sampled_prefix = true; break; }
        let mut reader = input::Reader::new(input::open(source.path.to_str().ok_or("non-UTF8 input path")?)?);
        loop {
            let mut batch = Vec::new();
            let mut batch_bases = 0;
            let first_ordinal = source.records + 1;
            while batch_bases < memory_policy.effective_batch_bases {
                if args.max_records > 0 && totals.records + batch.len() as u64 >= args.max_records {
                    sampled_prefix = true; break;
                }
                match reader.next_record().map_err(|e| format!("{}: {e}", source.path.display()))? {
                    None => { source.validated_eof = true; break; }
                    Some(record) => {
                        if record.seq.len() > scanner.config.max_record_bases {
                            return Err(format!("source {} record {} exceeds --max-record-bases; no complete result", source_zero + 1, source.records + batch.len() as u64 + 1));
                        }
                        batch_bases += record.seq.len();
                        batch.push(record);
                    }
                }
            }
            let results = scan_batch(&scanner, &batch, args.threads)?;
            if results.len() != batch.len() { return Err("scanner batch cardinality mismatch".into()); }
            for (index, (record, result)) in batch.iter().zip(&results).enumerate() {
                emit_record(record, result, source_zero + 1, first_ordinal + index as u64, args.export_anchors, &mut out, &mut totals)?;
            }
            source.records += batch.len() as u64;
            source.bases += batch_bases as u64;
            if source.validated_eof || sampled_prefix { break; }
        }
        if sampled_prefix { break; }
    }
    if totals.records == 0 { return Err("no input sequence records; no complete result".into()); }
    let analysis_seconds = analysis_started.elapsed().as_secs_f64();
    out.arrays.flush().map_err(io_error)?;
    out.anchors.flush().map_err(io_error)?;
    out.arrays.get_ref().sync_all().map_err(io_error)?;
    out.anchors.get_ref().sync_all().map_err(io_error)?;
    drop(out);
    let verification_started = Instant::now();
    for source in &mut sources { verify_source_after(source)?; }
    verify_source_after(&mut definition_source)?;
    if fingerprint(&executable)? != binary_fingerprint { return Err("binary changed during analysis".into()); }
    let validated_eof = sources.iter().all(|s| s.validated_eof);
    let complete_input = validated_eof && !sampled_prefix;
    if !complete_input && !sampled_prefix { return Err("requested full input did not reach validated EOF".into()); }
    let arrays = args.output.join("arrays.jsonl");
    let anchors = args.output.join("anchors.fasta");
    fs::rename(&arrays_partial, &arrays).map_err(io_error)?;
    fs::rename(&anchors_partial, &anchors).map_err(io_error)?;
    let array_fingerprint = fingerprint(&arrays)?;
    let anchor_fingerprint = fingerprint(&anchors)?;
    let verification_seconds = verification_started.elapsed().as_secs_f64();
    let summary = json!({
        "schema_version":1, "program":"scan_irregular_arrays", "model_version":"0.1.0",
        "native_package_version":env!("CARGO_PKG_VERSION"), "status":"complete",
        "complete_input":complete_input, "validated_eof":validated_eof, "sampled_prefix":sampled_prefix,
        "validated_requested_scope":true, "export_truncated":false, "resource_limit_hit":false,
        "inputs":sources, "definitions_source":definition_source, "definitions":scanner.definitions,
        "config":scanner.config, "threads":args.threads, "max_records":args.max_records,
        "export_anchors":args.export_anchors, "anchor_min_bp":200, "anchor_max_bp":1000,
        "memory_policy":memory_policy, "records":totals.records, "bases":totals.bases,
        "arrays":totals.arrays, "anchors":totals.anchors, "counts":totals,
        "binary":{"path":executable, "sha256":binary_fingerprint.sha256,"size_bytes":binary_fingerprint.size_bytes},
        "outputs":{"arrays.jsonl":array_fingerprint,"anchors.fasta":anchor_fingerprint},
        "timing":{"started_unix_seconds":started_unix_seconds,"wall_seconds":started.elapsed().as_secs_f64(),
            "preflight_seconds":preflight_seconds,"analysis_seconds":analysis_seconds,
            "postflight_verification_seconds":verification_seconds,"rss":"measure externally; memory limit is virtual address space"},
        "metric_definitions":{
            "modeled_span":"read-coordinate span explained by the supplied edit model; not a homogeneous or exact repeat length",
            "exact_unit_bases":"sum of exact token bases in one representative parse; not a contiguous exact tract",
            "full_record":"one array span covers >=90% of its source record; not proof of a physical chromosome end",
            "orientation_gate":"record/group/orientation with >=min_exact_anchors literal unit hits; unscored is not an absence call",
            "sums":"array sums can count overlapping model groups/orientations; not union-covered genome bases",
            "anchors":"adjacent original-strand normalized input bases, stopping at every reported array union or ambiguity",
            "scope":"optional model evidence only; no TRM/catalog/adequacy/absence promotion"
        }
    });
    let partial = args.output.join("summary.json.partial");
    let mut writer = create(&partial)?;
    serde_json::to_writer_pretty(&mut writer, &summary).map_err(|e| e.to_string())?;
    writer.write_all(b"\n").map_err(io_error)?;
    writer.flush().map_err(io_error)?;
    writer.get_ref().sync_all().map_err(io_error)?;
    drop(writer);
    fs::rename(partial, args.output.join("summary.json")).map_err(io_error)?;
    eprintln!("Completed {} records, {} bases, {} model arrays, {} anchors; complete_input={complete_input}", totals.records, totals.bases, totals.arrays, totals.anchors);
    Ok(())
}

fn main() {
    if let Err(error) = parse_arguments(std::env::args().skip(1)).and_then(|args| args.map_or(Ok(()), run)) {
        eprintln!("scan_irregular_arrays: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod exporter_tests {
    use super::*;

    #[test]
    fn defaults_and_explicit_prefix() {
        let args = parse_arguments(["a", "b", "--definitions", "d", "--output-dir", "o", "--max-records", "12", "--export-anchors"].map(str::to_string)).unwrap().unwrap();
        assert_eq!(args.inputs.len(), 2);
        assert_eq!(args.threads, 2);
        assert_eq!(args.max_records, 12);
        assert!(args.export_anchors);
        assert_eq!(args.config.min_span, 100);
    }

    #[test]
    fn invalid_options_rejected() {
        for extra in [vec!["--threads", "0"], vec!["--threads", "2", "--threads", "3"], vec!["--batch-bases", "0"], vec!["--memory-limit-mib", "0"], vec!["--export-anchors", "--export-anchors"], vec!["--max-records", "-1"]] {
            let mut args = vec!["a", "--definitions", "d", "--output-dir", "o"];
            args.extend(extra);
            assert!(parse_arguments(args.into_iter().map(str::to_string)).is_err());
        }
    }

    #[test]
    fn stable_numeric_ids_and_full_record_threshold() {
        assert_eq!(array_id(2, 9, 3), "s002_r000000000009_a000003");
        assert!(full_record(900, 1000));
        assert!(!full_record(899, 1000));
        assert!(full_record(usize::MAX, usize::MAX));
    }
}
