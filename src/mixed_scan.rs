//! Streaming, explicitly targeted mixed-array measurements. No biological call.
use crate::mixed::{MixedDefinition, MixedScanner, MixedTract};
use crate::{evidence, input, mixed_capture, Options, Source, Summary};
use serde::Serialize;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Arc, Mutex};
use std::time::{Instant, UNIX_EPOCH};

#[derive(Clone, Serialize)]
struct ArrayExample {
    input_index: usize,
    record_id: String,
    record_length: usize,
    #[serde(flatten)]
    tract: MixedTract,
}

#[derive(Default, Serialize)]
struct Group {
    name: String,
    units: Vec<String>,
    records: u64,
    mixed_records: u64,
    pure_only_records: u64,
    nonqualifying_records: u64,
    total_mixed_arrays: u64,
    full_mixed_records: u64,
    entire_mixed_records: u64,
    repeat_bases: u64,
    max_array_bp: usize,
    max_any_array_bp: usize,
    long_array_records: u64,
    segmentation_ambiguous_records: u64,
    ambiguous_orientation_records: u64,
    g_right_boundary_records: u64,
    c_left_boundary_records: u64,
    g_left_boundary_records: u64,
    c_right_boundary_records: u64,
    ambiguous_boundary_records: u64,
    internal_array_records: u64,
    clean_two_flank_records: u64,
    left_g_end_records: u64,
    left_c_end_records: u64,
    right_g_end_records: u64,
    right_c_end_records: u64,
    long_expected_end_records: u64,
    examples: Vec<ArrayExample>,
}

impl Group {
    fn new(definition: &MixedDefinition) -> Self {
        Self {
            name: definition.id.clone(),
            units: definition.units.clone(),
            ..Self::default()
        }
    }

    fn trim_examples(&mut self, limit: usize) {
        self.examples.sort_by(|a, b| {
            (a.input_index, &a.record_id, a.record_length, &a.tract).cmp(&(
                b.input_index,
                &b.record_id,
                b.record_length,
                &b.tract,
            ))
        });
        self.examples.truncate(limit);
    }

    fn observe(
        &mut self,
        all: &[&MixedTract],
        record: &input::Record,
        input_index: usize,
        opts: &Options,
    ) {
        if all.is_empty() {
            return;
        }
        self.records += 1;
        self.max_any_array_bp = self
            .max_any_array_bp
            .max(all.iter().map(|t| t.end - t.start).max().unwrap_or(0));
        self.segmentation_ambiguous_records +=
            u64::from(all.iter().any(|t| t.segmentation_ambiguous));
        self.ambiguous_orientation_records += u64::from(all.iter().any(|t| t.strand.is_none()));
        let mut mixed: Vec<&MixedTract> = all.iter().copied().filter(|t| t.is_mixed).collect();
        if mixed.is_empty() {
            if all.iter().any(|t| t.composition == "pure") {
                self.pure_only_records += 1;
            } else {
                self.nonqualifying_records += 1;
            }
            return;
        }
        mixed.sort_unstable();
        let mut blockers = all.to_vec();
        blockers.sort_by_key(|t| t.start);
        let mut max_end = 0;
        let prefix_ends: Vec<usize> = blockers
            .iter()
            .map(|t| {
                max_end = max_end.max(t.end);
                max_end
            })
            .collect();
        let overlaps = |start: usize, end: usize| {
            let before = blockers.partition_point(|t| t.start < end);
            before > 0 && prefix_ends[before - 1] > start
        };
        self.mixed_records += 1;
        self.total_mixed_arrays += mixed.len() as u64;
        let n = record.seq.len();
        let full = mixed
            .iter()
            .any(|t| (t.end - t.start) as f64 >= n as f64 * 0.9);
        self.full_mixed_records += u64::from(full);
        self.entire_mixed_records += u64::from(mixed.iter().any(|t| t.start == 0 && t.end == n));
        self.long_array_records +=
            u64::from(mixed.iter().any(|t| t.end - t.start >= opts.cfg.min_tract));
        self.internal_array_records += u64::from(
            mixed
                .iter()
                .any(|t| t.start > opts.evidence.end_slack && n - t.end > opts.evidence.end_slack),
        );
        let mut boundaries = [false; 6]; // left G/C/?; right G/C/?
        let mut ends = [false; 4]; // left G/C; right G/C
        let mut long_expected = false;
        let mut two_flanks = false;
        let mut covered_end = 0;
        let flank = opts.evidence.flank_bases;
        for (tract_index, t) in mixed.iter().enumerate() {
            let span = t.end - t.start;
            self.repeat_bases += t.end.saturating_sub(covered_end.max(t.start)) as u64;
            covered_end = covered_end.max(t.end);
            self.max_array_bp = self.max_array_bp.max(span);
            let strand = match t.strand {
                Some(true) => 0,
                Some(false) => 1,
                None => 2,
            };
            if strand < 2 {
                if t.start <= opts.cfg.end_window {
                    ends[strand] = true;
                }
                if n - t.end <= opts.cfg.end_window {
                    ends[2 + strand] = true;
                }
                long_expected |= span >= opts.cfg.min_tract
                    && ((strand == 1 && t.start <= opts.cfg.end_window)
                        || (strand == 0 && n - t.end <= opts.cfg.end_window));
            }
            let clean_left = t.start >= flank
                && evidence::usable_flank(&record.seq[t.start - flank..t.start])
                && !overlaps(t.start - flank, t.start);
            let clean_right = n - t.end >= flank
                && evidence::usable_flank(&record.seq[t.end..t.end + flank])
                && !overlaps(t.end, t.end + flank);
            two_flanks |= clean_left && clean_right;
            if !full {
                if t.start <= opts.evidence.end_slack && clean_right {
                    boundaries[strand] = true;
                }
                if n - t.end <= opts.evidence.end_slack && clean_left {
                    boundaries[3 + strand] = true;
                }
            }
            if tract_index < opts.mixed_max_examples {
                self.examples.push(ArrayExample {
                    input_index,
                    record_id: record.id.clone(),
                    record_length: n,
                    tract: (*t).clone(),
                });
            }
        }
        self.g_left_boundary_records += u64::from(boundaries[0]);
        self.c_left_boundary_records += u64::from(boundaries[1]);
        self.g_right_boundary_records += u64::from(boundaries[3]);
        self.c_right_boundary_records += u64::from(boundaries[4]);
        self.ambiguous_boundary_records += u64::from(boundaries[2] || boundaries[5]);
        self.left_g_end_records += u64::from(ends[0]);
        self.left_c_end_records += u64::from(ends[1]);
        self.right_g_end_records += u64::from(ends[2]);
        self.right_c_end_records += u64::from(ends[3]);
        self.long_expected_end_records += u64::from(long_expected);
        self.clean_two_flank_records += u64::from(two_flanks);
        self.trim_examples(opts.mixed_max_examples);
    }

    fn merge(&mut self, other: Self, example_limit: usize) {
        macro_rules! add { ($($field:ident),+) => { $(self.$field += other.$field;)+ }; }
        add!(
            records,
            mixed_records,
            pure_only_records,
            nonqualifying_records,
            total_mixed_arrays,
            full_mixed_records,
            entire_mixed_records,
            repeat_bases,
            long_array_records,
            segmentation_ambiguous_records,
            ambiguous_orientation_records,
            g_right_boundary_records,
            c_left_boundary_records,
            g_left_boundary_records,
            c_right_boundary_records,
            ambiguous_boundary_records,
            internal_array_records,
            clean_two_flank_records,
            left_g_end_records,
            left_c_end_records,
            right_g_end_records,
            right_c_end_records,
            long_expected_end_records
        );
        self.max_array_bp = self.max_array_bp.max(other.max_array_bp);
        self.max_any_array_bp = self.max_any_array_bp.max(other.max_any_array_bp);
        self.examples.extend(other.examples);
        self.trim_examples(example_limit);
    }
}

struct WorkerResult {
    groups: Vec<Group>,
    qc: evidence::QualitySummary,
    capture: Option<mixed_capture::Buffer>,
}

type Batch = Vec<(usize, u64, input::Record)>;

fn worker(
    rx: Arc<Mutex<mpsc::Receiver<Batch>>>,
    index: Arc<MixedScanner>,
    opts: Options,
    cancelled: Arc<AtomicBool>,
) -> Result<WorkerResult, String> {
    let mut groups: Vec<Group> = index.definitions.iter().map(Group::new).collect();
    let mut qc = evidence::QualitySummary::default();
    let mut capture = opts
        .mixed_capture
        .enabled
        .then(|| mixed_capture::Buffer::new(opts.mixed_capture.clone()));
    while !cancelled.load(Ordering::Relaxed) {
        let Ok(batch) = rx.lock().map_err(|_| "worker queue lock failed")?.recv() else {
            break;
        };
        for (input_index, record_ordinal, record) in batch {
            if cancelled.load(Ordering::Relaxed) {
                return Ok(WorkerResult {
                    groups,
                    qc,
                    capture,
                });
            }
            qc.observe(&record.seq);
            let hits = index.scan(&record.seq);
            for group in &mut groups {
                let tracts: Vec<_> = hits
                    .iter()
                    .filter(|t| t.definition_id == group.name)
                    .collect();
                group.observe(&tracts, &record, input_index, &opts);
            }
            if let Some(buffer) = &mut capture {
                buffer.observe(
                    input_index,
                    record_ordinal,
                    record,
                    hits,
                    &index.definitions,
                );
            }
        }
    }
    Ok(WorkerResult {
        groups,
        qc,
        capture,
    })
}

pub(crate) fn run(opts: Options) -> Result<(), String> {
    let definitions_path = opts
        .mixed_definitions
        .as_ref()
        .ok_or("mixed definitions missing")?;
    let definitions: Vec<MixedDefinition> = serde_json::from_reader(BufReader::new(
        File::open(definitions_path).map_err(|e| format!("{definitions_path}: {e}"))?,
    ))
    .map_err(|e| format!("{definitions_path}: {e}"))?;
    let index = Arc::new(MixedScanner::new(definitions, opts.mixed_config.clone())?);
    let started = Instant::now();
    let mut sources = Vec::new();
    for path in &opts.paths {
        let metadata = if path == "-" {
            None
        } else {
            Some(fs::metadata(path).map_err(|e| format!("{path}: {e}"))?)
        };
        sources.push(Source {
            path: path.clone(),
            bytes: metadata.as_ref().map(|m| m.len()),
            modified_unix_seconds: metadata
                .and_then(|m| m.modified().ok())
                .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
                .map(|d| d.as_secs()),
        });
    }
    let (tx, rx) = mpsc::sync_channel::<Batch>(opts.threads * 2);
    let rx = Arc::new(Mutex::new(rx));
    let cancelled = Arc::new(AtomicBool::new(false));
    let mut handles = Vec::new();
    for _ in 0..opts.threads {
        let receiver = Arc::clone(&rx);
        let scanner = Arc::clone(&index);
        let stop = Arc::clone(&cancelled);
        let worker_options = opts.clone();
        handles.push(
            std::thread::Builder::new()
                .spawn(move || worker(receiver, scanner, worker_options, stop))
                .map_err(|e| {
                    format!("cannot start mixed worker (check memory/thread limits): {e}")
                })?,
        );
    }
    drop(rx);
    let mut summary = Summary::default();
    let parse_result: Result<(), String> = (|| {
        let mut batch = Vec::new();
        let mut batch_bases = 0;
        'files: for (input_index, path) in opts.paths.iter().enumerate() {
            let source = input::open(path)?;
            let mut reader = if opts.mixed_capture.enabled {
                input::Reader::with_quality(source, opts.mixed_capture.max_record_bases)
            } else {
                input::Reader::new(source)
            };
            let mut record_ordinal = 0_u64;
            loop {
                if cancelled.load(Ordering::Relaxed) {
                    return Err("mixed analysis cancelled".into());
                }
                if opts.max_records > 0 && summary.records >= opts.max_records {
                    summary.sampled_prefix = true;
                    break 'files;
                }
                let Some(record) = reader.next_record().map_err(|e| format!("{path}: {e}"))? else {
                    break;
                };
                summary.records += 1;
                record_ordinal += 1;
                summary.bases += record.seq.len() as u64;
                summary.max_record_bases = summary.max_record_bases.max(record.seq.len());
                batch_bases += record.seq.len();
                batch.push((input_index, record_ordinal, record));
                if batch_bases >= opts.batch_bases {
                    tx.send(std::mem::take(&mut batch))
                        .map_err(|_| "all mixed workers stopped")?;
                    batch_bases = 0;
                }
            }
        }
        if !batch.is_empty() {
            tx.send(batch).map_err(|_| "all mixed workers stopped")?;
        }
        Ok(())
    })();
    if parse_result.is_err() {
        cancelled.store(true, Ordering::Relaxed);
    }
    drop(tx);
    let mut groups: Vec<Group> = index.definitions.iter().map(Group::new).collect();
    let mut qc = evidence::QualitySummary::default();
    let mut worker_error = None;
    let mut capture = opts
        .mixed_capture
        .enabled
        .then(|| mixed_capture::Buffer::new(opts.mixed_capture.clone()));
    for handle in handles {
        match handle
            .join()
            .map_err(|_| "mixed worker panicked".to_string())
            .and_then(|r| r)
        {
            Ok(part) => {
                qc.merge(part.qc);
                for (combined, local) in groups.iter_mut().zip(part.groups) {
                    combined.merge(local, opts.mixed_max_examples);
                }
                if let (Some(combined), Some(local)) = (&mut capture, part.capture) {
                    combined.merge(local);
                }
            }
            Err(error) => worker_error = Some(error),
        }
    }
    if let Some(error) = worker_error {
        return Err(error);
    }
    parse_result?;
    if summary.records == 0 {
        return Err("no sequence records found".into());
    }
    summary.call = "mixed_array_sequence_evidence_only".into();
    summary.elapsed_seconds = started.elapsed().as_secs_f64();
    summary.megabases_per_second = summary.bases as f64 / 1e6 / summary.elapsed_seconds;
    let mut report = serde_json::json!({
        "schema_version": 1, "program": "trm", "version": env!("CARGO_PKG_VERSION"),
        "analysis_kind": "targeted_mixed_arrays", "command": std::env::args().collect::<Vec<_>>(),
        "inputs": sources, "threads": opts.threads, "summary": summary, "evidence_qc": qc.json(),
        "memory_policy": opts.memory_policy,
        "definitions_source": definitions_path,
        "mixed_arrays": {
            "schema_version": 1, "mode": opts.cfg.mode,
            "settings": {"min_array_bp": index.config.min_span,
                "min_copies_per_unit": index.config.min_component_copies,
                "min_distinct_units": 2, "long_array_bp": opts.cfg.min_tract,
                "full_fraction": 0.9, "boundary_flank_bp": opts.evidence.flank_bases,
                "boundary_end_slack_bp": opts.evidence.end_slack, "end_window_bp": opts.cfg.end_window,
                "max_examples": opts.mixed_max_examples, "exact_only": true,
                "partials_at_record_edges": index.config.allow_edge_partials},
            "definitions": index.definitions,
            "groups": groups,
            "method": {
                "grammar": "Literal supplied G-strand primitive units or their reverse complements; exact contiguous parses only, no distance joining or homogeneous seed.",
                "coordinates": "Zero-based half-open start/end include compatible partial units at record edges only; core_start/core_end contain complete units only; when core_bounds_ambiguous, they represent one deterministic valid parse, while unit-count bounds span all tied parses and need not sum to that representative core length.",
                "counts": "unit_counts and unit_counts_max are marginal lower/upper bounds over tied parses, not a jointly attained composition; mixed requires >=2 guaranteed component counts above threshold.",
                "full_repeat": "Any one mixed array covers >=90% of the record; full records receive no clean-junction credit. entire_mixed_records separately requires 100%.",
                "orientation": "G equals the supplied unit spelling; expected geometry is C-left or G-right. Ambiguous orientation receives no directional credit.",
                "flank": "At least boundary_flank_bp contiguous ACGT, >=3 distinct bases, no >=85% lag agreement at periods 1..min(32,flank/3), and no overlapping reported same-group array.",
                "pure_only_records": "No qualifying mixed array and at least one reported pure tract; other subthreshold/ambiguous tracts may coexist. nonqualifying_records has neither qualifying mixed nor pure tract.",
                "record_counts": "Group records and direction quadrants can overlap. Mates, sequence duplicates, and contigs are not independent molecules. repeat_bases is the within-group interval union of mixed arrays."
            }
        },
        "warnings": [
            "This is an explicit unit-grammar structural test, not de novo motif discovery or a telomere call.",
            "Exact matching may split or miss arrays with substitutions, indels, unlisted units or insufficient complete copies; no hit is not absence or ALT.",
            "Record-edge clipping and read-boundary geometry do not establish physical chromosome termini; internal arrays remain sequence positives.",
            "Existing homogeneous motif evidence, classifications and data-adequacy estimates are not replaced.",
            "Definitions must supply compatible G-strand phases; independent unit rotations change the allowed junction grammar."
        ]
    });
    if let Some(buffer) = capture {
        report["mixed_arrays"]["capture"] =
            buffer.json(summary.sampled_prefix, opts.threads, opts.batch_bases);
    }
    write(&opts.output, &report)?;
    eprintln!(
        "{} records, {:.2} Mb, {} mixed grammar groups; {:.3}s ({:.2} Mb/s)",
        summary.records,
        summary.bases as f64 / 1e6,
        index.definitions.len(),
        summary.elapsed_seconds,
        summary.megabases_per_second
    );
    Ok(())
}

fn write(prefix: &str, report: &serde_json::Value) -> Result<(), String> {
    if let Some(parent) = std::path::Path::new(prefix)
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let json_path = format!("{prefix}.json");
    let tsv_path = format!("{prefix}.tsv");
    let suffix = format!(".{}.part", std::process::id());
    let json_temp = format!("{json_path}{suffix}");
    let tsv_temp = format!("{tsv_path}{suffix}");
    let mut json_created = false;
    let mut tsv_created = false;
    let result = (|| -> std::io::Result<()> {
        let json_file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&json_temp)?;
        json_created = true;
        let mut json = BufWriter::new(json_file);
        serde_json::to_writer_pretty(&mut json, report)?;
        writeln!(json)?;
        json.flush()?;
        let tsv_file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&tsv_temp)?;
        tsv_created = true;
        let mut tsv = BufWriter::new(tsv_file);
        let fields = [
            "name",
            "units",
            "records",
            "mixed_records",
            "pure_only_records",
            "full_mixed_records",
            "entire_mixed_records",
            "max_array_bp",
            "max_any_array_bp",
            "long_array_records",
            "repeat_bases",
            "c_left_boundary_records",
            "g_right_boundary_records",
            "g_left_boundary_records",
            "c_right_boundary_records",
            "internal_array_records",
            "clean_two_flank_records",
            "segmentation_ambiguous_records",
            "ambiguous_orientation_records",
        ];
        writeln!(tsv, "{}", fields.join("\t"))?;
        for group in report["mixed_arrays"]["groups"]
            .as_array()
            .into_iter()
            .flatten()
        {
            let cells: Vec<String> = fields
                .iter()
                .map(|f| match &group[*f] {
                    serde_json::Value::String(s) => s.replace(['\t', '\n', '\r'], " "),
                    serde_json::Value::Array(a) => a
                        .iter()
                        .filter_map(|v| v.as_str())
                        .collect::<Vec<_>>()
                        .join(","),
                    value => value.to_string(),
                })
                .collect();
            writeln!(tsv, "{}", cells.join("\t"))?;
        }
        tsv.flush()?;
        fs::rename(&tsv_temp, &tsv_path)?;
        fs::rename(&json_temp, &json_path)?;
        Ok(())
    })();
    if result.is_err() {
        if json_created {
            let _ = fs::remove_file(json_temp);
        }
        if tsv_created {
            let _ = fs::remove_file(tsv_temp);
        }
    }
    result.map_err(|e| format!("writing mixed output: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{detect, mixed, prior, support};

    fn options() -> Options {
        Options {
            action: "mixed".into(),
            prior: prior::PriorOptions::default(),
            support: support::SupportOptions::default(),
            cfg: detect::Config::default(),
            paths: vec![],
            output: String::new(),
            threads: 1,
            max_records: 0,
            batch_bases: 1000,
            memory_limit_mib: None,
            memory_policy: None,
            paired_layout: None,
            targets: None,
            evidence: evidence::Options::default(),
            mixed_definitions: Some("test.json".into()),
            mixed_config: mixed::MixedConfig::default(),
            mixed_max_examples: 20,
            mixed_capture: mixed_capture::Options::default(),
            performance_profile: false,
            quality_diagnostics: false,
            quality_max_record_bases: 100_000,
        }
    }

    fn measure(sequences: &[Vec<u8>]) -> Group {
        let opts = options();
        let definition = MixedDefinition {
            id: "test".into(),
            units: vec!["TTGGGTTAC".into(), "ACTTGGGTTAC".into()],
        };
        let scanner = MixedScanner::new(vec![definition], opts.mixed_config.clone()).unwrap();
        let mut group = Group::new(&scanner.definitions[0]);
        for (i, seq) in sequences.iter().enumerate() {
            let record = input::Record {
                id: format!("r{i}"),
                seq: seq.clone(),
                quality: None,
            };
            let hits = scanner.scan(seq);
            group.observe(&hits.iter().collect::<Vec<_>>(), &record, 0, &opts);
        }
        group
    }

    #[test]
    fn full_mixed_reads_do_not_get_boundary_credit() {
        let array = b"TTGGGTTACACTTGGGTTAC".repeat(8);
        let group = measure(&[array.clone(), crate::motif::reverse_complement(&array)]);
        assert_eq!(group.mixed_records, 2);
        assert_eq!(group.full_mixed_records, 2);
        assert_eq!(group.entire_mixed_records, 2);
        assert_eq!(
            group.g_right_boundary_records + group.c_left_boundary_records,
            0
        );
        assert_eq!(
            group.g_left_boundary_records + group.c_right_boundary_records,
            0
        );
    }

    #[test]
    fn clean_boundary_polarity_and_internal_counterevidence_are_separate() {
        let flank = b"GACTATCGATGACCTAGTCACTACGTCAATGC";
        assert!(evidence::usable_flank(&flank[..20]));
        assert!(evidence::usable_flank(&flank[flank.len() - 20..]));
        let array = b"TTGGGTTACACTTGGGTTAC".repeat(5);
        let right_g = [flank.as_slice(), array.as_slice()].concat();
        let left_c = crate::motif::reverse_complement(&right_g);
        let left_g = [array.as_slice(), flank.as_slice()].concat();
        let right_c = crate::motif::reverse_complement(&left_g);
        let internal = [flank.as_slice(), array.as_slice(), flank.as_slice()].concat();
        let group = measure(&[right_g, left_c, left_g, right_c, internal]);
        assert_eq!(group.mixed_records, 5);
        assert_eq!(group.full_mixed_records, 0);
        assert_eq!(group.g_right_boundary_records, 1);
        assert_eq!(group.c_left_boundary_records, 1);
        assert_eq!(group.g_left_boundary_records, 1);
        assert_eq!(group.c_right_boundary_records, 1);
        assert_eq!(group.clean_two_flank_records, 1);
        assert_eq!(group.internal_array_records, 1);
    }

    #[test]
    fn separated_arrays_and_separate_records_do_not_make_a_long_tract() {
        let array = b"TTGGGTTACACTTGGGTTAC".repeat(15);
        let separated = [array.as_slice(), b"NNNN", array.as_slice()].concat();
        let group = measure(&[separated, array.clone()]);
        assert_eq!(group.max_array_bp, 300);
        assert_eq!(group.repeat_bases, 900);
        assert_eq!(group.total_mixed_arrays, 3);
        assert_eq!(group.long_array_records, 0);
        assert_eq!(group.full_mixed_records, 1);
    }

    #[test]
    fn preexisting_temporary_path_is_neither_overwritten_nor_removed() {
        let directory = std::env::temp_dir().join(format!("trm-mixed-temp-{}", std::process::id()));
        fs::create_dir_all(&directory).unwrap();
        let prefix = directory.join("report");
        let collision = format!("{}.json.{}.part", prefix.display(), std::process::id());
        fs::write(&collision, b"preserve input").unwrap();
        assert!(write(prefix.to_str().unwrap(), &serde_json::json!({})).is_err());
        assert_eq!(fs::read(&collision).unwrap(), b"preserve input");
        assert!(!prefix.with_extension("json").exists());
        fs::remove_file(collision).unwrap();
        fs::remove_dir(directory).unwrap();
    }
}
