mod candidate;
mod detect;
mod evidence;
mod flank_period;
mod input;
mod memory;
mod mixed;
mod mixed_capture;
mod mixed_scan;
mod motif;
mod paired;
mod prior;
mod quality;
mod support;

use detect::{Config, Tract};
use input::profile;
use serde::Serialize;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Arc, Mutex};
use std::time::{Instant, UNIX_EPOCH};

const HELP: &str = concat!(
    "trm ",
    env!("CARGO_PKG_VERSION"),
    " — WGS telomeric motif candidate discovery

Usage: trm scan INPUT [INPUT ...] --output PREFIX [OPTIONS]
       trm evidence INPUT [INPUT ...] --targets TARGETS.json --output PREFIX [OPTIONS]
       trm mixed INPUT [INPUT ...] --mixed-definitions GROUPS.json --output PREFIX [OPTIONS]
       trm reclassify RESULT.json --output NEW_PREFIX [PRIOR OPTIONS]
       trm --help | --version

Inputs: FASTA / FASTQ, optionally gzip (magic detected); '-' reads stdin.
All inputs in one invocation MUST come from the same sample/species.

  --threads N          Worker threads [available CPUs, capped at 8]
  --min-k N            Minimum primitive period [3]
  --max-k N            Maximum primitive period, <=64 [32]
  --seed-bases N       Minimum exact tandem seed length [24]
  --min-copies N       Minimum exact seed copies, >=2 [3]
  --min-identity F     Substitution-only extension identity, 0.5..1 [0.85]
  --min-tract N        Single-record long-tract threshold in bp [500]
  --min-reads N        Minimum supporting records to report [2]
  --end-window N       Record-end proximity in bp [50]
  --mode reads|assembly   Assembly end proximity is evidence, not proof of termini
  --max-records N      Analyze first N records across inputs (0 = all) [0]
  --max-candidates N   Fail explicitly above this per-worker/global limit [100000]
  --batch-bases N      Approximate bases per work batch [4000000]
  --memory-limit-mib N Linux virtual-address-space cap, NOT an RSS budget;
                       reduces batch size if needed; allocation failure is fatal
  --paired-layout fr  Scan/evidence: exactly two synchronized FASTQ mates;
                       validate IDs and annotate pair orientation (opt in)
  --output PREFIX     Writes PREFIX.json, PREFIX.tsv, PREFIX.rejected.tsv
  --junction-flank-bases N  Minimum clean, non-short-period flank [20; >=20]
  --junction-end-slack N    Maximum bases outside a repeat at read boundary [3]
  --flank-period-annotation on|off  Describe candidate-period flank signals;
                       scan/evidence only, no change to calls or filtering [off]
  --performance-profile on|off  Buffer/batch producer wall-clock attribution;
                       scan/evidence only; no per-record timers [off]
  --quality-diagnostics off|phred33  Explicitly declare FASTQ encoding and
                       summarize whole-read/union-tract quality; no filter [off]
  --quality-max-record-bases N  Retain complete quality only up to this length
                       when diagnostics enabled [100000; <=1000000]
  --targets JSON      Evidence mode: motifs array or saved primary candidates
Evidence searches only requested families across full reads; zero-hit targets
are retained. It adds boundary evidence without replacing saved discovery data.

Mixed-array options (mixed only; exact, explicitly supplied unit grammar):
  --mixed-definitions JSON   Array of {id, units:[same-G-strand primitive DNA]}
  --mixed-min-span N         Minimum contiguous exact span [24]
  --mixed-min-component-copies N  Guaranteed copies of each of >=2 units [2]
  --mixed-max-examples N     Bounded deterministic examples per group [20]
  --mixed-capture on|off     Capture positive mixed SEQ/QUAL in JSON [off]
  --mixed-capture-max-records N  Global captured record cap [1000; <=1000000]
  --mixed-capture-max-bases N    Captured sequence base cap [1000000; <=1e9]
  --mixed-capture-max-record-bases N  Omit whole larger records [100000; <=1e7]
Capture retains the first eligible positive prefix in input order, validates the
entire analyzed stream, and never captures mates automatically. SEQ is uppercase;
FASTA quality is absent. Capture limits require --mixed-capture on.
Mixed mode does not use homogeneous tandem seeds or infer a novel unit grammar.
Pure-component arrays are recorded separately; sequence evidence is not a TRM call.

Prior options (scan and reclassify):
  --prior-policy canonical|telomere|discovery  Reference / variant / all candidates [canonical]
  --prior-max-edits N     Maximum circular/strand-invariant edit distance [2]
  --prior-max-distance F  Also require edits / max(primitive lengths) <= F [0.34]
  --related-motif DNA     Known motif of a close relative; repeat for alternatives
  --lineage TEXT          Label for explicitly supplied related motifs

Background support (scan and reclassify; opt in with an explicit model):
  --background-model JSON  Negative-control profile from calibrate_support.py
  --support-policy filter|annotate  Filter primary candidates or annotate [annotate]
  --support-min-enrichment F  Minimum read-bound and base-burden enrichment [2]
Fractions normalize sequencing depth, not genomic coverage. Background scores
are screening heuristics, not probabilities/FDR or proof of terminal location.

Canonical policy accepts exact established or explicitly supplied families;
other repeats remain available as raw evidence for novel-candidate assessment.
Legacy telomere policy permits reference edit neighbours and orders prior tiers.
Reclassify reads saved JSON evidence only; original scan timings stay separate.

Candidates are tandem repeats, not confirmed telomeres. Short reads cannot
establish a contiguous 500 bp tract. No candidate is not evidence of ALT.
Indel-rich reads and complex/irregular repeats need independent validation.
"
);

#[derive(Clone)]
struct Options {
    action: String,
    prior: prior::PriorOptions,
    support: support::SupportOptions,
    cfg: Config,
    paths: Vec<String>,
    output: String,
    threads: usize,
    max_records: u64,
    batch_bases: usize,
    memory_limit_mib: Option<u64>,
    memory_policy: Option<memory::MemoryPolicy>,
    paired_layout: Option<String>,
    targets: Option<String>,
    evidence: evidence::Options,
    mixed_definitions: Option<String>,
    mixed_config: mixed::MixedConfig,
    mixed_max_examples: usize,
    mixed_capture: mixed_capture::Options,
    performance_profile: bool,
    quality_diagnostics: bool,
    quality_max_record_bases: usize,
}

fn parse() -> Result<Option<Options>, String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|s| s == "--help" || s == "-h") {
        print!("{HELP}");
        return Ok(None);
    }
    if args == ["--version"] {
        println!("trm {}", env!("CARGO_PKG_VERSION"));
        return Ok(None);
    }
    if !["scan", "reclassify", "evidence", "mixed"].contains(&args[0].as_str()) {
        return Err("expected 'scan', 'evidence', 'mixed' or 'reclassify'; see --help".into());
    }
    let mut opts = Options {
        action: args[0].clone(),
        prior: prior::PriorOptions::default(),
        support: support::SupportOptions::default(),
        cfg: Config::default(),
        paths: Vec::new(),
        output: String::new(),
        threads: std::thread::available_parallelism().map_or(1, |n| n.get().min(8)),
        max_records: 0,
        batch_bases: 4_000_000,
        memory_limit_mib: None,
        memory_policy: None,
        paired_layout: None,
        targets: None,
        evidence: evidence::Options::default(),
        mixed_definitions: None,
        mixed_config: mixed::MixedConfig::default(),
        mixed_max_examples: 20,
        mixed_capture: mixed_capture::Options::default(),
        performance_profile: false,
        quality_diagnostics: false,
        quality_max_record_bases: 100_000,
    };
    let mut capture_limits_specified = false;
    let mut quality_limit_specified = false;
    let mut i = 1;
    while i < args.len() {
        let key = &args[i];
        if key == "--" {
            opts.paths.extend(args[i + 1..].iter().cloned());
            break;
        }
        if !key.starts_with('-') || key == "-" {
            opts.paths.push(key.clone());
            i += 1;
            continue;
        }
        let value = args
            .get(i + 1)
            .ok_or_else(|| format!("missing value for {key}"))?;
        let num = || {
            value
                .parse::<usize>()
                .map_err(|_| format!("invalid integer for {key}: {value}"))
        };
        if opts.action == "reclassify"
            && ![
                "--output",
                "-o",
                "--prior-policy",
                "--prior-max-edits",
                "--prior-max-distance",
                "--related-motif",
                "--lineage",
                "--background-model",
                "--support-policy",
                "--support-min-enrichment",
            ]
            .contains(&key.as_str())
        {
            return Err(format!("{key} is a scan-only or unknown option; reclassification cannot change read evidence"));
        }
        if opts.action == "evidence"
            && (key.starts_with("--prior-")
                || key.starts_with("--support-")
                || [
                    "--related-motif",
                    "--lineage",
                    "--background-model",
                    "--min-k",
                    "--max-k",
                ]
                .contains(&key.as_str()))
        {
            return Err(format!("{key} does not apply to targeted evidence; targets define the motif families and no prior/background classification is performed"));
        }
        if opts.action == "mixed"
            && ![
                "--mixed-definitions",
                "--mixed-min-span",
                "--mixed-min-component-copies",
                "--mixed-max-examples",
                "--mixed-capture",
                "--mixed-capture-max-records",
                "--mixed-capture-max-bases",
                "--mixed-capture-max-record-bases",
                "--threads",
                "-t",
                "--mode",
                "--output",
                "-o",
                "--max-records",
                "--batch-bases",
                "--memory-limit-mib",
                "--min-tract",
                "--end-window",
                "--junction-flank-bases",
                "--junction-end-slack",
            ]
            .contains(&key.as_str())
        {
            return Err(format!(
                "{key} does not apply to exact mixed-array analysis"
            ));
        }
        if opts.action != "mixed" && key.starts_with("--mixed-") {
            return Err(format!("{key} requires mixed mode"));
        }
        match key.as_str() {
            "--mixed-definitions" => opts.mixed_definitions = Some(value.clone()),
            "--mixed-min-span" => opts.mixed_config.min_span = num()?,
            "--mixed-min-component-copies" => opts.mixed_config.min_component_copies = num()?,
            "--mixed-max-examples" => opts.mixed_max_examples = num()?,
            "--mixed-capture" => {
                opts.mixed_capture.enabled = match value.as_str() {
                    "on" => true,
                    "off" => false,
                    _ => return Err("--mixed-capture must be on or off".into()),
                };
            }
            "--mixed-capture-max-records" => {
                opts.mixed_capture.max_records = num()?;
                capture_limits_specified = true;
            }
            "--mixed-capture-max-bases" => {
                opts.mixed_capture.max_bases = num()?;
                capture_limits_specified = true;
            }
            "--mixed-capture-max-record-bases" => {
                opts.mixed_capture.max_record_bases = num()?;
                capture_limits_specified = true;
            }
            "--targets" => opts.targets = Some(value.clone()),
            "--junction-flank-bases" => opts.evidence.flank_bases = num()?,
            "--junction-end-slack" => opts.evidence.end_slack = num()?,
            "--flank-period-annotation" => {
                opts.evidence.flank_period_annotation = match value.as_str() {
                    "on" => true,
                    "off" => false,
                    _ => return Err("--flank-period-annotation must be on or off".into()),
                };
            }
            "--performance-profile" => {
                opts.performance_profile = match value.as_str() {
                    "on" => true,
                    "off" => false,
                    _ => return Err("--performance-profile must be on or off".into()),
                };
            }
            "--quality-diagnostics" => {
                opts.quality_diagnostics = match value.as_str() {
                    "phred33" => true,
                    "off" => false,
                    _ => return Err("--quality-diagnostics must be off or phred33 (explicit encoding declaration)".into()),
                };
            }
            "--quality-max-record-bases" => {
                quality_limit_specified = true;
                opts.quality_max_record_bases = num()?;
            }
            "--prior-policy" => opts.prior.policy = value.clone(),
            "--prior-max-edits" => opts.prior.max_edits = num()?,
            "--prior-max-distance" => {
                opts.prior.max_distance =
                    value.parse().map_err(|_| "invalid --prior-max-distance")?
            }
            "--related-motif" => opts.prior.related_motifs.push(value.clone()),
            "--lineage" => opts.prior.lineage = Some(value.clone()),
            "--background-model" => opts.support.background_model = Some(value.clone()),
            "--support-policy" => opts.support.policy = value.clone(),
            "--support-min-enrichment" => {
                opts.support.min_enrichment = value
                    .parse()
                    .map_err(|_| "invalid --support-min-enrichment")?
            }
            "--threads" | "-t" => opts.threads = num()?,
            "--min-k" => opts.cfg.min_k = num()?,
            "--max-k" => opts.cfg.max_k = num()?,
            "--seed-bases" => opts.cfg.seed_bases = num()?,
            "--min-copies" => opts.cfg.min_copies = num()?,
            "--min-tract" => opts.cfg.min_tract = num()?,
            "--min-identity" => {
                opts.cfg.min_identity = value.parse().map_err(|_| "invalid --min-identity")?
            }
            "--end-window" => opts.cfg.end_window = num()?,
            "--min-reads" => opts.cfg.min_reads = num()? as u64,
            "--max-candidates" => opts.cfg.max_candidates = num()?,
            "--max-records" => opts.max_records = num()? as u64,
            "--batch-bases" => opts.batch_bases = num()?,
            "--memory-limit-mib" => {
                opts.memory_limit_mib =
                    Some(value.parse().map_err(|_| "invalid --memory-limit-mib")?)
            }
            "--paired-layout" => opts.paired_layout = Some(value.clone()),
            "--mode" => opts.cfg.mode = value.clone(),
            "--output" | "-o" => opts.output = value.clone(),
            _ => return Err(format!("unknown option {key}")),
        }
        i += 2;
    }
    if opts.paths.is_empty() || opts.output.is_empty() {
        return Err("inputs and --output PREFIX are required".into());
    }
    if quality_limit_specified && !opts.quality_diagnostics {
        return Err("--quality-max-record-bases requires --quality-diagnostics phred33".into());
    }
    if !(1..=1_000_000).contains(&opts.quality_max_record_bases) {
        return Err("--quality-max-record-bases must be 1..1000000".into());
    }
    if (opts.action == "mixed") != opts.mixed_definitions.is_some() {
        return Err("--mixed-definitions is required exclusively for mixed mode".into());
    }
    if opts.mixed_max_examples > 10_000 {
        return Err("--mixed-max-examples must be <=10000".into());
    }
    if capture_limits_specified && !opts.mixed_capture.enabled {
        return Err("mixed capture limits require --mixed-capture on".into());
    }
    opts.mixed_capture.validate()?;
    opts.prior.validate()?;
    opts.support.validate()?;
    opts.evidence.validate()?;
    if (opts.action == "evidence") != opts.targets.is_some() {
        return Err("--targets is required exclusively for evidence mode".into());
    }
    if opts.action == "evidence" && opts.support.background_model.is_some() {
        return Err(
            "evidence mode reports read evidence only; background models belong to scan/reclassify"
                .into(),
        );
    }
    if opts.action == "reclassify" && (opts.paths.len() != 1 || opts.paths[0] == "-") {
        return Err("reclassify requires exactly one saved result JSON path".into());
    }
    if opts.paths.iter().filter(|p| p.as_str() == "-").count() > 1 {
        return Err("stdin may only be specified once".into());
    }
    if !(1..=256).contains(&opts.threads) {
        return Err("--threads must be 1..256".into());
    }
    if opts.cfg.min_k < 1 || opts.cfg.min_k > opts.cfg.max_k || opts.cfg.max_k > 64 {
        return Err("require 1 <= min-k <= max-k <= 64".into());
    }
    if opts.cfg.min_copies < 2
        || opts.cfg.min_copies > 1_000_000
        || opts.cfg.seed_bases == 0
        || opts.cfg.seed_bases > 100_000_000
    {
        return Err("min-copies must be 2..1000000 and seed-bases 1..100000000".into());
    }
    if !opts.cfg.min_identity.is_finite() || !(0.5..=1.0).contains(&opts.cfg.min_identity) {
        return Err("--min-identity must be between 0.5 and 1".into());
    }
    if opts.cfg.min_tract == 0
        || opts.cfg.min_reads == 0
        || opts.cfg.max_candidates == 0
        || opts.batch_bases == 0
    {
        return Err("tract/reads/candidates/batch limits must be positive".into());
    }
    if !["reads", "assembly"].contains(&opts.cfg.mode.as_str()) {
        return Err("--mode must be reads or assembly".into());
    }
    if let Some(layout) = &opts.paired_layout {
        if layout != "fr" || !["scan", "evidence"].contains(&opts.action.as_str()) {
            return Err("--paired-layout accepts only fr in scan/evidence mode".into());
        }
        if opts.cfg.mode != "reads" || opts.paths.len() != 2 || opts.paths.iter().any(|p| p == "-")
        {
            return Err("paired FR requires --mode reads and exactly two FASTQ paths (R1 R2), without stdin".into());
        }
        if opts.max_records % 2 != 0 {
            return Err("paired --max-records must be even; a pair is never split".into());
        }
        if fs::canonicalize(&opts.paths[0]).ok().is_some()
            && fs::canonicalize(&opts.paths[0]).ok() == fs::canonicalize(&opts.paths[1]).ok()
        {
            return Err("paired inputs must be two different files".into());
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            if let (Ok(first), Ok(second)) =
                (fs::metadata(&opts.paths[0]), fs::metadata(&opts.paths[1]))
            {
                if first.dev() == second.dev() && first.ino() == second.ino() {
                    return Err(
                        "paired inputs must be two different files, not hard-link aliases".into(),
                    );
                }
            }
        }
    }
    Ok(Some(opts))
}

#[derive(Serialize, Clone)]
struct Example {
    input_index: usize,
    record: String,
    start: usize,
    end: usize,
    strand: String,
    identity: f64,
}

#[derive(Default, Clone, Serialize)]
struct Candidate {
    motif: String,
    reverse_complement: String,
    period: usize,
    distinct_bases: usize,
    reads: u64,
    repeat_bases: u64,
    max_tract: usize,
    long_tract_reads: u64,
    forward_reads: u64,
    reverse_reads: u64,
    ambiguous_strand_reads: u64,
    left_end_reads: u64,
    right_end_reads: u64,
    left_forward_reads: u64,
    left_reverse_reads: u64,
    right_forward_reads: u64,
    right_reverse_reads: u64,
    long_left_end_reads: u64,
    long_right_end_reads: u64,
    opposed_end_records: u64,
    full_repeat_reads: u64,
    status: String,
    rank_score: f64,
    examples: Vec<Example>,
    junction_evidence: evidence::JunctionEvidence,
    #[serde(skip_serializing_if = "Option::is_none")]
    paired_evidence: Option<paired::PairStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    quality_diagnostics: Option<Box<quality::CandidateQuality>>,
}

fn record_hit(
    c: &mut Candidate,
    motif: &[u8],
    tracts: &[Tract],
    record: &input::Record,
    input_index: usize,
    cfg: &Config,
    evidence_options: &evidence::Options,
) {
    c.reads += 1;
    let n = record.seq.len();
    c.junction_evidence
        .observe_with_motif(&record.seq, tracts, evidence_options, motif);
    let has = |predicate: &dyn Fn(&Tract) -> bool| tracts.iter().any(predicate) as u64;
    c.forward_reads += has(&|t| t.strand == Some(true));
    c.reverse_reads += has(&|t| t.strand == Some(false));
    c.ambiguous_strand_reads += has(&|t| t.strand.is_none());
    let left = |t: &Tract| t.start <= cfg.end_window;
    let right = |t: &Tract| n - t.end <= cfg.end_window;
    c.left_end_reads += has(&left);
    c.right_end_reads += has(&right);
    c.left_forward_reads += has(&|t| left(t) && t.strand == Some(true));
    c.left_reverse_reads += has(&|t| left(t) && t.strand == Some(false));
    c.right_forward_reads += has(&|t| right(t) && t.strand == Some(true));
    c.right_reverse_reads += has(&|t| right(t) && t.strand == Some(false));
    c.long_left_end_reads += has(&|t| left(t) && t.end - t.start >= cfg.min_tract);
    c.long_right_end_reads += has(&|t| right(t) && t.end - t.start >= cfg.min_tract);
    c.long_tract_reads += has(&|t| t.end - t.start >= cfg.min_tract);
    c.opposed_end_records += tracts.iter().any(|a| {
        left(a)
            && a.strand.is_some()
            && tracts
                .iter()
                .any(|b| right(b) && a.end <= b.start && b.strand.is_some() && a.strand != b.strand)
    }) as u64;
    c.full_repeat_reads += has(&|t| (t.end - t.start) as f64 >= n as f64 * 0.9);
    let mut covered_end = 0;
    for t in tracts {
        c.repeat_bases += t.end.saturating_sub(covered_end.max(t.start)) as u64;
        covered_end = covered_end.max(t.end);
        c.max_tract = c.max_tract.max(t.end - t.start);
        if c.examples.len() < 3 {
            c.examples.push(Example {
                input_index,
                record: record.id.clone(),
                start: t.start,
                end: t.end,
                strand: match t.strand {
                    Some(true) => "+",
                    Some(false) => "-",
                    None => ".",
                }
                .into(),
                identity: t.identity,
            });
        }
    }
}

fn merge(a: &mut Candidate, b: Candidate) {
    if let Some(incoming) = b.quality_diagnostics {
        a.quality_diagnostics
            .get_or_insert_with(Default::default)
            .merge(*incoming);
    }
    if let Some(pairs) = b.paired_evidence {
        a.paired_evidence
            .get_or_insert_with(Default::default)
            .merge(pairs);
    }
    a.junction_evidence.merge(b.junction_evidence);
    a.reads += b.reads;
    a.repeat_bases += b.repeat_bases;
    a.max_tract = a.max_tract.max(b.max_tract);
    a.long_tract_reads += b.long_tract_reads;
    a.forward_reads += b.forward_reads;
    a.reverse_reads += b.reverse_reads;
    a.ambiguous_strand_reads += b.ambiguous_strand_reads;
    a.left_end_reads += b.left_end_reads;
    a.right_end_reads += b.right_end_reads;
    a.left_forward_reads += b.left_forward_reads;
    a.left_reverse_reads += b.left_reverse_reads;
    a.right_forward_reads += b.right_forward_reads;
    a.right_reverse_reads += b.right_reverse_reads;
    a.long_left_end_reads += b.long_left_end_reads;
    a.long_right_end_reads += b.long_right_end_reads;
    a.opposed_end_records += b.opposed_end_records;
    a.full_repeat_reads += b.full_repeat_reads;
    a.examples.extend(b.examples);
    a.examples.sort_by(|x, y| {
        (x.input_index, &x.record, x.start).cmp(&(y.input_index, &y.record, y.start))
    });
    a.examples.truncate(3);
}

enum ScanUnit {
    Single(usize, input::Record, bool),
    Pair(input::Record, input::Record),
}
type ScanBatch = Vec<ScanUnit>;
type Counts = HashMap<String, Candidate>;

struct WorkerResult {
    counts: Counts,
    qc: evidence::QualitySummary,
    quality: Option<quality::RecordQuality>,
}

fn worker(
    rx: Arc<Mutex<mpsc::Receiver<ScanBatch>>>,
    cfg: Config,
    cancelled: Arc<AtomicBool>,
    evidence_options: evidence::Options,
    targets: Option<Arc<detect::TargetIndex>>,
    quality_diagnostics: bool,
) -> Result<WorkerResult, String> {
    let mut counts: Counts = HashMap::new();
    let mut qc = evidence::QualitySummary::default();
    let mut quality = quality_diagnostics.then(quality::RecordQuality::default);
    loop {
        if cancelled.load(Ordering::Relaxed) {
            break;
        }
        let batch = rx.lock().map_err(|_| "worker queue lock failed")?.recv();
        let Ok(batch) = batch else {
            break;
        };
        for unit in batch {
            if cancelled.load(Ordering::Relaxed) {
                return Ok(WorkerResult {
                    counts,
                    qc,
                    quality,
                });
            }
            let mut observe = |input_index, record: &input::Record, is_fastq: bool| {
                qc.observe(&record.seq);
                if let Some(summary) = &mut quality {
                    summary.observe(record.quality.as_deref(), is_fastq);
                }
                let discovered = match &targets {
                    Some(index) => index.discover(&record.seq, &cfg),
                    None => detect::discover(&record.seq, &cfg),
                };
                for (key, tracts) in &discovered {
                    if !counts.contains_key(key) && counts.len() >= cfg.max_candidates {
                        cancelled.store(true, Ordering::Relaxed);
                        return Err(
                            "candidate limit exceeded; increase --max-candidates or --seed-bases"
                                .to_string(),
                        );
                    }
                    record_hit(
                        counts.entry(key.clone()).or_default(),
                        key.as_bytes(),
                        tracts,
                        record,
                        input_index,
                        &cfg,
                        &evidence_options,
                    );
                    if quality_diagnostics {
                        counts
                            .get_mut(key)
                            .unwrap()
                            .quality_diagnostics
                            .get_or_insert_with(Default::default)
                            .observe(record.quality.as_deref(), is_fastq, tracts);
                    }
                }
                Ok(discovered)
            };
            match unit {
                ScanUnit::Single(input_index, record, is_fastq) => {
                    observe(input_index, &record, is_fastq)?;
                }
                ScanUnit::Pair(r1, r2) => {
                    let first = observe(0, &r1, true)?;
                    let second = observe(1, &r2, true)?;
                    for key in first
                        .keys()
                        .chain(second.keys().filter(|key| !first.contains_key(*key)))
                    {
                        let left = first.get(key).map_or(&[][..], Vec::as_slice);
                        let right = second.get(key).map_or(&[][..], Vec::as_slice);
                        counts
                            .get_mut(key)
                            .unwrap()
                            .paired_evidence
                            .get_or_insert_with(Default::default)
                            .observe(left, r1.seq.len(), right, r2.seq.len());
                    }
                }
            }
        }
    }
    Ok(WorkerResult {
        counts,
        qc,
        quality,
    })
}

#[derive(Serialize)]
struct Source {
    path: String,
    bytes: Option<u64>,
    modified_unix_seconds: Option<u64>,
}
#[derive(Default, Serialize)]
struct Summary {
    records: u64,
    bases: u64,
    max_record_bases: usize,
    candidate_families_before_filter: usize,
    reported_candidates: usize,
    sampled_prefix: bool,
    elapsed_seconds: f64,
    megabases_per_second: f64,
    call: String,
}
#[derive(Serialize)]
struct Report {
    schema_version: u32,
    program: &'static str,
    version: &'static str,
    command: Vec<String>,
    inputs: Vec<Source>,
    config: Config,
    threads: usize,
    summary: Summary,
    warnings: Vec<String>,
    candidates: Vec<Candidate>,
}

fn scan_reader(
    opts: &Options,
    path: &str,
    paired: bool,
    profile: &mut Option<profile::Producer>,
) -> Result<input::Reader, String> {
    let source = if let Some(profile) = profile {
        let metrics = profile::input_metrics(path);
        let source = input::open_profiled(path, metrics.clone())?;
        profile.inputs.push(metrics);
        source
    } else {
        input::open(path)?
    };
    let mut reader = if paired {
        input::Reader::paired(source)
    } else {
        input::Reader::new(source)
    };
    if opts.quality_diagnostics {
        reader.retain_quality(opts.quality_max_record_bases);
    }
    Ok(reader)
}

fn send_scan_batch(
    tx: &mpsc::SyncSender<ScanBatch>,
    batch: ScanBatch,
    profile: &mut Option<profile::Producer>,
) -> Result<(), String> {
    if let Some(profile) = profile {
        let start = Instant::now();
        tx.send(batch).map_err(|_| "all workers stopped")?;
        profile.queue_send_wall_seconds += start.elapsed().as_secs_f64();
        profile.batches_sent += 1;
    } else {
        tx.send(batch).map_err(|_| "all workers stopped")?;
    }
    Ok(())
}

fn read_pairs(
    opts: &Options,
    tx: &mpsc::SyncSender<ScanBatch>,
    cancelled: &AtomicBool,
    summary: &mut Summary,
    profile: &mut Option<profile::Producer>,
) -> Result<(), String> {
    let mut first = scan_reader(opts, &opts.paths[0], true, profile)?;
    let mut second = scan_reader(opts, &opts.paths[1], true, profile)?;
    let mut batch = Vec::new();
    let mut batch_bases = 0;
    loop {
        if cancelled.load(Ordering::Relaxed) {
            return Err("paired analysis cancelled after worker failure".into());
        }
        if opts.max_records > 0 && summary.records >= opts.max_records {
            summary.sampled_prefix = true;
            break;
        }
        let left = first
            .next_record()
            .map_err(|e| format!("{}: {e}", opts.paths[0]))?;
        let right = second
            .next_record()
            .map_err(|e| format!("{}: {e}", opts.paths[1]))?;
        let (r1, r2) = match (left, right) {
            (None, None) => break,
            (Some(r1), Some(r2)) => (r1, r2),
            _ => return Err(format!("paired inputs have unequal record counts after {} complete pairs; no rematching or orphan dropping is performed", summary.records / 2)),
        };
        if !first.is_fastq() || !second.is_fastq() {
            return Err("--paired-layout fr requires two FASTQ inputs, not FASTA".into());
        }
        paired::validate_pair_ids(first.pair_header(), second.pair_header())
            .map_err(|e| format!("pair {}: {e}", summary.records / 2 + 1))?;
        summary.records += 2;
        let bases = r1.seq.len() + r2.seq.len();
        summary.bases += bases as u64;
        summary.max_record_bases = summary.max_record_bases.max(r1.seq.len()).max(r2.seq.len());
        batch_bases += bases;
        batch.push(ScanUnit::Pair(r1, r2));
        if batch_bases >= opts.batch_bases {
            send_scan_batch(tx, std::mem::take(&mut batch), profile)?;
            batch_bases = 0;
        }
    }
    if !batch.is_empty() {
        send_scan_batch(tx, batch, profile)?;
    }
    Ok(())
}

fn run(mut opts: Options) -> Result<(), String> {
    protect_inputs(&opts)?;
    if opts.action == "reclassify" {
        return reclassify(opts);
    }
    let policy = if opts.quality_diagnostics {
        memory::prepare_with_payload(opts.memory_limit_mib, opts.threads, opts.batch_bases, 1)?
    } else {
        memory::prepare(opts.memory_limit_mib, opts.threads, opts.batch_bases)?
    };
    opts.batch_bases = policy.effective_batch_bases;
    if opts.memory_limit_mib.is_some() {
        eprintln!("Memory limit: Linux virtual address space, not RSS; effective batch {} bases. Resource exhaustion is a failed analysis, never a biological negative.", opts.batch_bases);
    }
    opts.memory_policy = Some(policy);
    if opts.action == "mixed" {
        return mixed_scan::run(opts);
    }
    let targets = opts
        .targets
        .as_ref()
        .map(|path| {
            evidence::load_targets(path)
                .map(|families| Arc::new(detect::TargetIndex::new(families, &opts.cfg)))
        })
        .transpose()?;
    support::validate_config(
        &serde_json::to_value(&opts.cfg).map_err(|e| e.to_string())?,
        &opts.support,
    )?;
    let start = Instant::now();
    let mut sources = Vec::new();
    for p in &opts.paths {
        let md = if p == "-" {
            None
        } else {
            Some(fs::metadata(p).map_err(|e| format!("{p}: {e}"))?)
        };
        sources.push(Source {
            path: p.clone(),
            bytes: md.as_ref().map(|m| m.len()),
            modified_unix_seconds: md
                .and_then(|m| m.modified().ok())
                .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
                .map(|d| d.as_secs()),
        });
    }
    let (tx, rx) = mpsc::sync_channel::<ScanBatch>(opts.threads * 2);
    let rx = Arc::new(Mutex::new(rx));
    let cancelled = Arc::new(AtomicBool::new(false));
    let mut handles = Vec::new();
    for _ in 0..opts.threads {
        let receiver = Arc::clone(&rx);
        let cfg = opts.cfg.clone();
        let stop = Arc::clone(&cancelled);
        let evidence_options = opts.evidence.clone();
        let target_index = targets.clone();
        let quality_diagnostics = opts.quality_diagnostics;
        handles.push(
            std::thread::Builder::new()
                .spawn(move || {
                    worker(
                        receiver,
                        cfg,
                        stop,
                        evidence_options,
                        target_index,
                        quality_diagnostics,
                    )
                })
                .map_err(|e| format!("cannot start worker (check memory/thread limits): {e}"))?,
        );
    }
    drop(rx);
    let mut summary = Summary::default();
    let paired_enabled = opts.paired_layout.is_some();
    let mut producer_profile = opts.performance_profile.then(profile::Producer::new);
    let parse_result: Result<(), String> = if paired_enabled {
        read_pairs(&opts, &tx, &cancelled, &mut summary, &mut producer_profile)
    } else {
        (|| {
            let mut batch = Vec::new();
            let mut batch_bases = 0;
            'files: for (input_index, path) in opts.paths.iter().enumerate() {
                let mut reader = scan_reader(&opts, path, false, &mut producer_profile)?;
                loop {
                    if cancelled.load(Ordering::Relaxed) {
                        return Err("analysis cancelled after worker failure".into());
                    }
                    if opts.max_records > 0 && summary.records >= opts.max_records {
                        summary.sampled_prefix = true;
                        break 'files;
                    }
                    let Some(record) = reader.next_record().map_err(|e| format!("{path}: {e}"))?
                    else {
                        break;
                    };
                    summary.records += 1;
                    summary.bases += record.seq.len() as u64;
                    summary.max_record_bases = summary.max_record_bases.max(record.seq.len());
                    batch_bases += record.seq.len();
                    batch.push(ScanUnit::Single(input_index, record, reader.is_fastq()));
                    if batch_bases >= opts.batch_bases {
                        send_scan_batch(&tx, std::mem::take(&mut batch), &mut producer_profile)?;
                        batch_bases = 0;
                    }
                }
            }
            if !batch.is_empty() {
                send_scan_batch(&tx, batch, &mut producer_profile)?;
            }
            Ok(())
        })()
    };
    if let Some(profile) = &mut producer_profile {
        profile.finish();
    }
    if parse_result.is_err() {
        cancelled.store(true, Ordering::Relaxed);
    }
    drop(tx);
    let mut counts: Counts = HashMap::new();
    let mut worker_error = None;
    let mut qc = evidence::QualitySummary::default();
    let mut quality = opts
        .quality_diagnostics
        .then(quality::RecordQuality::default);
    for h in handles {
        match h
            .join()
            .map_err(|_| "worker panicked".to_string())
            .and_then(|r| r)
        {
            Ok(part) => {
                qc.merge(part.qc);
                if let (Some(summary), Some(incoming)) = (&mut quality, part.quality) {
                    summary.merge(incoming);
                }
                for (key, c) in part.counts {
                    if !counts.contains_key(&key) && counts.len() >= opts.cfg.max_candidates {
                        worker_error = Some("global candidate limit exceeded; increase --max-candidates or --seed-bases".into());
                        break;
                    }
                    merge(counts.entry(key).or_default(), c);
                }
            }
            Err(e) => worker_error = Some(e),
        }
    }
    if let Some(e) = worker_error {
        return Err(e);
    }
    parse_result?;
    if summary.records == 0 {
        return Err("no sequence records found".into());
    }
    summary.candidate_families_before_filter = counts.len();
    if let Some(index) = &targets {
        for family in &index.families {
            counts.entry(family.clone()).or_default();
        }
    }
    let mut candidates = Vec::new();
    for (key, mut c) in counts {
        if targets.is_none() && c.reads < opts.cfg.min_reads {
            continue;
        }
        if paired_enabled && c.paired_evidence.is_none() {
            c.paired_evidence = Some(paired::PairStats::default());
        }
        if opts.quality_diagnostics && c.quality_diagnostics.is_none() {
            c.quality_diagnostics = Some(Box::default());
        }
        if opts.evidence.flank_period_annotation {
            c.junction_evidence
                .init_flank_period_annotation(key.as_bytes());
        }
        c.period = key.len();
        c.distinct_bases = b"ACGT"
            .iter()
            .filter(|b| key.as_bytes().contains(b))
            .count();
        c.reverse_complement =
            String::from_utf8(motif::reverse_complement(key.as_bytes())).unwrap();
        c.motif = key;
        c.status =
            if opts.cfg.mode == "assembly" && c.long_left_end_reads + c.long_right_end_reads > 0 {
                "terminal_repeat_candidate"
            } else if c.long_tract_reads > 0 {
                "long_tract_supported"
            } else {
                "repeat_candidate"
            }
            .into();
        let gc_skew = (c.motif.bytes().filter(|&b| b == b'G').count() as f64
            - c.motif.bytes().filter(|&b| b == b'C').count() as f64)
            .abs()
            / c.period as f64;
        // A transparent prioritization heuristic, not a calibrated probability.
        c.rank_score = c.repeat_bases as f64
            * (1.0 + c.max_tract as f64).ln()
            * (1.0 + gc_skew)
            * (1.0 + (c.long_tract_reads as f64).ln_1p());
        if opts.cfg.mode == "assembly" {
            c.rank_score *= 1.0 + (c.long_left_end_reads + c.long_right_end_reads) as f64;
        }
        candidates.push(c);
    }
    candidates.sort_by(|a, b| {
        let terminal_priority = |c: &Candidate| {
            if opts.cfg.mode != "assembly" {
                (false, false, false)
            } else {
                (
                    c.long_left_end_reads + c.long_right_end_reads > 0,
                    c.left_end_reads + c.right_end_reads > 0,
                    c.opposed_end_records > 0,
                )
            }
        };
        terminal_priority(b)
            .cmp(&terminal_priority(a))
            .then_with(|| {
                b.rank_score
                    .total_cmp(&a.rank_score)
                    .then(a.motif.cmp(&b.motif))
            })
    });
    summary.reported_candidates = candidates.len();
    summary.call = if candidates.is_empty() {
        "no_call"
    } else {
        "candidates_only"
    }
    .into();
    summary.elapsed_seconds = start.elapsed().as_secs_f64();
    summary.megabases_per_second = summary.bases as f64 / 1e6 / summary.elapsed_seconds;
    let mut warnings = vec![
        "Tandem-repeat candidates are not proof of telomeric location; satellites and internal repeats can be reported.".into(),
        "No call does not establish motif absence, ALT, or telomere maintenance mechanism.".into(),
        "Strands refer to the canonical motif in the input sequence, not an inferred chromosome 5-prime/3-prime end.".into(),
        "Reads are counted as records; mates and PCR duplicates are not independent molecules or deduplicated.".into(),
        "Extension tolerates substitutions only; indels, irregular telomeres, or insufficient exact seeds can cause missed/fragmented tracts.".into(),
    ];
    if opts.cfg.mode == "reads" {
        warnings.push(
            "Read-end proximity is descriptive and does not establish a chromosome end.".into(),
        );
    }
    if summary.max_record_bases < opts.cfg.min_tract {
        warnings.push("All reads are shorter than min-tract; contiguous tract support cannot be established from these reads.".into());
    }
    if summary.sampled_prefix {
        warnings.push("Only an input prefix was analyzed; this is not random sampling and may miss rare telomeric reads.".into());
    }
    let report = Report {
        schema_version: 1,
        program: "trm",
        version: env!("CARGO_PKG_VERSION"),
        command: std::env::args().collect(),
        inputs: sources,
        config: opts.cfg,
        threads: opts.threads,
        summary,
        warnings,
        candidates,
    };
    let mut value = serde_json::to_value(&report).map_err(|e| e.to_string())?;
    value["memory_policy"] =
        serde_json::to_value(&opts.memory_policy).map_err(|e| e.to_string())?;
    value["paired_end"] = serde_json::json!({
        "measured": paired_enabled,
        "layout": opts.paired_layout,
        "validated_pairs": if paired_enabled { Some(report.summary.records / 2) } else { None },
        "scope": if paired_enabled && report.summary.sampled_prefix { "synchronized_prefix_only" } else if paired_enabled { "complete_inputs" } else { "not_measured" },
        "method": "Optional synchronized FR mates; exact normalized header IDs are checked without rematching. Opposite canonical-family strands among both-repeat mates are consistency evidence. One-mate support is neutral; ambiguous strands are separate. R1/R2 are not chromosome strands. Counts remain records, not deduplicated molecules; no alignment, insert-size, gap, overlap or telomere-location inference. No change to classification or adequacy."
    });
    value["evidence_qc"] = qc.json();
    if let Some(profile) = producer_profile {
        value["performance_profile"] =
            profile.json(report.summary.records, report.summary.bases)?;
    }
    if let Some(quality) = quality {
        value["evidence_qc"]["quality_scores_assessed"] = (quality.quality_records > 0).into();
        value["quality_diagnostics"] = serde_json::json!({
            "schema_version": 1,
            "encoding": "phred33",
            "encoding_provenance": "explicit_user_declaration_not_inferred",
            "max_record_bases": opts.quality_max_record_bases,
            "memory_plan_extra_bytes_per_sequence_base": 1,
            "input_scope": if report.summary.sampled_prefix { "analyzed_prefix" } else { "complete_inputs" },
            "input_records": quality,
            "filtering_applied": false,
            "variant_discriminating_bases_assessed": false,
            "method": "Descriptive Phred+33 histograms; user declares encoding. All complete retained FASTQ quality is summarized globally, and supporting whole reads plus the union of each candidate's detected tract coordinates are summarized per family. Scores do not change repeat detection, ranking, support, or adequacy. FASTA has no quality; longer FASTQ records are entirely omitted from quality retention, not assigned synthetic scores. All records still undergo normal complete format validation. Candidate histograms overlap between families, records may be duplicate molecules, and these are not phase-aware variant-discriminating-base scores.",
            "memory": "Opt-in quality retention adds up to one quality byte per sequence base for eligible queued/worker records, plus per-candidate histograms. The record-length cap is not a total memory budget; batch queues, worker count, parser, candidate tables, allocator and output also consume memory. The existing Linux address-space cap still applies; exhaustion is a failed analysis."
        });
    }
    value["junction_options"] =
        serde_json::to_value(&opts.evidence).map_err(|err| err.to_string())?;
    value["junction_evidence_method"] = serde_json::json!({
        "orientation": "Forward/reverse refer to the canonical motif, never an inferred G-rich strand. Map raw left/right quadrants to a declared reference orientation downstream.",
        "flank": "At least flank_bases contiguous A/C/G/T, >=3 distinct bases, and no >=85% lag agreement at periods 1..min(32,flank_bases/3); no overlapping same-family tract.",
        "full_repeat": "Any tract covers >=90% of the record; such records receive no junction credit.",
        "interpretation": "Read-boundary repeat plus flank geometry is not proof of a physical chromosome end; fragmentation, library bias and interstitial repeats can mimic it.",
        "validation": {
            "schema_version": 1,
            "read_signatures": "Bottom-128 FNV-1a/128 hashes of min(read, IUPAC reverse complement), hexadecimal numeric order. Only junction-bearing records are profiled. Read-hash collisions remain possible; these are sequence signatures, not independent molecules.",
            "flank_contexts": "Nearest clean 20bp, A=0/C=1/G=2/T=3 packed into a collision-free 40-bit identity. Left junction uses its right flank; right junction uses the reverse complement of its left flank, physically orienting both repeat-left. Retain bottom-128 by FNV-1a/128 rank then exact packed identity.",
            "representative_flank": "At most 80 adjacent A/C/G/T bases, stopped by ambiguity or another same-family tract, normalized repeat-left; longest then lexicographically smallest representative. Bases beyond the identity20 are sequence context, not a claim of nonrepeat sequence or mapping uniqueness.",
            "quadrants": "quadrant_counts preserve input-record orientation. normalized_quadrant_counts reverse L/F with R/R and L/R with R/F when normalizing a read, or when physically reversing a right flank. Ambiguous counts collapse left/right ambiguity per record for read signatures, per event for flank contexts.",
            "counts": "Retained entry counts are exact for their signatures across merges. observed_records/events include discarded keys. distinct_lower_bound is retained count plus one when saturated, not a diversity estimate. Saturated omitted distinct keys are at most observed count minus sum of retained entry counts. No distinct-context count establishes loci, chromosome ends, karyotype or molecule independence.",
            "clean_two_flank_repeat_records": "Read contains a tract bounded on both sides by the configured clean flank window (>=20bp), with no overlapping same-family tract in either flank. Count once per record, including >=90% repeat-covered records if they have both clean flanks. Existing internal_repeat_reads and all v0.4 count rules are unchanged.",
        },
    });
    if opts.evidence.flank_period_annotation {
        value["flank_period_annotation_method"] = serde_json::json!({
            "schema_version": 1,
            "scope": "Every previously eligible junction tract-side event, not a retained sketch or a count of unique reads/molecules. Existing counters, ranking, filtering and TRM calls are unchanged.",
            "window": "max(20,3*primitive candidate period), at most 192 bp; adjacent A/C/G/T only, stopped by read boundary, ambiguous base or another detected same-family tract. A shorter window is insufficient, not negative.",
            "lag": "At least 85% agreement between bases k positions apart within the complete fixed window.",
            "family": "At least 85% ungapped identity to the best periodic template over all phases and both strands of the candidate family. Phase-free compatibility does not establish continuous repeat extension.",
            "classes": "insufficient_window / neither / lag_only / family_only / both partition observed events. neither is not proof of nonrepetitive sequence, uniqueness or subtelomere; signals are descriptive, not calibrated probabilities.",
            "orientation": "Each class preserves raw left/right and canonical forward/reverse/ambiguous quadrants; no inferred G-rich direction.",
            "absence": "Missing annotation is unmeasured. An enabled zero-hit target has measured zero events, not evidence of telomere absence. measured_all_events=false marks mixing with unannotated junction-bearing evidence."
        });
    }
    if let Some(index) = targets {
        value["analysis_kind"] = "targeted_evidence".into();
        value["schema_version"] = 4.into();
        value["targets_source"] = opts.targets.clone().into();
        value["target_motifs"] = serde_json::json!(index.families);
        value["summary"]["call"] = "targeted_evidence_only".into();
        value["candidate_evidence"] = value.as_object_mut().unwrap().remove("candidates").unwrap();
        value["rejected_candidates"] = serde_json::json!([]);
        let elapsed = start.elapsed().as_secs_f64();
        value["summary"]["elapsed_seconds"] = elapsed.into();
        value["summary"]["megabases_per_second"] =
            (report.summary.bases as f64 / 1e6 / elapsed).into();
        write_reports(&opts.output, &value)?;
        eprintln!(
            "{} records, {:.2} Mb, {} targeted families; {:.3}s ({:.2} Mb/s)",
            report.summary.records,
            report.summary.bases as f64 / 1e6,
            index.families.len(),
            elapsed,
            report.summary.bases as f64 / 1e6 / elapsed
        );
        return Ok(());
    }
    value["analysis_kind"] = "de_novo_discovery".into();
    prior::apply(&mut value, &opts.prior)?;
    support::apply(&mut value, &opts.support)?;
    let elapsed = start.elapsed().as_secs_f64();
    value["summary"]["elapsed_seconds"] = elapsed.into();
    value["summary"]["megabases_per_second"] = (report.summary.bases as f64 / 1e6 / elapsed).into();
    write_reports(&opts.output, &value)?;
    eprintln!(
        "{} records, {:.2} Mb, {} eligible / {} repeat candidates; {:.3}s ({:.2} Mb/s)",
        report.summary.records,
        report.summary.bases as f64 / 1e6,
        value["summary"]["eligible_telomere_candidates"],
        report.candidates.len(),
        elapsed,
        report.summary.bases as f64 / 1e6 / elapsed
    );
    Ok(())
}

fn protect_inputs(opts: &Options) -> Result<(), String> {
    for ext in ["json", "tsv", "rejected.tsv"] {
        let canonical_output = fs::canonicalize(format!("{}.{}", opts.output, ext)).ok();
        if canonical_output.is_some()
            && opts
                .paths
                .iter()
                .chain(opts.support.background_model.iter())
                .chain(opts.targets.iter())
                .chain(opts.mixed_definitions.iter())
                .any(|p| fs::canonicalize(p).ok() == canonical_output)
        {
            return Err("output would overwrite an input file".into());
        }
    }
    Ok(())
}

fn reclassify(opts: Options) -> Result<(), String> {
    let start = Instant::now();
    let source = &opts.paths[0];
    let reader = BufReader::new(File::open(source).map_err(|e| format!("{source}: {e}"))?);
    let mut report: serde_json::Value =
        serde_json::from_reader(reader).map_err(|e| format!("{source}: {e}"))?;
    let source_version = report.get("version").cloned();
    prior::apply(&mut report, &opts.prior)?;
    support::apply(&mut report, &opts.support)?;
    report["reclassification"] = serde_json::json!({
        "program_version": env!("CARGO_PKG_VERSION"),
        "source_result": fs::canonicalize(source).map_err(|e| e.to_string())?,
        "source_scan_version": source_version,
        "command": std::env::args().collect::<Vec<_>>(),
        "elapsed_seconds": start.elapsed().as_secs_f64(),
        "fastq_reread": false,
        "scan_metrics_preserved": true,
    });
    write_reports(&opts.output, &report)?;
    eprintln!(
        "Reclassified {} repeat candidates: {} eligible; {:.3}s (saved JSON only)",
        report["summary"]["raw_reported_candidates"],
        report["summary"]["eligible_telomere_candidates"],
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

fn write_candidate_tsv(path: &str, candidates: &serde_json::Value) -> std::io::Result<()> {
    let fields = [
        "motif",
        "reverse_complement",
        "period",
        "reads",
        "repeat_bases",
        "max_tract",
        "long_tract_reads",
        "forward_reads",
        "reverse_reads",
        "left_end_reads",
        "right_end_reads",
        "status",
        "rank_score",
        "candidate_eligible",
    ];
    let prior_fields = [
        "class",
        "eligible",
        "confidence_tier",
        "reference_evidence",
        "matched_reference",
        "edit_distance",
        "normalized_distance",
        "source_lineage",
        "explanation",
    ];
    let support_fields = [
        "status",
        "eligible",
        "reads_per_million",
        "repeat_bases_per_million",
        "full_repeat_fraction",
        "read_enrichment_lower",
        "repeat_base_enrichment",
        "reasons",
    ];
    let junction_fields = [
        "full_repeat_records",
        "junction_reads",
        "left_forward_junction_reads",
        "left_reverse_junction_reads",
        "left_ambiguous_junction_reads",
        "right_forward_junction_reads",
        "right_reverse_junction_reads",
        "right_ambiguous_junction_reads",
        "internal_repeat_reads",
    ];
    let pair_fields = [
        "supporting_pairs",
        "both_mates_repeat_pairs",
        "r1_only_repeat_pairs",
        "r2_only_repeat_pairs",
        "opposite_strand_pairs",
        "same_strand_pairs",
        "ambiguous_strand_pairs",
        "both_full_repeat_pairs",
        "exactly_one_full_repeat_pairs",
    ];
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}",
        fields.join("\t"),
        prior_fields
            .iter()
            .map(|name| format!("prior_{name}"))
            .collect::<Vec<_>>()
            .join("\t"),
        support_fields
            .iter()
            .map(|name| format!("support_{name}"))
            .collect::<Vec<_>>()
            .join("\t"),
        junction_fields
            .iter()
            .map(|name| format!("junction_{name}"))
            .collect::<Vec<_>>()
            .join("\t"),
        pair_fields
            .iter()
            .map(|name| format!("pair_{name}"))
            .collect::<Vec<_>>()
            .join("\t")
    )?;
    if let Some(rows) = candidates.as_array() {
        for row in rows {
            let values = fields
                .iter()
                .map(|name| &row[name])
                .chain(prior_fields.iter().map(|name| &row["telomere_prior"][name]))
                .chain(support_fields.iter().map(|name| &row["support"][name]))
                .chain(
                    junction_fields
                        .iter()
                        .map(|name| &row["junction_evidence"][name]),
                )
                .chain(pair_fields.iter().map(|name| &row["paired_evidence"][name]))
                .map(|v| match v {
                    serde_json::Value::Null => String::new(),
                    serde_json::Value::String(s) => s.replace(['\t', '\n', '\r'], " "),
                    _ => v.to_string(),
                })
                .collect::<Vec<_>>();
            writeln!(out, "{}", values.join("\t"))?;
        }
    }
    out.flush()
}

fn write_reports(prefix: &str, report: &serde_json::Value) -> Result<(), String> {
    let path = Path::new(prefix);
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let suffix = format!(".{}.part", std::process::id());
    let json_path = format!("{prefix}.json");
    let tsv_path = format!("{prefix}.tsv");
    let rejected_path = format!("{prefix}.rejected.tsv");
    let json_tmp = format!("{json_path}{suffix}");
    let tsv_tmp = format!("{tsv_path}{suffix}");
    let rejected_tmp = format!("{rejected_path}{suffix}");
    let result = (|| -> std::io::Result<()> {
        let mut json = BufWriter::new(File::create(&json_tmp)?);
        serde_json::to_writer_pretty(&mut json, report)?;
        writeln!(json)?;
        json.flush()?;
        let rows = if report["analysis_kind"] == "targeted_evidence" {
            &report["candidate_evidence"]
        } else {
            &report["candidates"]
        };
        write_candidate_tsv(&tsv_tmp, rows)?;
        write_candidate_tsv(&rejected_tmp, &report["rejected_candidates"])?;
        fs::rename(&tsv_tmp, &tsv_path)?;
        fs::rename(&rejected_tmp, &rejected_path)?;
        fs::rename(&json_tmp, &json_path)?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(json_tmp);
        let _ = fs::remove_file(tsv_tmp);
        let _ = fs::remove_file(rejected_tmp);
    }
    result.map_err(|e| format!("writing output: {e}"))
}

fn main() {
    if let Err(e) = parse().and_then(|opts| opts.map_or(Ok(()), run)) {
        eprintln!("trm: {e}");
        std::process::exit(1);
    }
}
