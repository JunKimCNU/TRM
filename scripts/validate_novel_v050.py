#!/usr/bin/env python3
"""Predeclared synthetic evidence challenges, not biological FDR calibration.

Prepare fixtures before inspecting v0.5 real-data results; then compare frozen
v0.4/v0.5 full discovery and assessment on exactly the same files.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SEED = 20260923
UNKNOWN = 'ACGATCG'


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def canonical(s):
    return min(t[i:]+t[:i] for t in (s,rc(s)) for i in range(len(t)))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def write_tsv(path, rows):
    with path.open('w') as target:
        writer = csv.DictWriter(target, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def dna(rng, n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def repeat(motif, n, phase=0):
    return (motif*((n+phase)//len(motif)+2))[phase:phase+n]


def terminal(rng, motif=UNKNOWN, ends=4, windows=8, full=80):
    """Unique subtelomere blocks + >=600-bp arrays, deterministic150-bp windows."""
    reads = [repeat(motif,150,i%len(motif)) for i in range(full)]
    reads = [r if i%2 else rc(r) for i,r in enumerate(reads)]
    for _ in range(ends):
        subtelomere = dna(rng,200)
        for index in range(windows):
            tract = 45 + index*7
            read = subtelomere[-(150-tract):]+repeat(motif,tract)
            reads.extend([read,rc(read)])
    return reads


def prepare(out):
    if out.exists():
        raise ValueError('Fixture directory exists; preserve prior fixtures and use another path')
    (out/'inputs').mkdir(parents=True)
    rng = random.Random(SEED)
    positive = terminal(rng)
    duplicate = positive[:80] + [positive[80]]*100 + [rc(positive[80])]*100
    one_flank = terminal(rng, ends=1, windows=12)
    # PCR/error copies retain one near-identical anchor family despite raw uniqueness.
    error_copies = duplicate[:80]
    template = positive[80]
    for i in range(20):
        chars = list(template)
        position = 105-20+i  # closest20 nonrepeat bases in the45-bp-tract read
        chars[position] = next(b for b in 'ACGT' if b!=chars[position])
        read = ''.join(chars)
        error_copies += [read,rc(read)]*5
    internal = [dna(rng,35)+repeat(UNKNOWN,80)+dna(rng,35) for _ in range(24)]
    noisy = []
    for read in positive:
        noisy.append(''.join(rng.choice([b for b in 'ACGT' if b!=base]) if rng.random()<0.003 else base for base in read))
    opposite = []
    for _ in range(32):
        s = repeat(UNKNOWN,95)+dna(rng,55)
        opposite += [s,rc(s)]
    fixtures = [
        ('terminal_multiflank',positive,UNKNOWN,'terminal_model','provisional_expected',''),
        ('terminal_noisy',noisy,UNKNOWN,'terminal_model_with_substitutions','sensitivity_challenge',''),
        ('duplicate_junctions',duplicate,UNKNOWN,'one_sequence_cloned','strong_promotion_must_not_follow_raw_counts',''),
        ('single_flank',one_flank,UNKNOWN,'one_boundary_or_shared_subtelomere','diversity_unresolved_not_absence',''),
        ('error_clones',error_copies,UNKNOWN,'near_identical_error_copies','strong_promotion_must_not_follow_raw_counts',''),
        ('internal_bridge',positive+internal,UNKNOWN,'terminal_and_internal_copies_or_internal_array','resolve_internal_counterevidence',''),
        ('internal_only',positive[:80]+internal,UNKNOWN,'internal_repeat_without_clean_terminal_junction','no_geometry_nomination',''),
        ('mixed_orientation',positive+opposite,UNKNOWN,'mixed_boundaries','no_strong_promotion',''),
        ('full_only',positive[:80],UNKNOWN,'full_repeat_no_boundaries','unresolved_not_absence',''),
        ('palindromic',terminal(rng,'AATT'),'AATT','strand_ambiguous','no_geometry_nomination',''),
        ('canonical_ttagac',terminal(rng,'TTAGAC'),'TTAGAC','known_terminal_model','known_observation_preserved',''),
        ('canonical_ttaggt',terminal(rng,'TTAGGT'),'TTAGGT','known_terminal_model','known_observation_preserved',''),
        ('irregular_context',positive,UNKNOWN,'same_reads_with_known_irregular_context','primary_context_before_novel','irregular'),
        ('rna_context',positive,UNKNOWN,'same_reads_with_RNA_metadata_conflict','library_review_before_novel','RNA PolyA'),
        ('hic_context',positive,UNKNOWN,'same_reads_with_HiC_metadata','library_review_before_novel','Arima Hi-C'),
        ('saturated_contexts',terminal(rng,ends=160,windows=2),UNKNOWN,'many_contexts_beyond_sketch_capacity','bounds_and_incomplete_evidence',''),
        ('random_background',[dna(rng,150) for _ in range(3000)],UNKNOWN,'no_introduced_array','no_introduced_signal',''),
    ]
    rows, metadata, assignments, declared = [], [], [], []
    for sid,seqs,motif,truth,expectation,context in fixtures:
        # Equal genome-like nuisance background; no motif enrichment filtering.
        seqs = list(seqs)+[dna(rng,150) for _ in range(500)]
        random.Random(SEED+sum(map(ord,sid))).shuffle(seqs)
        path = (out/'inputs'/f'{sid}.fastq').resolve()
        with path.open('w') as target:
            for i,seq in enumerate(seqs):
                target.write(f'@artificial_{i}\n{seq}\n+\n'+ 'I'*len(seq)+'\n')
        rows.append(dict(sample_id=sid,path1=str(path)))
        metadata.append(dict(sample_id=sid,scientific_name='Artificial evidence fixture',tax_id='2759',
            lineage='Eukaryota',status='synthetic',library_strategy='WGS',library_source='GENOMIC',
            library_layout='SINGLE',library_selection='RANDOM',library_construction_protocol=context if context!='irregular' else '',
            sample_accession='synthetic_'+sid,metagenome='false'))
        assignments.append(dict(sample_id=sid,scientific_name='Artificial evidence fixture',tax_id='2759',lineage='Eukaryota',
            reference_scope='species' if context=='irregular' else '',reference_status='out_of_scope' if context=='irregular' else 'not_curated',
            repeat_architecture='irregular' if context=='irregular' else '',related_motifs='',reference_ids='',reference_evidence='',
            prior_source='',mapping_status='synthetic_context_fixture'))
        declared.append(dict(sample_id=sid,focal_motif=motif,canonical_family=canonical(motif),generating_scenario=truth,
            predeclared_expectation=expectation,records=len(seqs),path=str(path),sha256=digest(path)))
    write_tsv(out/'manifest.tsv',rows)
    write_tsv(out/'metadata.tsv',metadata)
    write_tsv(out/'assignments.tsv',assignments)
    write_json(out/'design.json',dict(seed=SEED,created_utc=datetime.now(timezone.utc).isoformat(),
        generator_sha256=digest(Path(__file__)),fixtures=declared,
        nonidentifiability_control={'input_fixture':'terminal_multiflank','alternative_truth':'An internal array with only one boundary orientation represented can emit these identical short reads.',
            'expectation':'No short-read classifier can distinguish identical observations; provisional must never mean physically confirmed.'},
        limitations=['Synthetic sequence windows and context regression fixtures, not biological sensitivity/specificity or FDR calibration.',
            'Reads include deliberate reverse-complement copies and repeated full-repeat sequences; these are not independent biological replicates.',
            'No genome coverage, physical terminal localization, ALT or literature-wide novelty inferred.',
            'True terminal sequence may share one flank or coexist with interstitial copies; review holds are evidence limits, not false biological absence.']))


def command(argv, log):
    start = time.monotonic()
    result = subprocess.run(argv,cwd=ROOT,text=True,capture_output=True,timeout=600)
    log.parent.mkdir(parents=True,exist_ok=True)
    log.write_text(result.stdout+'\n'+result.stderr)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}); inspect {log}')
    return dict(argv=argv,wall_seconds=time.monotonic()-start,returncode=result.returncode)


def execute(fixtures, out, baseline, updated):
    if out.exists():
        raise ValueError('Preserve previous validation outputs; use a new output directory')
    out.mkdir(parents=True)
    design=json.loads((fixtures/'design.json').read_text())
    for row in design['fixtures']:
        if digest(Path(row['path']))!=row['sha256']:
            raise ValueError('Predeclared fixture changed')
    runs={}
    for label,workspace in [('v040',baseline),('v050',updated)]:
        dest=out/label
        dest.mkdir()
        binary=workspace/'trm'
        if not binary.exists(): binary=workspace/'results/release_v040_20260921/trm'
        scan=[sys.executable,str(workspace/'scripts/batch.py'),str(fixtures/'manifest.tsv'),'--binary',str(binary),
            '--output-dir',str(dest/'scans'),'--jobs','4','--threads','2','--cpu-budget','8','--record-resources','--prior-policy','canonical']
        assessment=[sys.executable,str(workspace/'scripts/report_assessment.py'),'--source',str(dest/'scans'),
            '--metadata',str(fixtures/'metadata.tsv'),'--assignments',str(fixtures/'assignments.tsv'),
            '--catalog',str(workspace/'data/telomere_motifs.tsv'),'--taxon-references',str(workspace/'data/taxon_motif_references.tsv'),
            '--output',str(dest/'assessed'),'--comparison-output',str(dest/'comparison')]
        if label=='v050': assessment+=['--previous-assessment',str(out/'v040/assessed')]
        runs[label]=dict(binary_sha256=digest(binary),scan=command(scan,dest/'scan.log'),
            assessment=command(assessment,dest/'assessment.log'))
    rows=[]
    for fixture in design['fixtures']:
        sid,family=fixture['sample_id'],fixture['canonical_family']
        old=json.loads((out/'v040/assessed/samples'/f'{sid}.json').read_text())
        new=json.loads((out/'v050/assessed/samples'/f'{sid}.json').read_text())
        focal_old=next((r for r in old['candidates'] if r['canonical_family']==family),{})
        focal_new=next((r for r in new['candidates'] if r['canonical_family']==family),{})
        before=json.loads((out/'v040/scans'/sid/'result.json').read_text())
        after=json.loads((out/'v050/scans'/sid/'result.json').read_text())
        def measurements(raw):
            return {c['motif']:{k:c[k] for k in ('reads','full_repeat_reads','repeat_bases','max_tract')}
                for c in raw['candidates']+raw['rejected_candidates']}
        rows.append(dict(sample_id=sid,scenario=fixture['generating_scenario'],expectation=fixture['predeclared_expectation'],
            motif=fixture['focal_motif'],previous_classification=old['classification'],updated_classification=new['classification'],
            previous_focal_status=focal_old.get('candidate_status','not_observed'),updated_focal_status=focal_new.get('candidate_status','not_observed'),
            geometry_nomination=focal_new.get('geometry_nomination',False),stronger_provisional=focal_new.get('novel_candidate',False),
            discovery_measurements_unchanged=measurements(before)==measurements(after),
            confirmation_readiness=json.dumps(focal_new.get('confirmation_readiness',{}),ensure_ascii=False)))
    write_tsv(out/'comparison.tsv',rows)
    write_json(out/'summary.json',dict(runs=runs,fixtures=len(rows),all_discovery_measurements_unchanged=all(r['discovery_measurements_unchanged'] for r in rows),
        synthetic_only=True,biological_FDR_estimated=False,design_sha256=digest(fixtures/'design.json'),results=rows))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--output',type=Path)
    p.add_argument('--baseline-workspace',type=Path)
    p.add_argument('--updated-workspace',type=Path)
    args=p.parse_args()
    if args.prepare_only:
        prepare(args.fixtures.resolve())
    else:
        if not all((args.output,args.baseline_workspace,args.updated_workspace)):
            p.error('Execution requires output and both frozen workspaces')
        execute(args.fixtures.resolve(),args.output.resolve(),args.baseline_workspace.resolve(),args.updated_workspace.resolve())


if __name__=='__main__':
    main()
