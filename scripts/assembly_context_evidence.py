"""Optional, sealed same-source assembly annotations; never changes TRM calls."""
from copy import deepcopy
import csv
import hashlib
import io
from pathlib import Path
import re

from evidence_io import EvidenceIO

SCHEMA_VERSION = 1
SCOPE = 'conditional_saved_mixed_reads_and_same_source_assembly_context'
MAX_BYTES = 8 * 1024 * 1024
MAX_TOTAL_BYTES = 64 * 1024 * 1024
MAX_SAMPLES = 100
MAX_ROWS = 10000
SHA = re.compile(r'^[0-9a-f]{64}$')
FIELDS = ('assembly_context_status', 'assembly_context_read_display',
          'assembly_context_extension_display', 'assembly_context_origin_display',
          'assembly_context_quality_display', 'assembly_context_limit_display')
SUMMARY_FIELDS = ('read_mixed_intervals_selected', 'read_mixed_intervals_preserved',
    'oriented_boundary_events', 'oriented_boundary_hsp_rows', 'oriented_boundary_pairs',
    'scaffolds_reviewed', 'supported_assembly_intervals', 'outward_edge_events',
    'outward_edge_intervals', 'outward_edge_pairs', 'max_supported_assembly_bp',
    'original_fastq_http403', 'normalized_header_http206', 'alternative_same_read_quality_comparisons',
    'source_origin_status', 'compartment_status', 'fcs_status', 'quality_origin_status',
    'physical_terminal_status', 'main_trm_status', 'novel_trm_status', 'independent_validation_status')
SCAFFOLD_FIELDS = ('assembly_accession', 'target', 'length_bp', 'source_organism_label',
    'source_taxon', 'biosample', 'sequence_role', 'assembly_unit', 'assigned_compartment',
    'submitter_voucher_agrees', 'fcs_status', 'fcs_action', 'host_origin_confirmed',
    'compartment_confirmed', 'physical_terminal_confirmed', 'independent_biological_replication')
INTERVAL_FIELDS = ('assembly_accession', 'target', 'grammar_interval_id', 'strand',
    'start0', 'end0', 'length_bp', 'grammar_class', 'left_stop', 'right_stop',
    'outward_side', 'outward_residual_bp', 'outward_scaffold_edge', 'same_source_assembly',
    'physical_terminal_confirmed', 'independent_biological_replication')
EVENT_FIELDS = ('event_id', 'fullread_hsp_id', 'pair_query_id', 'read_query_id', 'target',
    'g_strand_in_reference', 'saved_interval_start0', 'saved_interval_end0', 'saved_interval_bp',
    'grammar_interval_id', 'extended_start0', 'extended_end0', 'extended_exact_grammar_bp',
    'added_left_bp', 'added_right_bp', 'outward_residual_bp', 'outward_stop', 'scaffold_length_bp',
    'retained_qualifying_fullread_hsps', 'retained_qualifying_fullread_target_scaffolds',
    'known_repeat_overlap_in_exterior_bp', 'clean_exterior_confirmed', 'search_capped',
    'same_hsp_inward_mate_witness', 'boundary_outward_direction_consistent',
    'physical_terminal_confirmed', 'independent_biological_replication')
LIMITATIONS = [
    'Previously selected, motif-conditioned saved reads are a conditional subset, not a whole-library detection or absence denominator.',
    'Preserved read intervals, supported assembly intervals and physical chromosome termini are different evidence levels.',
    'The assembly uses the same source specimen and short-read data; extension is not independent biological validation.',
    'Events, mates, pairs, HSP alternatives and overlapping intervals must not be added as independent molecules.',
    'Source organism labels and an FCS report without action rows do not prove scaffold host or nuclear origin.',
    'Exact mixed-grammar lengths and scaffold edges do not confirm physical telomeres, main TRM or novel TRM; no universal 500 bp rule is applied.',
    'Constant Q30 origin and error-rate calibration remain unresolved; no alternative same-read quality was recovered.',
    'The reviewed TTAGG/TTAGGG mixture is separate from the software representative motif, especially in B. rubropictus.',
    'Heuristic search caps remain visible; source hashes establish artifact integrity, not biological truth.',
    'This optional annotation changes no native candidates, counts, ranking, classifications or adequacy.',
]
ABSENT_LIMITATIONS = [
    'No optional assembly context evidence was provided for this sample; all additional measurements remain unassessed.',
    'Missing optional evidence does not indicate zero observations, absence, failed detection or a biological conclusion.',
    'This optional annotation changes no native candidates, counts, ranking, classifications or adequacy.',
]


class EvidenceError(ValueError):
    pass


IO = EvidenceIO(EvidenceError, 'assembly context evidence: ')
SOURCE_ROLES = ('assessment_completion', 'prior_assessment', 'matched_completion',
    'read_summary', 'oriented_contexts', 'origin_completion', 'scaffold_annotations',
    'fcs_review', 'grammar_intervals', 'event_extensions', 'quality_summary')
JSON_ROLES = {'assessment_completion', 'prior_assessment', 'matched_completion', 'origin_completion'}
MEASUREMENTS = ('read_mixed_intervals_selected', 'read_mixed_intervals_preserved',
    'original_fastq_http403', 'normalized_header_http206', 'alternative_same_read_quality_comparisons')
CONFIRMATIONS = ('host_origin', 'compartment', 'physical_terminal', 'main_trm', 'novel_trm',
                 'independent_biological_replication', 'quality_calibrated')
_IDENTITY = ('scientific_name', 'tax_id', 'biosample')
_PAYLOAD_FIELDS = ('schema_version', 'scope', 'sample_id', 'identity', 'original_native_sha256',
    'representative_motif', 'grammar_units', 'same_source_assembly', 'measurements',
    'scaffolds', 'intervals', 'events', 'confirmations', 'review_notes', 'source_artifacts', 'payload_sha256')


def payload_sha256(document):
    return IO.payload_sha256(document)


def _keys(value, keys, where):
    if not isinstance(value, dict) or set(value) != set(keys):
        IO.fail(where + ' has missing or unsupported fields')


def _text(value, where, maximum=4096):
    if not isinstance(value, str) or not value.strip() or len(value) > maximum:
        IO.fail(where + ' must be bounded nonempty text')
    return value


def _int(value, where, maximum=1000000000):
    if type(value) is not int or not 0 <= value <= maximum:
        IO.fail(where + ' must be a bounded nonnegative integer')
    return value


def _bool(value, where, expected=None):
    if type(value) is not bool or (expected is not None and value is not expected):
        IO.fail(where + ' must be the required boolean')
    return value


def _sha(value):
    if not isinstance(value, str) or not SHA.fullmatch(value):IO.fail('invalid SHA256')
    return value


def _tsv_int(value):
    if not isinstance(value, str) or not re.fullmatch(r'0|[1-9][0-9]*', value):IO.fail('invalid source TSV integer')
    return _int(int(value), 'source TSV integer')


def _tsv_bool(value):
    if value not in ('True', 'False'):IO.fail('invalid source TSV boolean')
    return value == 'True'


def _stamp(path):
    try:
        s=path.stat()
        return (str(path.resolve(strict=True)), s.st_dev, s.st_ino, s.st_size,
                s.st_mtime_ns, s.st_ctime_ns, s.st_mode)
    except OSError as error:
        raise EvidenceError('assembly context evidence: source unavailable: ' + str(path)) from error


def _json(data, path):
    return IO.decode_json(data, path, object_pairs_hook=IO.unique_object,
                          parse_constant=IO.reject_constant)[0]


def _table(data, path):
    try:
        reader=csv.DictReader(io.StringIO(data.decode('utf-8')), delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames))!=len(reader.fieldnames):IO.fail('duplicate/missing TSV columns')
        rows=list(reader)
        if len(rows)>MAX_ROWS or any(None in r or any(v is None for v in r.values()) for r in rows):IO.fail('source TSV shape or count invalid')
        return rows
    except (UnicodeError,csv.Error) as error:
        raise EvidenceError('assembly context evidence: invalid TSV: '+str(path)) from error


def _identity(document):
    if not isinstance(document,dict) or not isinstance(document.get('metadata') or {},dict):IO.fail('assessment/metadata must be objects')
    return dict(scientific_name=document.get('scientific_name'),tax_id=str(document.get('tax_id','')),
                biosample=(document.get('metadata') or {}).get('sample_accession'))


def _assessment(document, sample, identity, digest):
    if not isinstance(document,dict) or not isinstance(document.get('provenance') or {},dict):IO.fail('assessment/provenance must be objects')
    if document.get('analysis_complete') is not True or document.get('sample_id') != sample:
        IO.fail('source assessment incomplete or sample differs')
    if document.get('native_scan_version') != '0.8.5':IO.fail('only reviewed native0.8.5 assessments supported')
    if _identity(document)!=identity:IO.fail('sample/taxon/BioSample identity differs')
    if (document.get('provenance') or {}).get('original_discovery_sha256')!=digest:
        IO.fail('original native SHA differs')


def _single(rows, sample, where, extra=None):
    selected=[r for r in rows if r.get('sample_id')==sample and (extra is None or extra(r))]
    if len(selected)!=1:IO.fail(where+' requires exactly one sample row')
    return selected[0]


def _derive(sample, sources):
    try:
        return _derive_checked(sample,sources)
    except (KeyError,TypeError,AttributeError) as error:
        raise EvidenceError('assembly context evidence: detailed source shape is invalid') from error


def _derive_checked(sample, sources):
    """Rebuild the small exported annotations from their sealed tabular sources."""
    prior=sources['prior_assessment'];identity=_identity(prior)
    native=(prior.get('provenance') or {}).get('original_discovery_sha256');_sha(native)
    _assessment(prior,sample,identity,native)
    read=_single(sources['read_summary'],sample,'gapped read summary',lambda r:r.get('method')=='gapped')
    quality=_single(sources['quality_summary'],sample,'quality review')
    fcs=_single(sources['fcs_review'],sample,'FCS review')
    if not _tsv_bool(fcs['report_present']) or _tsv_int(fcs['report_data_rows'])!=0 or _tsv_bool(fcs['positive_scaffold_taxonomy_report_available']):
        IO.fail('reviewed FCS no-action scope differs')
    if quality['quality_origin']!='unresolved' or _tsv_bool(quality['alternative_same_read_quality_compared']):IO.fail('reviewed quality scope differs')
    if (not _tsv_bool(quality['all_observed_bases_q30']) or
        _tsv_int(quality['local_head_complete_reads'])<1 or _tsv_int(quality['local_head_bases'])<1 or
        _tsv_int(quality['normalized_decoded_fastq_records'])!=0):
        IO.fail('reviewed constant-Q30 prefix evidence scope differs')
    measurement=dict(read_mixed_intervals_selected=_tsv_int(read['strict_mixed_intervals']),
        read_mixed_intervals_preserved=_tsv_int(read['identical_observed_mixed_intervals']),
        original_fastq_http403=_tsv_int(quality['original_submitted_http403']),
        normalized_header_http206=_tsv_int(quality['normalized_http206']),
        alternative_same_read_quality_comparisons=0)
    scaffolds=[]
    for r in sources['scaffold_annotations']:
        if r['sample_id']!=sample:continue
        if r['biosample']!=identity['biosample'] or r['fcs_report_status']!='present_no_action_rows':IO.fail('scaffold identity/FCS scope differs')
        for key in ('host_origin_independently_confirmed','physical_terminal_confirmed','independent_biological_replication'):
            if _tsv_bool(r[key]):IO.fail('source confirmation outside reviewed scope')
        scaffolds.append(dict(assembly_accession=r['assembly_accession'],target=r['scaffold_accession_version'],length_bp=_tsv_int(r['length_bp']),
            source_organism_label=r['genbank_source_organism_label'],source_taxon=r['genbank_taxon_dbxref'],biosample=r['biosample'],
            sequence_role=r['assembly_sequence_role'],assembly_unit=r['assembly_unit'],assigned_compartment='unresolved',
            submitter_voucher_agrees=_tsv_bool(r['submitter_seqid_voucher_prefix_agrees']),fcs_status=r['fcs_report_status'],
            fcs_action=r['fcs_scaffold_action'],host_origin_confirmed=False,compartment_confirmed=False,
            physical_terminal_confirmed=False,independent_biological_replication=False))
    bytarget={r['target']:r for r in scaffolds}
    if len(bytarget)!=len(scaffolds):IO.fail('duplicate source scaffold')
    oriented={};extensions={}
    for role,target in (('oriented_contexts',oriented),('event_extensions',extensions)):
        for r in sources[role]:
            if r['sample_id']!=sample:continue
            key=(r['event_id'],r['fullread_hsp_id'])
            if key in target:IO.fail('duplicate source event/HSP')
            target[key]=r
    if set(oriented)!=set(extensions):IO.fail('oriented event/extension HSP binding differs')
    intervals=[];events=[];needed={r['grammar_interval_id'] for r in extensions.values()}
    for r in sources['grammar_intervals']:
        if r['sample_id']!=sample or r['grammar_interval_id'] not in needed:continue
        if r['target'] not in bytarget:IO.fail('interval target not reviewed')
        a,b=_tsv_int(r['start0']),_tsv_int(r['end0']);length=bytarget[r['target']]['length_bp']
        strand=r['strand'];residual=length-b if strand=='+' else a
        intervals.append(dict(assembly_accession=r['assembly_accession'],target=r['target'],grammar_interval_id=r['grammar_interval_id'],
            strand=strand,start0=a,end0=b,length_bp=_tsv_int(r['length_bp']),grammar_class=r['grammar_class'],left_stop=r['left_stop'],right_stop=r['right_stop'],
            outward_side='right' if strand=='+' else 'left',outward_residual_bp=residual,outward_scaffold_edge=residual==0,
            same_source_assembly=True,physical_terminal_confirmed=False,independent_biological_replication=False))
    if {r['grammar_interval_id'] for r in intervals}!=needed or len(intervals)!=len(needed):IO.fail('supported interval source identity differs')
    for key,r in sorted(extensions.items()):
        o=oriented[key]
        for k in ('pair_query_id','target'):
            if r[k]!=o[k]:IO.fail('event/HSP source geometry identity differs')
        for current,previous in (('saved_interval_start0','mixed_interval_start0'),('saved_interval_end0','mixed_interval_end0'),('scaffold_length_bp','scaffold_length_bp')):
            if _tsv_int(r[current])!=_tsv_int(o[previous]):IO.fail('event/HSP source coordinates differ')
        if r['g_strand_in_reference']!=('+' if o['repeat_outward_side']=='right' else '-'):
            IO.fail('saved outward orientation differs')
        if not o['inward_mate_witness_hsp_ids']:IO.fail('same-HSP mate witness missing')
        for key in ('physical_terminal_confirmed','independent_evidence'):
            if _tsv_bool(r[key]):IO.fail('source extension confirmation outside reviewed scope')
        row={k:r[k] for k in ('event_id','fullread_hsp_id','pair_query_id','target','g_strand_in_reference','grammar_interval_id','outward_stop')}
        for k in ('saved_interval_start0','saved_interval_end0','saved_interval_bp','extended_start0','extended_end0','extended_exact_grammar_bp','added_left_bp','added_right_bp','outward_residual_bp','scaffold_length_bp'):
            row[k]=_tsv_int(r[k])
        row.update(read_query_id=o['read_query_id'],retained_qualifying_fullread_hsps=_tsv_int(o['retained_qualifying_fullread_hsps']),
            retained_qualifying_fullread_target_scaffolds=_tsv_int(o['retained_qualifying_fullread_target_scaffolds']),
            known_repeat_overlap_in_exterior_bp=_tsv_int(o['known_repeat_overlap_in_exterior_bp']),clean_exterior_confirmed=_tsv_bool(o['clean_exterior_confirmed']),
            search_capped=any(_tsv_bool(o[k]) for k in ('fullread_retention_truncated','fullread_target_search_cap_possibly_reached','fullread_per_target_search_cap_possibly_reached')),
            same_hsp_inward_mate_witness=True,boundary_outward_direction_consistent=True,physical_terminal_confirmed=False,independent_biological_replication=False)
        events.append(row)
    if len({r['event_id'] for r in events})!=_tsv_int(read['oriented_same_hsp_inward_mate_join_events']) or len({r['pair_query_id'] for r in events})!=_tsv_int(read['oriented_same_hsp_inward_mate_join_pairs']):
        IO.fail('oriented event/pair counts differ from saved mapping summary')
    return dict(identity=identity,original_native_sha256=native,representative_motif=prior.get('detected_motif'),measurements=measurement,
                scaffolds=sorted(scaffolds,key=lambda r:r['target']),intervals=sorted(intervals,key=lambda r:r['grammar_interval_id']),events=events)


def _validate_payload(doc):
    _keys(doc,_PAYLOAD_FIELDS,'payload')
    if type(doc['schema_version']) is not int or doc['schema_version']!=SCHEMA_VERSION or doc['scope']!=SCOPE:IO.fail('unsupported payload schema/scope')
    if payload_sha256(doc)!=_sha(doc['payload_sha256']):IO.fail('payload seal differs')
    _text(doc['sample_id'],'sample ID');_sha(doc['original_native_sha256'])
    _keys(doc['identity'],_IDENTITY,'identity')
    for k,v in doc['identity'].items():_text(v,'identity.'+k)
    if not re.fullmatch('[ACGT]{1,64}',_text(doc['representative_motif'],'representative motif')):IO.fail('invalid representative motif')
    if doc['grammar_units']!=['TTAGG','TTAGGG']:IO.fail('unsupported reviewed grammar')
    _bool(doc['same_source_assembly'],'same source assembly',True)
    _keys(doc['confirmations'],CONFIRMATIONS,'confirmations')
    for k,v in doc['confirmations'].items():_bool(v,k,False)
    if not isinstance(doc['review_notes'],list) or len(doc['review_notes'])>16:IO.fail('review notes bound exceeded')
    for note in doc['review_notes']:_text(note,'review note',2000)
    _keys(doc['measurements'],MEASUREMENTS,'measurements')
    for k,v in doc['measurements'].items():_int(v,k)
    m=doc['measurements']
    if m['read_mixed_intervals_preserved']>m['read_mixed_intervals_selected'] or m['alternative_same_read_quality_comparisons']!=0:IO.fail('measurement scope/count differs')
    bool_fields={'submitter_voucher_agrees','host_origin_confirmed','compartment_confirmed','physical_terminal_confirmed','independent_biological_replication',
        'outward_scaffold_edge','same_source_assembly','clean_exterior_confirmed','search_capped','same_hsp_inward_mate_witness','boundary_outward_direction_consistent'}
    int_fields={'length_bp','start0','end0','outward_residual_bp','saved_interval_start0','saved_interval_end0','saved_interval_bp',
        'extended_start0','extended_end0','extended_exact_grammar_bp','added_left_bp','added_right_bp','scaffold_length_bp',
        'retained_qualifying_fullread_hsps','retained_qualifying_fullread_target_scaffolds','known_repeat_overlap_in_exterior_bp'}
    for kind,fields in (('scaffolds',SCAFFOLD_FIELDS),('intervals',INTERVAL_FIELDS),('events',EVENT_FIELDS)):
        if not isinstance(doc[kind],list) or not 1<=len(doc[kind])<=MAX_ROWS:IO.fail(kind+' rows outside supported bounds')
        for row in doc[kind]:
            _keys(row,fields,kind+' row')
            for k,v in row.items():
                if k in bool_fields:_bool(v,k)
                elif k in int_fields:_int(v,k)
                else:_text(v,k)
            for key in ('physical_terminal_confirmed','independent_biological_replication'):_bool(row[key],key,False)
    scaffolds={r['target']:r for r in doc['scaffolds']}
    if len(scaffolds)!=len(doc['scaffolds']):IO.fail('duplicate scaffold')
    for row in scaffolds.values():
        if row['biosample']!=doc['identity']['biosample'] or row['length_bp']==0:IO.fail('scaffold sample/length differs')
        if row['assigned_compartment']!='unresolved' or row['fcs_status']!='present_no_action_rows':IO.fail('unsupported origin/FCS claim')
        _bool(row['host_origin_confirmed'],'host confirmation',False);_bool(row['compartment_confirmed'],'compartment confirmation',False)
    intervals={r['grammar_interval_id']:r for r in doc['intervals']}
    if len(intervals)!=len(doc['intervals']):IO.fail('duplicate grammar interval')
    for row in intervals.values():
        scaffold=scaffolds.get(row['target'])
        if scaffold is None or row['assembly_accession']!=scaffold['assembly_accession']:IO.fail('interval scaffold binding differs')
        if row['strand'] not in ('+','-') or not 0<=row['start0']<row['end0']<=scaffold['length_bp'] or row['length_bp']!=row['end0']-row['start0']:IO.fail('interval coordinates differ')
        if row['grammar_class']!='mixed_grammar_neither_pure_phase_compatible':IO.fail('supported interval must retain mixed grammar')
        if row['outward_side']!=('right' if row['strand']=='+' else 'left'):IO.fail('interval outward orientation differs')
        expected=scaffold['length_bp']-row['end0'] if row['strand']=='+' else row['start0']
        if row['outward_residual_bp']!=expected or row['outward_scaffold_edge']!=(expected==0):IO.fail('interval edge distance differs')
        _bool(row['same_source_assembly'],'same source interval',True)
    seen=set();event_bindings={}
    for row in doc['events']:
        key=(row['event_id'],row['fullread_hsp_id'])
        if key in seen:IO.fail('duplicate event/HSP')
        seen.add(key)
        binding=(row['pair_query_id'],row['read_query_id'])
        if row['event_id'] in event_bindings and event_bindings[row['event_id']]!=binding:IO.fail('event read/pair binding differs')
        event_bindings[row['event_id']]=binding
        interval=intervals.get(row['grammar_interval_id'])
        if interval is None or row['target']!=interval['target'] or row['g_strand_in_reference']!=interval['strand']:IO.fail('event interval binding differs')
        a,b=row['saved_interval_start0'],row['saved_interval_end0']
        if not interval['start0']<=a<b<=interval['end0'] or row['saved_interval_bp']!=b-a:IO.fail('saved interval outside assembly interval')
        expected=dict(extended_start0=interval['start0'],extended_end0=interval['end0'],extended_exact_grammar_bp=interval['length_bp'],
            added_left_bp=a-interval['start0'],added_right_bp=interval['end0']-b,outward_residual_bp=interval['outward_residual_bp'],scaffold_length_bp=scaffolds[row['target']]['length_bp'])
        if any(row[k]!=v for k,v in expected.items()):IO.fail('event extension arithmetic differs')
        if row['outward_stop']!=(interval['right_stop'] if interval['strand']=='+' else interval['left_stop']):IO.fail('outward stopping reason differs')
        for field in ('same_hsp_inward_mate_witness','boundary_outward_direction_consistent'):_bool(row[field],field,True)
        if row['retained_qualifying_fullread_hsps']<1 or row['retained_qualifying_fullread_target_scaffolds']<1:IO.fail('supporting HSP count missing')
    if set(intervals)!={r['grammar_interval_id'] for r in doc['events']} or set(scaffolds)!={r['target'] for r in doc['events']}:IO.fail('unused/missing supported intervals/scaffolds')
    if len(event_bindings)>m['read_mixed_intervals_preserved']:IO.fail('selected events exceed preserved interval denominator')


def _bind_sources(sources, paths, hashes):
    groups={'assessment_completion':('prior_assessment',),
            'matched_completion':('read_summary','oriented_contexts'),
            'origin_completion':('scaffold_annotations','fcs_review','grammar_intervals','event_extensions','quality_summary')}
    for receipt_role,roles in groups.items():
        receipt=sources[receipt_role]
        if not isinstance(receipt,dict) or receipt.get('status')!='complete' or not isinstance(receipt.get('artifacts_sha256'),dict):IO.fail('source completion receipt invalid')
        base=paths[receipt_role].parent
        for role in roles:
            try:key=str(paths[role].resolve().relative_to(base.resolve()))
            except ValueError:IO.fail('source outside its completion receipt root')
            if receipt['artifacts_sha256'].get(key)!=hashes[role]:IO.fail('source differs from completed artifact seal')


def load_evidence(path):
    if path is None:return None
    path=Path(path).absolute();cache={};stamps={};total=0
    def read(source):
        nonlocal total
        source=Path(source).absolute();key=str(source)
        if key not in cache:
            before=_stamp(source);data=IO.read_bytes(source,max_bytes=MAX_BYTES)
            after=_stamp(source)
            if before!=after:IO.fail('source changed around bounded read')
            total+=len(data)
            if total>MAX_TOTAL_BYTES:IO.fail('total evidence bytes exceed bound')
            cache[key]=(data,hashlib.sha256(data).hexdigest());stamps[key]=after
        return cache[key]
    raw,digest=read(path);manifest=_json(raw,path)
    _keys(manifest,('schema_version','scope','samples','payload_sha256'),'manifest')
    if type(manifest['schema_version']) is not int or manifest['schema_version']!=SCHEMA_VERSION or manifest['scope']!=SCOPE:IO.fail('manifest schema/scope unsupported')
    if payload_sha256(manifest)!=_sha(manifest['payload_sha256']):IO.fail('manifest seal differs')
    if not isinstance(manifest['samples'],list) or not 1<=len(manifest['samples'])<=MAX_SAMPLES:IO.fail('manifest sample count outside bound')
    result=dict(manifest_sha256=digest,samples={},paths=[],stamps=stamps)
    for row in manifest['samples']:
        _keys(row,('sample_id','identity','original_native_sha256','evidence'),'manifest sample')
        sample=_text(row['sample_id'],'sample ID')
        if sample in result['samples']:IO.fail('duplicate sample ID')
        _keys(row['identity'],_IDENTITY,'manifest identity')
        for value in row['identity'].values():_text(value,'manifest identity')
        _sha(row['original_native_sha256']);_keys(row['evidence'],('path','sha256'),'evidence pin')
        payload_path=Path(_text(row['evidence']['path'],'payload path'));payload_path=payload_path if payload_path.is_absolute() else path.parent/payload_path
        data,actual=read(payload_path)
        if actual!=_sha(row['evidence']['sha256']):IO.fail('sample evidence hash differs')
        doc=_json(data,payload_path);_validate_payload(doc)
        for key in ('sample_id','identity','original_native_sha256'):
            if row[key]!=doc[key]:IO.fail('manifest/payload identity differs')
        pins=doc['source_artifacts']
        if not isinstance(pins,list) or len(pins)!=len(SOURCE_ROLES):IO.fail('source artifact count differs')
        sources={};paths={};hashes={};seen=set()
        for pin in pins:
            _keys(pin,('role','path','sha256'),'source pin');role=pin['role']
            if not isinstance(role,str) or role not in SOURCE_ROLES or role in sources:IO.fail('missing/duplicate/unsupported source role')
            source=Path(_text(pin['path'],'source path'));source=source if source.is_absolute() else payload_path.parent/source
            if str(source.resolve()) in seen:IO.fail('duplicate source path')
            seen.add(str(source.resolve()))
            expected_suffix='.json' if role in JSON_ROLES else '.tsv'
            if source.suffix!=expected_suffix:IO.fail('source format not supported')
            data,actual=read(source)
            if actual!=_sha(pin['sha256']):IO.fail('detailed source hash differs')
            sources[role]=_json(data,source) if role in JSON_ROLES else _table(data,source)
            paths[role]=source.absolute();hashes[role]=actual
        if set(sources)!=set(SOURCE_ROLES):IO.fail('source roles incomplete')
        _bind_sources(sources,paths,hashes)
        rebuilt=_derive(sample,sources)
        for key,expected in rebuilt.items():
            if doc[key]!=expected:IO.fail('exported '+key+' differs from detailed sealed sources')
        result['samples'][sample]=dict(document=doc,path=str(payload_path.absolute()),sha256=row['evidence']['sha256'],
            source_artifacts=[dict(role=role,path=str(paths[role]),sha256=hashes[role]) for role in SOURCE_ROLES])
    result['paths']=sorted(stamps)
    return result


def source_dependencies(path):
    return [] if path is None else load_evidence(path)['paths']


def _summary(doc):
    events=doc['events'];intervals=doc['intervals'];edge=[r for r in events if r['outward_residual_bp']==0]
    summary=dict(doc['measurements'])
    summary.update(oriented_boundary_events=len({r['event_id'] for r in events}),oriented_boundary_hsp_rows=len(events),
        oriented_boundary_pairs=len({r['pair_query_id'] for r in events}),scaffolds_reviewed=len(doc['scaffolds']),
        supported_assembly_intervals=len(intervals),outward_edge_events=len({r['event_id'] for r in edge}),
        outward_edge_intervals=len({r['grammar_interval_id'] for r in edge}),outward_edge_pairs=len({r['pair_query_id'] for r in edge}),
        max_supported_assembly_bp=max(r['length_bp'] for r in intervals),source_origin_status='not_confirmed',compartment_status='unresolved',
        fcs_status='present_no_action_rows_not_host_proof',quality_origin_status='unresolved_no_alternative_same_read_quality',
        physical_terminal_status='not_confirmed',main_trm_status='not_confirmed',novel_trm_status='not_confirmed',independent_validation_status='not_confirmed')
    return {key:summary[key] for key in SUMMARY_FIELDS}


def build_evidence(sample, document, loaded):
    result=dict(status='not_provided',display=dict(zip(FIELDS,('Assembly 추가 근거 미제공','저장 read 위치 미평가','Assembly 반복 연장 미평가',
        'Scaffold 기원·FCS 미평가','대체 품질 기원 미평가','별도 검토 자료가 제공되지 않음'))),
        summary={key:None for key in SUMMARY_FIELDS},scaffolds=[],intervals=[],events=[],scope=SCOPE,
        limitations=list(ABSENT_LIMITATIONS),changes_classification_or_adequacy=False,source_artifacts=[])
    if loaded is None or sample not in loaded['samples']:return result
    for path,stamp in loaded['stamps'].items():
        if _stamp(Path(path))!=stamp:IO.fail('loaded evidence is stale')
    row=loaded['samples'][sample];doc=row['document'];_validate_payload(doc)
    _assessment(document,sample,doc['identity'],doc['original_native_sha256'])
    if document.get('detected_motif')!=doc['representative_motif']:IO.fail('representative motif differs; overlay cannot replace it')
    summary=_summary(doc)
    display=dict(assembly_context_status='파일·표본 연결 무결성 검증됨 · 생물학적 확정 아님',
        assembly_context_read_display=f"선정 mixed {summary['read_mixed_intervals_selected']}구간 중 assembly 정렬에서 구간 전체 동일 보존 {summary['read_mixed_intervals_preserved']}; 방향 일치 경계 {summary['oriented_boundary_events']}건",
        assembly_context_extension_display=f"같은 원천 assembly {summary['supported_assembly_intervals']}구간, 최대 {summary['max_supported_assembly_bp']} bp; outward scaffold 끝 {summary['outward_edge_events']}event/{summary['outward_edge_intervals']}배열 · 물리적 말단 미확인",
        assembly_context_origin_display='등록 source label·FCS action 0행 · 숙주/핵/소기관 기원 미확정',
        assembly_context_quality_display=f"고정 Q30 기원 미해결 · 대체 same-read 품질 비교 {summary['alternative_same_read_quality_comparisons']}건",
        assembly_context_limit_display=f"대표 motif {doc['representative_motif']}와 TTAGG/TTAGGG 혼합 근거 별도 · 같은 원천, 독립 검증/physical/main/novel 미확인" + ('; '+'; '.join(doc['review_notes']) if doc['review_notes'] else ''))
    result.update(status='sources_validated',display=display,summary=summary,limitations=list(LIMITATIONS),
        scaffolds=deepcopy(doc['scaffolds']),intervals=deepcopy(doc['intervals']),events=deepcopy(doc['events']),
        source_artifacts=deepcopy(row['source_artifacts']),manifest_sha256=loaded['manifest_sha256'],evidence_sha256=row['sha256'],
        original_native_sha256=doc['original_native_sha256'],grammar_units=list(doc['grammar_units']),representative_motif=doc['representative_motif'],
        same_source_assembly=True,confirmations=deepcopy(doc['confirmations']),review_notes=list(doc['review_notes']))
    return result
