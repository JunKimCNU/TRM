#!/usr/bin/env python3
"""Join immutable TRM observations and library diagnostics into a review report."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from html import escape
import json
import os
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
LABELS = {
    'RNA_component_supported': 'RNA 성분의 강한 증거',
    'metadata_conflict_unresolved': 'RNA/DNA 미확정',
    'HiC_supported': 'Hi-C 정황 지지',
    'HiC_annotation_only': 'Hi-C 표기: 추가 확인 필요',
    'WGS_context_compatible': '일반 WGS 정황과 부합',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def load_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_tsv(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def interpretation(sequence, alignment=None, review=None):
    """Keep observation, identity review and genomic adequacy strictly separate."""
    if sequence['sample_id'] == 'ERR13623736' and alignment is not None and review is not None:
        if review.get('review_label') != 'RNA_supported' or review.get('library_purity_estimated') is not False:
            raise ValueError('JU765 reviewed RNA-component evidence/purity guard is not satisfied')
        count = alignment['counts']['strong_spliced_primary_records']
        denom = alignment['counts']['high_mapq_primary_mapped_records']
        if not (0 < count <= denom):
            raise ValueError('Invalid strong-splice count/denominator')
        fraction = count / denom
        if abs(fraction - alignment['high_mapq_strong_spliced_fraction']) > 1e-9:
            raise ValueError('Strong-splice fraction does not match its denominator')
        code = 'RNA_component_supported'
        proof = (f"강한 splice reads {count:,}/{denom:,} ({fraction:.1%}); "
                 f"정렬 ACGT 염기의 exon 비율 {alignment['high_mapq_aligned_reference_ACGT_fractions']['exonic']:.1%}.")
        impact = '이 run 단독으로 genomic WGS가 충분하다고 판단하거나 TRM 부재·novel telomeric sequence를 주장할 수 없습니다.'
        action = '동일 isolate/species의 독립적인 genomic DNA WGS를 확인하고 재검증합니다. RNA/DNA 혼합비와 실제 library kit는 미확정입니다.'
        evidence_level = 'strong_splice_and_annotation_evidence'
    elif sequence['role'] == 'RNA_PolyA_metadata_conflict':
        code = 'metadata_conflict_unresolved'
        proof = (f"WGS/GENOMIC/RANDOM 대 RNA PolyA 표기 충돌. 정확한 말단 SL hit "
                 f"{sequence['full_SL_anchored_records']}건; exact-isolate reference 미확인.")
        impact = '이전 novel TRM은 관측된 반복 후보로 보존합니다. 라이브러리 정체가 확인되기 전에는 genomic telomere 해석을 확정하지 않습니다.'
        action = '실제 protocol·molecular source를 확인하거나 matching genome/GFF 또는 독립 genomic DNA WGS를 확보합니다.'
        evidence_level = 'unresolved_metadata_conflict'
    elif sequence['role'] == 'HiC_protocol_annotation':
        markov = sequence['ligation_internal_enrichment_vs_markov1']
        iid = sequence['ligation_internal_enrichment_vs_iid']
        supported = markov is not None and iid is not None and markov > 1 and iid > 1
        code = 'HiC_supported' if supported else 'HiC_annotation_only'
        proof = (f"Hi-C - Arima v2 표기; GATCGATC 포함 reads {sequence['ligation_read_fraction']:.2%}; "
                 f"내부 motif/1차 Markov 기대치 {markov:.2f}×.") if markov is not None else 'Hi-C protocol 표기; 조성 보정 수치 미확정.'
        impact = '기존 TTAGGC 검출 관측은 유지합니다. Hi-C는 genomic DNA이지만 균일 WGS가 아니므로 미검출·genome size·검출력 판단에는 주의가 필요합니다.'
        action = '알려진 motif 검출 기록으로 보존합니다. 부재/novel 가능성을 판단해야 할 때 독립적인 일반 genomic WGS로 확인합니다.'
        evidence_level = 'protocol_and_internal_ligation_enrichment' if supported else 'protocol_annotation_only'
    elif sequence['role'] == 'ordinary_genomic_comparison':
        code = 'WGS_context_compatible'
        proof = (f"공개 WGS/GENOMIC 정황과 부합; GATCGATC 포함 reads {sequence['ligation_read_fraction']:.2%}; "
                 f"내부 motif/1차 Markov 기대치 {sequence['ligation_internal_enrichment_vs_markov1']:.2f}×.")
        impact = '기존 알려진 motif 검출 관측을 유지합니다. 이 비교로 library 순도나 genome coverage·telomere 검출력이 증명되지는 않습니다.'
        action = '비교용 genomic 자료로 활용합니다. 미검출을 해석할 경우 library와 sampling 적합성을 별도로 평가합니다.'
        evidence_level = 'metadata_and_descriptive_comparison_only'
    else:
        raise ValueError(f"Unreviewed sample role: {sequence['role']}")
    return {'library_review_code': code, 'library_review_label': LABELS[code],
            'evidence_level': evidence_level, 'key_evidence': proof,
            'TRM_interpretation_impact': impact, 'next_step': action}


def build_rows(report_dir, metadata_path, assessment_dir):
    signatures_path = report_dir / 'sequence_evidence.json'
    inspection_path = report_dir / 'prefix_inspection/inspection.json'
    review_path = report_dir / 'JU765_alignment/evidence_review.json'
    alignment_path = report_dir / 'JU765_alignment/ERR13623736_prefix/summary.json'
    references_path = report_dir / 'reference_inventory.tsv'
    signatures, inspection = read_json(signatures_path), read_json(inspection_path)
    review, alignment = read_json(review_path), read_json(alignment_path)
    metadata = {row['sample_id']: row for row in load_tsv(metadata_path)}
    if signatures['source_inspection_sha256'] != sha(inspection_path):
        raise ValueError('Sequence signatures used a different prefix inspection')
    prefix = Counter()
    prefix_files = Counter()
    for row in inspection['results']:
        prefix[row['sample_id']] += row['counts']['records']
        prefix_files[row['sample_id']] += 1
    if set(prefix) != {r['sample_id'] for r in signatures['samples']}:
        raise ValueError('Prefix and sequence-audit sample sets differ')
    ju_prefix = {row['extracted_sha256'] for row in inspection['results'] if row['sample_id'] == 'ERR13623736'}
    if ju_prefix != set(alignment['input_sha256']) or len(ju_prefix) != 2:
        raise ValueError('JU765 alignment does not use exactly the two inspected prefixes')
    if alignment['counts']['primary_records'] != prefix['ERR13623736']:
        raise ValueError('JU765 alignment/prefix read counts differ')
    dependencies = [signatures_path, inspection_path, review_path, alignment_path, references_path, metadata_path]
    rows = []
    order = {'RNA_component_supported': 0, 'metadata_conflict_unresolved': 1,
             'HiC_supported': 2, 'HiC_annotation_only': 2, 'WGS_context_compatible': 3}
    for seq in signatures['samples']:
        sid = seq['sample_id']
        old_path = assessment_dir / f'{sid}.json'
        old = read_json(old_path)
        meta = metadata[sid]
        if meta.get('assessment_sha256') and sha(old_path) != meta['assessment_sha256']:
            raise ValueError(f'Previous assessment hash changed for {sid}')
        if prefix[sid] != seq['prefix_records'] or prefix_files[sid] != seq['prefix_files']:
            raise ValueError(f'Prefix accounting mismatch for {sid}')
        detail_path = 'JU765_alignment/README.md' if sid == 'ERR13623736' else 'sequence_evidence_notes.md'
        reviewed = interpretation(seq, alignment if sid == 'ERR13623736' else None,
                                  review if sid == 'ERR13623736' else None)
        row = {'sample_id': sid, 'scientific_name': old['scientific_name'], 'tax_id': old['tax_id'],
               'previous_TRM_classification': old['classification'], 'previous_TRM_label': old['result'],
               'previous_detected_motif': old['detected_motif'] or '',
               'previous_full_read_records': old['input_summary']['records'],
               'previous_full_bases': old['input_summary']['bases'],
               'inspected_prefix_files': prefix_files[sid], 'inspected_prefix_read_records': prefix[sid],
               'inspected_prefix_bases': seq['prefix_bases'],
               'sampling': 'first_100000_records_per_available_file_not_random',
               'library_strategy': meta['library_strategy'], 'library_source': meta['library_source'],
               'library_selection': meta['library_selection'],
               'library_construction_protocol': meta['library_construction_protocol'],
               **reviewed,
               'full_SL_anchored_records': seq['full_SL_anchored_records'],
               'polyA_or_rc_polyT_ge20_records': seq['polyA_or_rc_polyT_ge20_records'],
               'ligation_GATCGATC_read_fraction': seq['ligation_read_fraction'],
               'ligation_internal_enrichment_vs_iid': seq['ligation_internal_enrichment_vs_iid'],
               'ligation_internal_enrichment_vs_markov1': seq['ligation_internal_enrichment_vs_markov1'],
               'strong_spliced_read_records': alignment['counts']['strong_spliced_primary_records'] if sid == 'ERR13623736' else '',
               'splice_denominator_high_mapq_primary_mapped_records': alignment['counts']['high_mapq_primary_mapped_records'] if sid == 'ERR13623736' else '',
               'strong_spliced_fraction': alignment['high_mapq_strong_spliced_fraction'] if sid == 'ERR13623736' else '',
               'exonic_fraction_of_high_mapq_aligned_ACGT_bases': alignment['high_mapq_aligned_reference_ACGT_fractions']['exonic'] if sid == 'ERR13623736' else '',
               'library_purity_estimated': False, 'molecular_mixture_proportion_estimated': False,
               'genomic_TRM_absence_established': False, 'genomic_detection_adequacy_established_by_this_audit': False,
               'hard_motif_filter_applied': False, 'previous_TRM_outputs_modified': False,
               'previous_assessment_path': os.path.relpath(old_path, report_dir), 'previous_assessment_sha256': sha(old_path),
               'metadata_source_url': meta['library_source_url'], 'detailed_evidence_path': detail_path,
               'SRA_genomic_search_url': 'https://www.ncbi.nlm.nih.gov/sra/?term=' + quote(
                   f'txid{old["tax_id"]}[Organism:exp] AND WGS[Strategy] AND genomic[Source]')}
        rows.append(row)
        dependencies.append(old_path)
    rows.sort(key=lambda row: (order[row['library_review_code']], row['sample_id']))
    return rows, signatures, inspection, alignment, review, dependencies


def render_html(rows, summary, alignment, review):
    esc = lambda value: escape(str(value), quote=True)
    cards = ''.join(f'<div class="card {code}"><strong>{summary["review_counts"].get(code, 0)}</strong><span>{esc(label)}</span></div>'
                    for code, label in LABELS.items() if code != 'HiC_annotation_only')
    rendered = []
    for row in rows:
        motif = row['previous_detected_motif'] or '대표 motif 미확정'
        sources = (f'<a href="{esc(row["previous_assessment_path"])}">이전 원본 JSON</a> · '
                   f'<a href="{esc(row["metadata_source_url"])}">공개 metadata</a> · '
                   f'<a href="{esc(row["detailed_evidence_path"])}">검증 상세</a>')
        next_link = f'<a href="{esc(row["SRA_genomic_search_url"])}">동일 TaxID의 SRA WGS 검색</a>'
        rendered.append(f'''<tr data-kind="{esc(row['library_review_code'])}">
<td><b>{esc(row['sample_id'])}</b><br>{esc(row['scientific_name'])}<small>검사: {row['inspected_prefix_files']} files · {row['inspected_prefix_read_records']:,} reads</small></td>
<td><span class="motif">{esc(motif)}</span><small>{esc(row['previous_TRM_label'])}</small></td>
<td><span class="badge {esc(row['library_review_code'])}">{esc(row['library_review_label'])}</span><small>{esc(row['library_construction_protocol'] or 'Protocol 상세 미기재')}</small></td>
<td>{esc(row['key_evidence'])}<details><summary>근거 파일 / 수치</summary><p>정확한 말단 SL: {row['full_SL_anchored_records']}건. 말단 A/T ≥20 bp: {row['polyA_or_rc_polyT_ge20_records']}건. 내부 ligation motif 기대치 대비: IID {row['ligation_internal_enrichment_vs_iid']:.2f}× / Markov {row['ligation_internal_enrichment_vs_markov1']:.2f}×.</p>{sources}</details></td>
<td>{esc(row['TRM_interpretation_impact'])}</td><td>{esc(row['next_step'])}<small>{next_link}</small></td></tr>''')
    options = ''.join(f'<option value="{esc(code)}">{esc(label)}</option>' for code, label in LABELS.items()
                      if summary['review_counts'].get(code))
    return f'''<!doctype html>
<html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>TRM 자료의 라이브러리 정체 검토 · 13 samples</title>
<style>
:root{{color-scheme:light;--ink:#172b3a;--muted:#556675;--line:#dce4ea;--blue:#176b94}}*{{box-sizing:border-box}}body{{margin:0;background:#f4f7fa;color:var(--ink);font:15px/1.65 system-ui,-apple-system,'Noto Sans KR',sans-serif}}main{{max-width:1540px;margin:auto;padding:32px 24px 60px}}h1{{font-size:30px;line-height:1.3;margin:0 0 14px}}h2{{font-size:21px;margin:0 0 12px}}p{{margin:8px 0 14px}}a{{color:#075b87;text-underline-offset:3px}}.eyebrow,small{{color:var(--muted)}}.lead{{max-width:1060px;font-size:17px}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:24px 0}}.card{{background:white;border:1px solid var(--line);border-top:4px solid #658491;border-radius:10px;padding:14px 18px}}.card strong{{display:block;font-size:29px}}.RNA_component_supported{{border-color:#d08862;background:#fff5ed}}.metadata_conflict_unresolved{{border-color:#d4b559;background:#fffbee}}.HiC_supported{{border-color:#78a5c0;background:#eef7fd}}.WGS_context_compatible{{border-color:#86aa99;background:#f0f8f3}}section{{background:white;border:1px solid var(--line);border-radius:12px;padding:22px;margin:20px 0}}.columns{{display:grid;grid-template-columns:1.3fr 1fr;gap:22px}}.notice{{border-left:4px solid #c88448;padding:8px 15px;background:#fff9f2}}.toolbar{{display:flex;flex-wrap:wrap;gap:12px;align-items:center;margin:18px 0}}input,select{{font:inherit;border:1px solid #bac9d4;border-radius:6px;padding:9px 12px}}input{{min-width:270px;flex:1}}.tablewrap{{overflow-x:auto}}table{{width:100%;border-collapse:collapse}}th{{text-align:left;background:#eef3f7;padding:12px;font-size:13px;position:sticky;top:0}}td{{vertical-align:top;padding:14px 12px;border-bottom:1px solid var(--line);min-width:150px}}td:first-child{{min-width:200px}}td:nth-child(2){{min-width:190px;max-width:220px}}td:nth-child(4){{min-width:260px}}td:nth-child(5),td:nth-child(6){{min-width:240px}}small{{display:block;font-size:12px;margin-top:7px}}.motif{{font-family:ui-monospace,monospace;overflow-wrap:anywhere;font-weight:600}}.badge{{display:inline-block;border-width:1px;border-style:solid;border-radius:6px;padding:4px 8px;font-size:13px;font-weight:650}}details{{margin-top:9px;font-size:13px}}summary{{cursor:pointer;color:var(--blue)}}.links{{display:flex;flex-wrap:wrap;gap:8px 20px}}.metrics{{font-size:20px;font-weight:700}}.foot{{font-size:13px;color:var(--muted)}}@media(max-width:850px){{main{{padding:20px 12px}}.cards{{grid-template-columns:repeat(2,1fr)}}.columns{{grid-template-columns:1fr}}h1{{font-size:25px}}section{{padding:16px}}}}@media print{{input,select,.toolbar{{display:none}}main{{max-width:none;padding:0}}section{{break-inside:avoid}}.tablewrap{{overflow:visible}}td{{font-size:10px;min-width:0!important}}}}
</style></head><body><main>
<p class="eyebrow">기존 FASTQ 검증 · 2026-09-22 · TRM 관측과 라이브러리 해석을 함께 보기</p>
<h1>반복서열은 관측됐지만,<br>genomic telomere 해석에는 라이브러리 확인이 필요합니다.</h1>
<p class="lead"><b>JU765는 RNA 성분의 강한 증거</b>가 있습니다. APS26·AF13·APS22는 RNA/DNA 여부가 아직 불분명합니다. 7개 Hi-C 표기는 서열 증거가 뒷받침합니다. 이 검토는 기존 motif와 후보 기록을 보존하면서 해석에 필요한 조건을 추가합니다.</p>
<div class="cards">{cards}</div>
<div class="notice"><b>범위:</b> 기존 13샘플의 {summary['prefix_file_count']}개 파일에서 앞부분 100,000 reads씩, 총 <b>{summary['prefix_read_records']:,} reads / {summary['prefix_bases']/1e6:.0f} Mb</b>를 검사했습니다. 무작위 표본이나 전체 파일 재분석이 아닙니다. 기존 TRM 결과와 motif 필터는 바꾸지 않았습니다.</div>
<div class="columns"><section><h2>JU765: RNA 성분을 지지하는 직접적인 정렬 증거</h2>
<p class="metrics">강한 splice reads {alignment['counts']['strong_spliced_primary_records']:,}건 · exon 비율 {alignment['high_mapq_aligned_reference_ACGT_fractions']['exonic']:.2%}</p>
<p>강한 splice reads는 고품질 primary 정렬 {alignment['counts']['high_mapq_primary_mapped_records']:,}건의 {alignment['high_mapq_strong_spliced_fraction']:.2%}입니다. {alignment['distinct_strong_junctions']:,}개 junction이 {review['strong_junction_scaffolds']:,}개 scaffold에 분포하고, {alignment['distinct_strong_annotation_matched_junctions']:,}개 junction은 annotation과 일치합니다.</p>
<p>이는 RNA를 포함한 transcript-enriched 자료를 강하게 지지합니다. <b>RNA 순도, RNA/DNA 혼합비, 사용한 kit를 산출한 것은 아닙니다.</b> 동일 named strain의 reference/GFF를 사용했지만 reference와 현재 자료의 BioSample은 다릅니다.</p>
<p><a href="JU765_alignment/README.md">정렬 방법·합성 대조군·한계 전체 보기</a> · <a href="JU765_alignment/ERR13623736_prefix/summary.json">정량 결과</a> · <a href="JU765_alignment/ERR13623736_prefix/junctions.tsv">Junction 표</a></p></section>
<section><h2>남은 세 샘플과 Hi-C를 읽는 방법</h2><p><b>APS26·AF13·APS22:</b> WGS/GENOMIC/RANDOM과 RNA PolyA 표기가 충돌합니다. 조사한 catalogue에서 matching reference를 찾지 못해 정렬로 확인하지 못했습니다. 후보를 제거하거나 DNA라고 단정하지 않습니다.</p>
<p><b>Hi-C 7샘플:</b> 내부 GATCGATC가 조성·1차 전이를 반영한 기대치보다 {summary['hic_markov_enrichment_min']:.2f}–{summary['hic_markov_enrichment_max']:.2f}배 많습니다. 일반 genomic 비교 자료는 {summary['genomic_markov_enrichment_min']:.2f}–{summary['genomic_markov_enrichment_max']:.2f}배입니다. Hi-C 정황을 지지하지만 균일한 WGS sampling을 뜻하지 않습니다.</p>
<p class="foot">말단 SL이 없거나 polyA/T·중복이 적다는 이유만으로 DNA라고 판단하지 않습니다. JU765도 정확한 말단 SL은 0건이었습니다. Hi-C 배수는 기술적 비교치이며 p-value가 아닙니다.</p></section></div>
<section><h2>샘플별 판단과 다음 조치</h2><p>“이전 TRM”은 전체 입력에서 이미 얻은 관측입니다. “이번 판단”은 제한된 prefix와 공개 protocol·reference 검토에 따른 라이브러리 해석입니다. 신규 TRM 재분류를 의미하지 않습니다.</p>
<div class="toolbar"><input id="search" type="search" aria-label="샘플 검색" placeholder="Accession, 종 이름, motif 검색"><select id="kind" aria-label="판단별 필터"><option value="">모든 판단</option>{options}</select><span id="count" aria-live="polite">{len(rows)} / {len(rows)} samples</span></div>
<div class="tablewrap"><table id="samples"><thead><tr><th>샘플 / 검사 범위</th><th>이전 TRM 관측</th><th>이번 라이브러리 판단</th><th>핵심 근거</th><th>TRM 해석에 미치는 영향</th><th>다음 단계</th></tr></thead><tbody>{''.join(rendered)}</tbody></table></div></section>
<section><h2>결과 파일과 해석의 경계</h2><div class="links"><a href="samples.tsv">통합 샘플 표 TSV</a><a href="comparison_previous.tsv">이전 결과 대비 해석 변경 TSV</a><a href="README.md">설명 문서</a><a href="sequence_evidence_notes.md">서열 signature 상세</a><a href="sequence_evidence.json">서열 원자료 JSON</a><a href="prefix_inspection/inspection.json">Prefix QC</a><a href="notes.md">Reference 조사</a><a href="reference_inventory.tsv">Reference 목록</a><a href="report_summary.json">재현·해시·범위 JSON</a></div>
<p>RNA 성분 검출은 RNA로만 이루어진 library라는 뜻이 아닙니다. 일반 WGS와 부합하는 정황도 순도나 telomere 검출력을 증명하지 않습니다. 어떤 샘플에서도 이 보고서만으로 TRM 부재·충분한 genomic sampling·ALT를 판정하지 않습니다. SRA 검색 링크는 후보 탐색용이며, WGS 표기 외에 protocol·read layout·BioSample/library의 독립성을 함께 확인해야 합니다.</p>
<p class="foot">원본 FASTQ와 기존 TRM 결과는 그대로 보존했습니다. 이 보고서 생성 과정에서 대체 FASTQ는 다운로드하지 않았습니다. 합성 정렬 대조군은 구현 확인용이며 실제 library의 민감도/특이도 또는 혼합비 교정 자료가 아닙니다.</p></section></main>
<script>const search=document.getElementById('search'),kind=document.getElementById('kind'),rows=[...document.querySelectorAll('#samples tbody tr')];function filter(){{const q=search.value.toLocaleLowerCase(),k=kind.value;let n=0;for(const r of rows){{const show=(!k||r.dataset.kind===k)&&r.textContent.toLocaleLowerCase().includes(q);r.hidden=!show;if(show)n++}}document.getElementById('count').textContent=n+' / '+rows.length+' samples'}}search.addEventListener('input',filter);kind.addEventListener('change',filter);</script></body></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-dir', type=Path, default=ROOT / 'results/library_identity_20260922')
    parser.add_argument('--metadata', type=Path, default=ROOT / 'results/genome_feasibility_20260921/samples.tsv')
    parser.add_argument('--assessments', type=Path, default=ROOT / 'results/wgs_assessed_v040/samples')
    args = parser.parse_args()
    args.report_dir = args.report_dir.resolve()
    rows, signatures, inspection, alignment, review, dependencies = build_rows(args.report_dir, args.metadata.resolve(), args.assessments.resolve())
    write_tsv(args.report_dir / 'samples.tsv', rows)
    comparisons = [{key: row[key] for key in ('sample_id', 'scientific_name', 'previous_TRM_classification',
                    'previous_TRM_label', 'previous_detected_motif', 'library_review_code', 'library_review_label',
                    'key_evidence', 'TRM_interpretation_impact', 'next_step', 'hard_motif_filter_applied',
                    'previous_TRM_outputs_modified', 'previous_assessment_path')} |
                   {'update_scope': 'Library interpretation only; previous motif and TRM classification files retained.'} for row in rows]
    write_tsv(args.report_dir / 'comparison_previous.tsv', comparisons)
    hic = [r['ligation_internal_enrichment_vs_markov1'] for r in rows if r['library_review_code'] == 'HiC_supported']
    genomic = [r['ligation_internal_enrichment_vs_markov1'] for r in rows if r['library_review_code'] == 'WGS_context_compatible']
    summary = {'schema_version': 1, 'sample_count': len(rows),
               'prefix_file_count': len(inspection['results']),
               'prefix_read_records': sum(r['counts']['records'] for r in inspection['results']),
               'prefix_bases': sum(r['counts']['bases'] for r in inspection['results']),
               'review_counts': dict(Counter(r['library_review_code'] for r in rows)),
               'hic_markov_enrichment_min': min(hic), 'hic_markov_enrichment_max': max(hic),
               'genomic_markov_enrichment_min': min(genomic), 'genomic_markov_enrichment_max': max(genomic),
               'scope': 'Interpretation update only; no original FASTQ, TRM measurement or production classification modified.',
               'replacement_fastq_downloaded_by_report': False,
               'hard_motif_filter_applied': False, 'purity_or_mixture_proportions_estimated': False,
               'full_file_rescan_performed_by_report': False,
               'genomic_detection_adequacy_established_by_this_audit': False,
               'script_sha256': sha(Path(__file__)),
               'dependencies': [{'path': str(p), 'sha256': sha(p)} for p in dependencies]}
    (args.report_dir / 'index.html').write_text(render_html(rows, summary, alignment, review))
    (args.report_dir / 'report_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
    table = '\n'.join(f"| {r['sample_id']} | {r['scientific_name']} | {r['previous_detected_motif'] or '미확정'} | {r['library_review_label']} |" for r in rows)
    readme = f'''# 기존 TRM 자료의 라이브러리 정체 검토

[통합 HTML 보고서](index.html) · [샘플별 TSV](samples.tsv) · [이전 결과와 해석 비교](comparison_previous.tsv)

**JU765는 RNA 성분의 강한 증거가 있습니다. APS26 / AF13 / APS22는 RNA/DNA 미확정으로 유지합니다.** Hi-C 표기 7건은 서열의 ligation motif 증가가 뒷받침합니다. 비교용 genomic 2건은 일반 WGS 정황과 부합하지만 순도나 telomere 검출력을 증명한 것은 아닙니다.

기존 13샘플의 **22개 파일에서 처음 100,000 reads씩**, 총 **2,200,000 reads / 327,000,000 bp**를 검사했습니다. 비무작위 prefix이며 mates/PCR copies는 독립 관측이 아닙니다. 기존 전체 입력의 TRM 후보·motif·classification 파일은 그대로 보존했습니다. 새 hard motif filter는 적용하지 않았습니다. 이 보고서와 비교 TSV는 **라이브러리에 따른 해석 조건을 추가**합니다.

| Run | 종 / isolate | 이전 관측 motif | 이번 라이브러리 판단 |
|---|---|---|---|
{table}

## 무엇이 달라졌나요?

- **JU765 / ERR13623736:** 강한 splice read {alignment['counts']['strong_spliced_primary_records']:,}/{alignment['counts']['high_mapq_primary_mapped_records']:,} ({alignment['high_mapq_strong_spliced_fraction']:.2%}), 정렬된 ACGT 염기 중 exon 비율 {alignment['high_mapq_aligned_reference_ACGT_fractions']['exonic']:.2%}, {alignment['distinct_strong_junctions']:,}개 junction / {review['strong_junction_scaffolds']:,}개 scaffold가 RNA 성분을 강하게 지지합니다. RNA 순도·RNA/DNA 혼합비·실제 kit는 미확정입니다. 이 run만으로 genomic WGS가 충분했다거나 TRM 부재·novel telomeric sequence가 있다고 해석하면 안 됩니다.
- **APS26 / AF13 / APS22:** 기존 novel 후보 서열은 남깁니다. 공개 WGS/GENOMIC/RANDOM 대 RNA PolyA 표기가 충돌하고, 조사한 catalogue에서 exact-isolate reference를 찾지 못했습니다. 짧은 말단 SL signature 음성, polyA/T 빈도, 중복률만으로는 RNA/DNA를 판정할 수 없습니다. JU765도 정확한 말단 SL hit가 0건이므로 이 음성 검사는 RNA 배제 근거가 되지 않습니다. 실제 protocol 확인 또는 matching reference·독립 genomic WGS가 필요합니다.
- **Hi-C 7건:** 내부 GATCGATC의 조성/1차 전이 기대치 대비 증가가 {min(hic):.2f}–{max(hic):.2f}×입니다. genomic 비교 자료는 {min(genomic):.2f}–{max(genomic):.2f}×입니다. Hi-C 정황을 지지하지만 Hi-C DNA는 균일한 WGS sampling이 아닙니다. 기존 TTAGGC 검출은 유지하며, 부재/novel 가능성 또는 genome-size/검출력은 별도 일반 WGS 자료로 평가해야 합니다.
- **ERR1233601 / DRR674033:** 기존 TTAGGT / TTAGAC 관측을 유지합니다. 이번 비교만으로 library 순도, genomic coverage 또는 telomere 검출력이 확정되지는 않습니다.

SRA 검색 링크는 [통합 샘플 TSV](samples.tsv)에 포함했습니다. accession 검색 후 WGS 표기만 믿지 말고 construction protocol·분자 source·paired layout·BioSample/library의 독립성을 확인해야 합니다. 이 보고서 생성 과정에서 대체 FASTQ를 다운로드하지 않았습니다.

## 근거와 한계

- [JU765 정렬 전체 보고서](JU765_alignment/README.md), [정량 JSON](JU765_alignment/ERR13623736_prefix/summary.json), [junction 표](JU765_alignment/ERR13623736_prefix/junctions.tsv), [독립 감사](JU765_alignment/validation_audit.json).
- [서열 signature 설명](sequence_evidence_notes.md), [샘플별 수치](sequence_evidence.tsv), [파일별 JSON](sequence_evidence.json), [원래 prefix QC](prefix_inspection/inspection.json).
- [공개 reference·protocol 조사](notes.md), [reference 목록](reference_inventory.tsv), [보고서 입력 SHA-256](report_summary.json).

JU765 reference와 현재 자료는 같은 named strain이지만 다른 BioSample이며 reference가 fragmented 상태입니다. MAPQ ≥30인 primary alignment에서 최소 20-bp intron, canonical splice motif, 양쪽 최소 10-bp exact ACGT anchor를 요구했습니다. 합성 DNA/RNA 대조군은 구현 확인용이고, RNA 대조군은 junction을 포함하도록 의도적으로 구성했으므로 혼합비·민감도/특이도 교정에 사용하지 않았습니다. Read placements나 pair names도 독립 원분자 수가 아닙니다.

Hi-C signature는 가능한 junction 중 GATCGATC 한 종류만 분석했습니다. IID/Markov 배수는 기술적 비교치이며 고차 genomic 반복이나 의존성을 완전히 모델링한 p-value가 아닙니다. 어떤 샘플에서도 이 보고서만으로 TRM 부재·ALT·충분한 genomic 검출력을 판정하지 않습니다.

재생성: `python3 scripts/report_library_identity.py`. 생성기는 원래 JSON/TSV를 읽고 이 보고서·통합표만 씁니다. 서열 audit과 JU765 정렬이 같은 prefix를 사용했는지 SHA/read 수를 확인하며, 이전 TRM JSON의 SHA도 검증합니다.
'''
    (args.report_dir / 'README.md').write_text(readme)
    print(json.dumps({key: summary[key] for key in ('sample_count', 'prefix_file_count', 'prefix_read_records', 'prefix_bases', 'review_counts')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
