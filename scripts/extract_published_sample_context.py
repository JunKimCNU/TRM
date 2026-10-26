#!/usr/bin/env python3
"""Extract exact-run context from Lim2023 Table S2 without changing TRM calls.

Uses only the Python standard library. The workbook is a local, archived copy;
this helper performs no download and never joins by species or inferred motif.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile

NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL_ID = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
FIELDS = ['sample_id', 'published_species', 'published_putative_trm',
          'published_classification', 'source', 'supplement_table', 'row_number',
          'same_input_as_publication', 'interpretation']
SOURCE = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10760449/supplementaryFiles'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_ids(path):
    with Path(path).open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    ids = [r['sample_id'] for r in rows]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r'[SED]RR\d+', s) for s in ids):
        raise ValueError(f'{path}: invalid or duplicate run IDs')
    return set(ids)


def worksheet_rows(path, name):
    with ZipFile(path) as archive:
        shared = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            shared = [''.join(e.itertext()) for e in ET.fromstring(
                archive.read('xl/sharedStrings.xml')).findall('s:si', NS)]
        sheets = ET.fromstring(archive.read('xl/workbook.xml')).findall('s:sheets/s:sheet', NS)
        matching = [s for s in sheets if s.attrib['name'] == name]
        if len(matching) != 1:
            raise ValueError(f'Expected one worksheet named {name}')
        links = {e.attrib['Id']: e.attrib['Target'] for e in ET.fromstring(
            archive.read('xl/_rels/workbook.xml.rels'))}
        target = links[matching[0].attrib[REL_ID]]
        internal = target.lstrip('/') if target.startswith('/') else 'xl/' + target
        if '..' in PurePosixPath(internal).parts:
            raise ValueError('Unexpected workbook relationship path')
        result = []
        for row in ET.fromstring(archive.read(internal)).findall('s:sheetData/s:row', NS):
            cells = {}
            for cell in row.findall('s:c', NS):
                value = cell.find('s:v', NS)
                inline = cell.find('s:is', NS)
                text = (value.text or '') if value is not None else (
                    ''.join(inline.itertext()) if inline is not None else '')
                if cell.attrib.get('t') == 's' and text:
                    text = shared[int(text)]
                cells[re.sub(r'\d', '', cell.attrib['r'])] = text
            result.append((int(row.attrib['r']), cells))
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, required=True)
    parser.add_argument('--download', type=Path, help='Original supplementaryFiles ZIP to hash')
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    ids = read_ids(args.manifest)
    if ids != read_ids(args.metadata):
        parser.error('Manifest and metadata run sets differ')
    sheet = worksheet_rows(args.workbook, 'Table S2')
    header = dict(sheet).get(2, {})
    expected = {'B': 'Species', 'E': 'SRA_accession', 'F': 'Putative_TRM', 'P': 'Classification'}
    if any(header.get(col) != value for col, value in expected.items()):
        parser.error('Table S2 column names differ from reviewed workbook')
    output = []
    publication_runs = set()
    for number, cells in sheet:
        if number <= 2 or not cells.get('B'):
            continue
        accessions = set(re.findall(r'[SED]RR\d+', cells.get('E', '')))
        publication_runs.update(accessions)
        for sample in sorted(ids & accessions):
            if sample == 'ERR161911' and cells['B'] == 'Gongylonema pulchrum':
                interpretation = (
                    '같은 ERR161911의 기존 분석은 60 bp로 자른 500만 R1 reads에서 '
                    '23-mer를 세어 5–7 bp 단위를 검색했고 ND(저빈도)로 남겼습니다. '
                    '현재 15 bp 후보는 그 검색 범위 밖이므로 당시 ND는 이 후보의 '
                    '기각·검증이 아닙니다. 같은 입력의 재분석이며 독립 검증으로 세지 않습니다.')
            else:
                interpretation = (
                    f'같은 {sample}에서 기존 방법의 TRM 미확정'
                    '(ND·저빈도)이 보고됐습니다. TRM 부재를 뜻하지 않으며, '
                    '같은 입력의 재분석은 독립 검증이 아닙니다.')
                if cells.get('P') != 'Low_copy_number_motifs':
                    interpretation = (
                        f'같은 {sample}의 기존 보고 분류는 {cells.get("P", "")}입니다. '
                        '논문의 해석을 기록한 것으로 새 검증 결과가 아니며, '
                        '같은 입력의 재분석은 독립 검증이 아닙니다.')
            output.append(dict(zip(FIELDS, [sample, cells['B'], cells.get('F', ''),
                cells.get('P', ''), SOURCE, 'Table S2', number, 'true', interpretation])))
    output.sort(key=lambda r: r['sample_id'])
    if len(output) != len({r['sample_id'] for r in output}):
        parser.error('A current run occurs in more than one publication row')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output)
    inputs = [args.workbook, args.manifest, args.metadata, Path(__file__)]
    if args.download:
        inputs.append(args.download)
    provenance = {'created_utc': datetime.now(timezone.utc).isoformat(),
        'publication': 'https://doi.org/10.1101/gr.278124.123', 'source_url': SOURCE,
        'download_member': 'Supplemental_Tables.xlsx', 'worksheet': 'Table S2',
        'join_rule': 'Exact run-accession token intersection only; no species or motif inference',
        'publication_run_count': len(publication_runs), 'manifest_run_count': len(ids),
        'matched_run_count': len(output), 'matched_run_ids': [r['sample_id'] for r in output],
        'input_sha256': {str(p): sha256(p) for p in inputs},
        'output_sha256': sha256(args.output), 'does_not_modify_thresholds_or_taxon_references': True}
    args.output.with_suffix('.provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps({'matched_run_ids': provenance['matched_run_ids'], 'count': len(output)}))


if __name__ == '__main__':
    main()
