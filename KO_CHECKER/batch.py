"""Read an experiment list CSV into KO_Checker requests for batch processing."""
from __future__ import annotations

import csv
import io
import re

from pipeline import validate_request

MAX_ROWS = 5000
MAX_ITEMS = 5000
# How rows with several target genes (e.g. double KO "Tet1 / Tet2") are handled.
MULTI_GENE = ('split', 'first', 'skip')
# Accepted column names (lower-case); the first one present is used.
COLUMNS = {
    'wt': ('wt_srx', 'wt'),
    'ko': ('ko_srx', 'ko'),
    'gene': ('target_symbols', 'gene', 'target_label'),
    'experiment_id': ('experiment_id', 'id'),
    'label': ('target_label',),
    'status': ('status',),
    'fastq_gb': ('total_fastq_gb',),
    'sra_gb': ('total_sra_gb',),
    'study': ('study',),
    'description': ('description',),
}
GENE_SEPARATOR = re.compile(r'[\s/,;、]+')


def number(value):
    try:
        return float(str(value).replace(',', ''))
    except ValueError:
        return None


def parse_csv(text):
    """Return rows of the CSV, finding the header line below any title lines.

    Each row has the CSV fields used here, the list of genes, and an error when
    the row cannot become a KO_Checker request.
    """
    if not isinstance(text, str) or not text.strip():
        raise ValueError('CSVが空です。')
    lines = list(csv.reader(io.StringIO(text.lstrip('﻿'))))
    header_at = None
    for line_number, line in enumerate(lines[:50]):
        names = {cell.strip().lower() for cell in line}
        if any(n in names for n in COLUMNS['wt']) and any(n in names for n in COLUMNS['ko']):
            header_at = line_number
            break
    if header_at is None:
        raise ValueError('CSVにWT_SRXとKO_SRXの列が見つかりません（先頭50行以内に列名の行が必要です）。')
    header = [cell.strip().lower() for cell in lines[header_at]]
    columns = {}
    for key, names in COLUMNS.items():
        columns[key] = next((header.index(n) for n in names if n in header), None)
    if columns['gene'] is None:
        raise ValueError('CSVに対象遺伝子の列（target_symbols / gene / target_label）が見つかりません。')
    rows = []
    for line in lines[header_at + 1:]:
        if not any(cell.strip() for cell in line):
            continue
        if len(rows) >= MAX_ROWS:
            raise ValueError(f'CSVの行数は{MAX_ROWS}行までです。')

        def cell(key):
            position = columns[key]
            return line[position].strip() if position is not None and position < len(line) else ''
        row = {key: cell(key) for key in COLUMNS}
        row['index'] = len(rows)
        row['experiment_id'] = row['experiment_id'] or f'行{header_at + 2 + len(rows)}'
        row['genes'] = list(dict.fromkeys(g for g in GENE_SEPARATOR.split(row.pop('gene')) if g))
        row['fastq_gb'], row['sra_gb'] = number(row['fastq_gb']), number(row['sra_gb'])
        row['error'] = ''
        try:
            if not row['genes']:
                raise ValueError('対象遺伝子が空欄です。')
            for gene in row['genes']:
                validate_request({'gene': gene, 'wt': row['wt'], 'ko': row['ko']})
        except ValueError as error:
            row['error'] = str(error)
        rows.append(row)
    if not rows:
        raise ValueError('CSVにデータ行がありません。')
    return rows


def build_items(rows, selected, options):
    """Expand selected CSV rows into one analysis per (experiment, gene)."""
    if not isinstance(selected, list) or not selected or not all(isinstance(i, int) and not isinstance(i, bool) for i in selected):
        raise ValueError('一括処理する行を1つ以上選択してください。')
    by_index = {row['index']: row for row in rows}
    items = []
    for index in dict.fromkeys(selected):
        row = by_index.get(index)
        if row is None:
            raise ValueError('選択した行がCSVにありません。CSVを読み込み直してください。')
        if row['error']:
            raise ValueError(f'{row["experiment_id"]}: {row["error"]}')
        genes = row['genes']
        if len(genes) > 1:
            mode = options.get('multi_gene', 'split')
            genes = [] if mode == 'skip' else genes[:1] if mode == 'first' else genes
        for gene in genes:
            request = validate_request({'gene': gene, 'wt': row['wt'], 'ko': row['ko'],
                                        'assembly': options.get('assembly', 'mm10'), 'threads': options.get('threads'),
                                        'confirm_replicates': options.get('confirm_replicates') is True,
                                        'ram_work': options.get('ram_work') is True})
            items.append({'experiment_id': row['experiment_id'], 'gene': gene, 'request': request,
                          'description': row['description'], 'study': row['study'],
                          'fastq_gb': row['fastq_gb'], 'job_id': None, 'cancelled': False})
    if not items:
        raise ValueError('解析する項目がありません（選択した行がすべて複数遺伝子のためスキップされました）。')
    if len(items) > MAX_ITEMS:
        raise ValueError(f'一括処理は{MAX_ITEMS}解析までです。')
    return items
