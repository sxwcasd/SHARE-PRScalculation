"""Annotate BCAC 313 PRS using GENCODE v19 (GRCh37).

Requires pandas and numpy. Set source_path, gtf_path, and output_path below
for your environment. Existing output files are protected from overwriting.
Overlapping genes and nearest intergenic genes are recorded separately.
"""

from collections import defaultdict
from pathlib import Path
import csv
import gzip
import hashlib
import io
import os
import re
import tempfile

import numpy as np
import pandas as pd

source_path = Path('/Users/leon/Desktop/Linghao/UofT/Khokha_lab/script/SHARE-PRScalculation/PRSmodels_CanRisk/BCAC_313_PRS.prs')
gtf_path = Path('/Users/leon/Desktop/Linghao/UofT/Khokha_lab/data/reference_files/gencode.v50lift37.annotation.gtf.gz')
output_path = source_path.with_name('BCAC_313_PRS_gene_annotated.prs')
if output_path.exists():
    raise FileExistsError(output_path)

source_bytes = source_path.read_bytes()
source_hash = hashlib.sha256(source_bytes).hexdigest()
header_offset = source_bytes.index(b'Chromosome,')
metadata = source_bytes[:header_offset]
reader = csv.reader(io.StringIO(source_bytes[header_offset:].decode('utf-8')))
original_columns = next(reader)
original_rows = list(reader)
assert len(original_rows) == 313
assert original_columns == [
    'Chromosome', 'Position', 'Reference_Allele', 'Effect_Allele',
    'Effect_Allele_Frequency', 'Log_Odds_Ratio',
]

gene_rows = []
with gzip.open(gtf_path, 'rt') as handle:
    description = next(handle)
    #assert 'GRCh37' in description and 'version 19' in description
    for line in handle:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if fields[2] != 'gene':
            continue
        attributes = dict(re.findall(r'(\S+) "([^"]*)";', fields[8]))
        gene_rows.append((fields[0], int(fields[3]), int(fields[4]), attributes['gene_name']))

genes_df = pd.DataFrame(gene_rows, columns=['chromosome', 'start', 'end', 'gene_name'])
genes_by_chromosome = dict(tuple(genes_df.groupby('chromosome', sort=False)))
annotated_rows = []
overlap_histogram = defaultdict(int)
nearest_ties = []
for row in original_rows:
    chromosome = 'chr' + row[0].removeprefix('chr')
    position = int(row[1])
    genes = genes_by_chromosome[chromosome]
    distances = np.maximum(np.maximum(genes['start'] - position, position - genes['end']), 0)
    overlapping_names = sorted(set(genes.loc[distances.eq(0), 'gene_name']))
    overlap_histogram[len(overlapping_names)] += 1
    nearest_names = []
    nearest_distance = ''
    if not overlapping_names:
        nearest_distance = int(distances.min())
        nearest_names = sorted(set(genes.loc[distances.eq(nearest_distance), 'gene_name']))
        if len(nearest_names) > 1:
            nearest_ties.append((chromosome, position, nearest_names, nearest_distance))
    annotated_rows.append(row + [
        ';'.join(overlapping_names), ';'.join(nearest_names), nearest_distance,
    ])


assert [row[:len(original_columns)] for row in annotated_rows] == original_rows
assert all(row[6] or row[7] for row in annotated_rows)

newline = '\r\n' if b'\r\n' in metadata else '\n'
annotation_metadata = newline.join([
    '#annotation_reference=GENCODE v19 GRCh37.p13 Ensembl 74',
    '#annotation_gtf_sha256=' + hashlib.sha256(gtf_path.read_bytes()).hexdigest(),
    '#annotation_scope=All gene types and the reported one-based variant anchor',
    '#annotation_overlap=Gene-body start and end are inclusive and multiple names use semicolons',
    '#annotation_nearest=Intergenic anchors use nearest gene body on the same chromosome and retain all ties',
    '#annotation_distance_unit=base pairs',
]) + newline
csv_buffer = io.StringIO(newline='')
writer = csv.writer(csv_buffer, lineterminator=newline)
writer.writerow(original_columns + ['gene_name', 'nearest_gene_name', 'nearest_gene_distance_bp'])
writer.writerows(annotated_rows)
output_bytes = annotation_metadata.encode('utf-8') + metadata + csv_buffer.getvalue().encode('utf-8')

# Validate the serialized body before publishing the new file.
new_header_offset = output_bytes.index(b'Chromosome,')
new_reader = csv.reader(io.StringIO(output_bytes[new_header_offset:].decode('utf-8')))
assert next(new_reader) == original_columns + ['gene_name', 'nearest_gene_name', 'nearest_gene_distance_bp']
assert [row[:len(original_columns)] for row in new_reader] == original_rows
assert output_bytes[len(annotation_metadata.encode('utf-8')):new_header_offset] == metadata

with tempfile.NamedTemporaryFile(dir=output_path.parent, prefix='.bcac313-annotation-', delete=False) as handle:
    temporary_path = Path(handle.name)
    handle.write(output_bytes)
try:
    if output_path.exists():
        raise FileExistsError(output_path)
    os.replace(temporary_path, output_path)
finally:
    temporary_path.unlink(missing_ok=True)
assert hashlib.sha256(source_path.read_bytes()).hexdigest() == source_hash
print('Output', output_path)
print('Variant count', len(annotated_rows))
print('Overlapping-gene count histogram', dict(sorted(overlap_histogram.items())))
print('Nearest-gene tie', nearest_ties)
print('Source SHA256', source_hash)
print('Output SHA256', hashlib.sha256(output_bytes).hexdigest())
