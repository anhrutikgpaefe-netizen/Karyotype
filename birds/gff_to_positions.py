#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NCBI RefSeq GFF3 (.gff.gz) -> <name>_positions.tsv and <name>.txt, the same layout parse_ncbi_v2.py writes. Standard library only.

Only 'gene' and 'pseudogene' features are used (one row per gene, symbol from the gene= attribute). Sequence names come from the
'region' features: chromosome=<n> for chromosomes, MT for the mitochondrion, otherwise the accession (unplaced scaffolds).
Usage: python gff_to_positions.py GFF.gz NAME OUT_DIR
"""
import csv, gzip, re, sys
from collections import defaultdict
from pathlib import Path


def attrs(s):
    d = {}
    for kv in s.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
    return d


def nk(v):
    return tuple((0, int(x)) if x.isdigit() else (1, x.lower()) for x in re.split(r'(\d+)', str(v)))


def main():
    src, name, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    out.mkdir(parents=True, exist_ok=True)
    seqname, build = {}, ''
    genes = []
    skipped = defaultdict(int)
    with gzip.open(src, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                m = re.match(r'#!genome-build-accession\s+NCBI_Assembly:(\S+)', line)
                if m:
                    build = m.group(1)
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            if p[2] == 'region':
                a = attrs(p[8])
                if a.get('genome') == 'mitochondrion':
                    seqname[p[0]] = 'MT'
                elif a.get('chromosome') and a['chromosome'] != 'Unknown':
                    seqname[p[0]] = a['chromosome']
                continue
            if p[2] not in ('gene', 'pseudogene'):
                continue
            a = attrs(p[8])
            sym = a.get('gene') or a.get('Name') or ''
            gid = re.search(r'GeneID:(\d+)', a.get('Dbxref', ''))
            if not sym:
                skipped['no_symbol'] += 1
                continue
            genes.append([seqname.get(p[0], p[0]), sym, int(p[3]), int(p[4]), 'plus' if p[6] == '+' else 'minus',
                          gid.group(1) if gid else '', build, p[0]])
    genes.sort(key=lambda r: (nk(r[0]), r[2], r[3], r[1]))
    counters = defaultdict(int)
    with open(out / f'{name}_positions.tsv', 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Chromosome', 'Symbol', 'Start', 'End', 'Strand', 'GeneID', 'Assembly', 'Gene_order'])
        for r in genes:
            counters[r[0]] += 1
            w.writerow(r[:7] + [counters[r[0]]])
    legacy = defaultdict(set)
    for r in genes:
        legacy[r[0]].add(r[1].upper())
    with open(out / f'{name}.txt', 'w', encoding='utf-8', newline='') as f:
        for chrom in sorted(legacy, key=nk):
            for s in sorted(legacy[chrom]):
                f.write(f'{s}\t{chrom}\n')
    loc = sum(1 for r in genes if re.match(r'^LOC\d+$', r[1]))
    print(f'{name}: {len(genes)} genes on {len(counters)} sequences ({loc} LOC placeholders); skipped {dict(skipped)}')


if __name__ == '__main__':
    main()
