#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sensitivity of units.py to --min-genes and --min-present (standard library only).
Runs units.py over a grid, writes units/sens/*, sensitivity.tsv and sensitivity.json.
Agreement with the baseline run is the adjusted Rand index over genes placed in both runs.
Usage: python sensitivity.py   (run from the project folder holding the species .txt files)"""
import csv, json, subprocess, sys
from collections import Counter
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPECIES = ['Balae=Balae.txt', 'Hippo=Hippo.txt', 'elephant=elephant.txt', 'cat=cat.txt',
           'meerkat=meerkat.txt', 'chimpanzee=hxx.txt', 'Dugong=Dugong.txt', 'inunguis=inunguis.txt']
BASE = (30, 5)
MIN_GENES = [10, 20, 30, 50, 100]
MIN_PRESENT = [4, 5, 6, 7, 8]
sens = HERE / 'sens'
sens.mkdir(exist_ok=True)


def run(mg, mp):
    out = sens / f'g{mg}_p{mp}'
    cmd = [sys.executable, '-I', str(HERE / 'units.py'), '--keep-loc', '--min-genes', str(mg),
           '--min-present', str(mp), '--out', str(out)]
    for s in SPECIES:
        cmd += ['--species', s]
    txt = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', check=True).stdout
    pool = int(txt.split('genes with >=')[1].split(':')[1].split()[0])
    genes = {r['Gene']: r['Unit_ID'] for r in csv.DictReader(open(f'{out}_genes.tsv', encoding='utf-8-sig'), delimiter='\t')}
    units = list(csv.DictReader(open(f'{out}_units.tsv', encoding='utf-8-sig'), delimiter='\t'))
    return pool, genes, units


def ari(a, b):
    keys = [g for g in a if g in b]
    n = len(keys)
    if n < 2:
        return float('nan')
    cont = Counter((a[g], b[g]) for g in keys)
    ra, rb = Counter(a[g] for g in keys), Counter(b[g] for g in keys)
    s = sum(comb(v, 2) for v in cont.values())
    sa, sb = sum(comb(v, 2) for v in ra.values()), sum(comb(v, 2) for v in rb.values())
    exp = sa * sb / comb(n, 2)
    mx = (sa + sb) / 2
    return (s - exp) / (mx - exp) if mx != exp else 1.0


results = {}
for mg in MIN_GENES:
    for mp in MIN_PRESENT:
        results[(mg, mp)] = run(mg, mp)
base_genes = results[BASE][1]
rows = []
for (mg, mp), (pool, genes, units) in results.items():
    rows.append({'min_genes': mg, 'min_present': mp, 'gene_pool': pool, 'units': len(units),
                 'placed': len(genes), 'coverage': round(len(genes) / pool, 4),
                 'ari_vs_baseline': round(ari(base_genes, genes), 4),
                 'smallest_unit': min(int(u['Genes']) for u in units)})
with open(HERE / 'sensitivity.tsv', 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader(); w.writerows(rows)
json.dump({'base': BASE, 'rows': rows}, open(HERE / 'sensitivity.json', 'w', encoding='utf-8'), ensure_ascii=False)
print('grid', len(rows), 'runs; baseline', BASE)
for r in rows:
    print(r['min_genes'], r['min_present'], 'pool', r['gene_pool'], 'units', r['units'], 'cov', r['coverage'], 'ARI', r['ari_vs_baseline'])
