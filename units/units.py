#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Conserved linkage units from joint chromosome signatures. Standard library only.

A gene's signature is the tuple of chromosomes it sits on in every species ('' = absent or
non-unique). Genes sharing a signature form a unit; there is no anchor species, so units
never overlap or duplicate. Outputs <out>_units.tsv, <out>_karyotypes.tsv, <out>_genes.tsv.
Usage: python units.py --species NAME=FILE ... [--min-genes 30] [--min-present 5] [--out mammal8]
"""
import argparse, csv, re
from collections import Counter, defaultdict
from pathlib import Path

JUNK_CHR = {'MT', 'MITO', 'MITOCHONDRION', 'UN', 'UNPLACED', 'UNLOCALIZED'}
LOC_RE = re.compile(r'^LOC\d+$')


def nk(v):
    return tuple((0, int(x)) if x.isdigit() else (1, x.lower()) for x in re.split(r'(\d+)', str(v)))


def read_legacy(path, keep_loc):
    loc = defaultdict(set)
    with open(path, encoding='utf-8-sig', newline='') as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 2:
                continue
            g, c = row[0].strip().upper(), row[1].strip()
            if not g or not c or c.upper() in JUNK_CHR or (not keep_loc and LOC_RE.match(g)):
                continue
            loc[g].add(c)
    return {g: next(iter(v)) for g, v in loc.items() if len(v) == 1}


def compatible(sig, unit):
    return all(a == b for a, b in zip(sig, unit) if a and b)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--species', action='append', required=True, metavar='NAME=FILE')
    p.add_argument('--min-genes', type=int, default=30, help='minimum genes for a unit (default 30)')
    p.add_argument('--min-present', type=int, default=5, help='gene must have a unique chromosome in at least this many species (default 5)')
    p.add_argument('--keep-loc', action='store_true')
    p.add_argument('--out', default='units_out')
    a = p.parse_args()
    specs = [tuple(s.split('=', 1)) for s in a.species]
    names = [n for n, _ in specs]
    loc = {n: read_legacy(f, a.keep_loc) for n, f in specs}
    genes = {g for n in names for g in loc[n]}
    sigs = {}
    for g in genes:
        s = tuple(loc[n].get(g, '') for n in names)
        if sum(1 for x in s if x) >= a.min_present:
            sigs[g] = s
    groups = defaultdict(list)
    for g, s in sigs.items():
        groups[s].append(g)
    # 1) exact signatures with enough genes become units
    core = sorted((s for s, gl in groups.items() if len(gl) >= a.min_genes), key=lambda s: (-len(groups[s]), s))
    members = {s: list(groups[s]) for s in core}
    # 2) remaining genes join a unit when exactly one unit is compatible on the species both have
    rest = [(g, s) for s, gl in groups.items() if len(gl) < a.min_genes for g in gl]
    leftover = []
    for g, s in rest:
        cand = [u for u in core if compatible(s, u)]
        if len(cand) == 1:
            members[cand[0]].append(g)
        else:
            leftover.append(g)
    order = sorted(core, key=lambda s: (-len(members[s]), s))
    uid = {u: f'U{i:02}' for i, u in enumerate(order, 1)}
    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)

    def write(path, head, rows):
        with open(path, 'w', encoding='utf-8-sig', newline='') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(head); w.writerows(rows)

    write(f'{outp}_units.tsv', ['Unit_ID', 'Genes', 'Exact_signature_genes', 'Species_present'] + [f'{n}_chr' for n in names],
          [[uid[u], len(members[u]), len(groups[u]), sum(1 for x in u if x)] + list(u) for u in order])
    unit_of = {g: uid[u] for u in order for g in members[u]}
    write(f'{outp}_genes.tsv', ['Gene', 'Unit_ID'], sorted(unit_of.items()))
    # karyotypes: which units sit on each chromosome of each species
    kt = []
    for i, n in enumerate(names):
        by_chr = defaultdict(Counter)
        for u in order:
            if u[i]:
                by_chr[u[i]][uid[u]] += len(members[u])
        for c in sorted(by_chr, key=nk):
            parts = by_chr[c]
            kt.append([n, c, len(parts), sum(parts.values()),
                       ';'.join(f'{k}({v})' for k, v in sorted(parts.items()))])
    write(f'{outp}_karyotypes.tsv', ['Species', 'Chromosome', 'Units_on_chromosome', 'Genes_in_units', 'Units(gene_count)'], kt)
    placed = sum(len(m) for m in members.values())
    print(f'species: {len(names)}  genes with >={a.min_present} species: {len(sigs)}  distinct signatures: {len(groups)}')
    print(f'units: {len(order)}  genes in units: {placed} ({placed / len(sigs):.1%})  unplaced: {len(leftover)}')
    print('outputs:', f'{outp}_units.tsv', f'{outp}_karyotypes.tsv', f'{outp}_genes.tsv')


if __name__ == '__main__':
    main()
