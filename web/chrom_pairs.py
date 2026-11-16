#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pairwise chromosome-correspondence tables ('partner_chr=gene_count') from per-species gene lists.
Standard library only.

Each species file is either the legacy two-column 'SYMBOL<tab>CHROMOSOME' list or a *_positions.tsv
(columns Chromosome, Symbol, ...). A gene is used when its symbol sits on exactly one chromosome in
both species; MT / unplaced are dropped. For a pair (A, B) the count for (chrA, chrB) is the number of
shared symbols located on chrA in A and chrB in B.

Gene-level data of every group of species.tsv, as JSON for build_chrom_page.py (the page itself counts
any pair, and any three species, from it):
    python chrom_pairs.py species.tsv pairs.json
One pair as an Excel-style wide table (row = chromosome of A, cells 'chrB=count', largest first):
    python chrom_pairs.py species.tsv --pair coriacea gar [--min 1] [--out coriacea_gar.tsv]

species.tsv columns: id, group, 中文名, 拉丁名, file (path relative to the project root, i.e. the
parent of this script's directory, or absolute).
"""
import argparse, csv, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
JUNK_CHR = {'MT', 'MITO', 'MITOCHONDRION', 'UN', 'UNPLACED', 'UNLOCALIZED'}
ALPHABET = ''.join(chr(i) for i in range(35, 127) if chr(i) not in '\\<&')   # one character per chromosome index


def nk(v):
    return tuple((0, int(x)) if x.isdigit() else (1, x.lower()) for x in re.split(r'(\d+)', str(v)) if x != '')


def read_genes(path):
    """symbol -> chromosome, for symbols found on exactly one chromosome."""
    loc = defaultdict(set)
    with open(path, encoding='utf-8-sig', newline='') as f:
        rows = csv.reader(f, delimiter='\t')
        head = next(rows, None)
        if head and 'Chromosome' in head and 'Symbol' in head:
            ci, si = head.index('Chromosome'), head.index('Symbol')
        else:
            ci, si = 1, 0
            rows = [head] + list(rows) if head else []
        for row in rows:
            if len(row) <= max(ci, si):
                continue
            g, c = row[si].strip().upper(), row[ci].strip()
            if g and c and c.upper() not in JUNK_CHR:
                loc[g].add(c)
    return {g: next(iter(v)) for g, v in loc.items() if len(v) == 1}


def read_species(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        sp = [r for r in csv.DictReader(f, delimiter='\t') if r.get('id')]
    for s in sp:
        p = Path(s['file'])
        s['path'] = p if p.is_absolute() else ROOT / p
    return sp


def links(ga, gb):
    """[(chrA, chrB, n)] sorted by chrA, then n descending."""
    cnt = Counter((ca, gb[g]) for g, ca in ga.items() if g in gb)
    return sorted(((a, b, n) for (a, b), n in cnt.items()), key=lambda r: (nk(r[0]), -r[2], nk(r[1])))


def wide(rows, minimum):
    by = defaultdict(list)
    for a, b, n in rows:
        if n >= minimum:
            by[a].append(f'{b}={n}')
    return [[a] + cells for a, cells in by.items()]


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('species', help='species.tsv')
    p.add_argument('json', nargs='?', help='output JSON with the gene-level data of every group')
    p.add_argument('--pair', nargs=2, metavar=('A', 'B'), help='print one pair as a wide table instead')
    p.add_argument('--min', type=int, default=1, help='with --pair: drop cells below this gene count')
    p.add_argument('--out', help='with --pair: write the table here (UTF-8 BOM TSV) instead of stdout')
    a = p.parse_args()
    sp = read_species(a.species)
    by_id = {s['id']: s for s in sp}

    if a.pair:
        for x in a.pair:
            if x not in by_id:
                sys.exit(f'unknown species id {x}; known: {", ".join(by_id)}')
        A, B = (by_id[x] for x in a.pair)
        table = [[A['中文名'], B['中文名']]] + wide(links(read_genes(A['path']), read_genes(B['path'])), a.min)
        if a.out:
            with open(a.out, 'w', encoding='utf-8-sig', newline='') as f:
                csv.writer(f, delimiter='\t').writerows(table)
            print('written', a.out, len(table) - 1, 'rows')
        else:
            sys.stdout.reconfigure(encoding='utf-8')
            csv.writer(sys.stdout, delimiter='\t', lineterminator='\n').writerows(table)
        return

    if not a.json:
        p.error('give an output JSON path or --pair A B')
    # Gene-level output: inside a group every species gets one string with one character per gene of the
    # group's gene list (symbols present in >= 2 species); the character is the chromosome's index in that
    # species' label list, '!' = gene absent or not unique. The page counts pairs and triples from these.
    genes = {s['id']: read_genes(s['path']) for s in sp}
    chrs = {i: sorted(set(g.values()), key=nk) for i, g in genes.items()}
    for i, c in chrs.items():
        if len(c) > len(ALPHABET):
            sys.exit(f'{i}: {len(c)} chromosome labels, more than the {len(ALPHABET)} the encoding supports')
    groups, enc = [], {}
    for s in sp:
        g = next((x for x in groups if x['id'] == s['group']), None)
        if g is None:
            g = dict(id=s['group'], label=s['group'], species=[]); groups.append(g)
        g['species'].append(s['id'])
    for g in groups:
        seen = Counter(sym for i in g['species'] for sym in genes[i])
        syms = sorted(sym for sym, n in seen.items() if n >= 2)
        g['genes'] = len(syms)
        for i in g['species']:
            code = {c: ALPHABET[k] for k, c in enumerate(chrs[i])}
            enc[i] = ''.join(code[genes[i][sym]] if sym in genes[i] else '!' for sym in syms)
        print(f'{g["id"]}: {len(g["species"])} species, {len(syms)} gene symbols shared by >= 2 species')
    data = dict(
        species={s['id']: [s['中文名'], s['拉丁名'], len(genes[s['id']]), len(chrs[s['id']])] for s in sp},
        groups=groups, alphabet=ALPHABET, chr=chrs, genes=enc)
    with open(a.json, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, separators=(',', ':'))
    print('written', a.json)


if __name__ == '__main__':
    main()
