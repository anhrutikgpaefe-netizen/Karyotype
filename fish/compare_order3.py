#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Three-species gene-order comparison inside corresponding chromosome segments. Standard library only.

1. Genes with a unique, placed locus in all three species (LOC<digits>, MT/Un dropped) are the shared genes.
2. Shared genes are grouped by the triple of chromosomes they sit on; a triple with >= --min-genes genes
   is one corresponding segment (no anchor species, segments never overlap).
3. In each segment, gene order is compared for the three pairs (A-B, A-C, B-C). A pair is collinear when its
   rank permutation forms one block (forward or fully reversed); jumps larger than --window ranks start a new
   block and runs shorter than --min-run genes are treated as noise and removed first.
4. Pattern: all three collinear / one species differs while the other two agree / no pair agrees.
   Private breakpoints per species are estimated from pairwise breakpoints (a=(AB+AC-BC)/2, ...).
Input: *_positions.tsv from parse_ncbi_v2.py (Chromosome, Symbol, Start, End columns).
Usage: python compare_order3.py --species A=A_positions.tsv --species B=B_positions.tsv --species C=C_positions.tsv --out order3_ABC
"""
import argparse, csv, re
from collections import defaultdict
from itertools import combinations
from pathlib import Path

JUNK_CHR = {'MT', 'MITO', 'MITOCHONDRION', 'UN', 'UNPLACED', 'UNLOCALIZED'}
LOC_RE = re.compile(r'^LOC\d+$')


def nk(v):
    return tuple((0, int(x)) if x.isdigit() else (1, x.lower()) for x in re.split(r'(\d+)', str(v)))


def read_positions(path, keep_loc):
    loci = defaultdict(list)
    with open(path, encoding='utf-8-sig', newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            g, c = r['Symbol'].strip().upper(), r['Chromosome'].strip()
            if not g or not c or c.upper() in JUNK_CHR or (not keep_loc and LOC_RE.match(g)):
                continue
            try:
                s, e = int(r['Start']), int(r['End'])
            except (ValueError, TypeError):
                continue
            loci[g].append((c, min(s, e), max(s, e)))
    return {g: v[0] for g, v in loci.items() if len(v) == 1}


def split_runs(d, w):
    runs, cur = [], [0]
    for i in range(1, len(d)):
        if abs(d[i] - d[i - 1]) > w:
            runs.append(cur); cur = []
        cur.append(i)
    runs.append(cur)
    return runs


def dense(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out = [0] * len(vals)
    for k, i in enumerate(order):
        out[i] = k
    return out


def blocks(r, w, min_run):
    """r = rank in species Y of each gene, listed in species X order. Returns (blocks of indices into r, noise count)."""
    keep = list(range(len(r)))
    for _ in range(20):
        d = dense([r[i] for i in keep])
        runs = split_runs(d, w)
        short = [x for x in runs if len(x) < min_run]
        if not short or len(short) == len(runs):
            break
        drop = {i for x in short for i in x}
        keep = [k for j, k in enumerate(keep) if j not in drop]
    d = dense([r[i] for i in keep])
    runs = split_runs(d, w)
    return [[keep[j] for j in x] for x in runs], len(r) - len(keep), d, runs


def orient(d, run):
    return '+' if len(run) < 2 or d[run[-1]] >= d[run[0]] else '-'


def classify(names, cons):
    ok = [p for p in combinations(range(3), 2) if cons[p]]
    if len(ok) == 3:
        return 'all_consistent', '三者顺序一致', ''
    if len(ok) == 1:
        odd = ({0, 1, 2} - set(ok[0])).pop()
        return 'one_differs', f'{names[odd]}与另两者不一致，另两者({names[ok[0][0]]},{names[ok[0][1]]})一致', names[odd]
    if len(ok) == 0:
        return 'all_differ', '三者顺序互不一致', ''
    centre = (set(ok[0]) & set(ok[1])).pop()
    return 'ambiguous', f'阈值边缘：{names[centre]}与其余两者各自一致，但那两者之间不一致', ''


def write(path, head, rows):
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(head); w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--species', action='append', required=True, metavar='NAME=POSITIONS_TSV', help='exactly three')
    p.add_argument('--min-genes', type=int, default=30, help='min shared genes for a segment (default 30)')
    p.add_argument('--window', type=int, default=3, help='max rank jump inside one block (default 3)')
    p.add_argument('--min-run', type=int, default=5, help='runs shorter than this are noise (default 5)')
    p.add_argument('--keep-loc', action='store_true')
    p.add_argument('--out', default='order3')
    a = p.parse_args()
    specs = [tuple(s.split('=', 1)) for s in a.species]
    if len(specs) != 3:
        p.error('exactly three --species required')
    names = [n for n, _ in specs]
    pos = [read_positions(f, a.keep_loc) for _, f in specs]
    shared = set(pos[0]) & set(pos[1]) & set(pos[2])
    groups = defaultdict(list)
    for g in shared:
        groups[tuple(pos[i][g][0] for i in range(3))].append(g)
    segs = sorted((s for s, gl in groups.items() if len(gl) >= a.min_genes), key=lambda s: (-len(groups[s]), s))
    seg_rows, block_rows = [], []
    pairs = list(combinations(range(3), 2))
    counts = defaultdict(int)
    for k, sig in enumerate(segs, 1):
        gl = groups[sig]
        order = [sorted(gl, key=lambda g, i=i: (pos[i][g][1], pos[i][g][2])) for i in range(3)]
        rank = [{g: j for j, g in enumerate(o)} for o in order]
        res = {}
        for x, y in pairs:
            r = [rank[y][g] for g in order[x]]
            bl, noise, d, runs = blocks(r, a.window, a.min_run)
            res[(x, y)] = (bl, noise)
            for b, run in enumerate(runs, 1):
                genes = [order[x][i] for i in bl[b - 1]]
                dd = [r[i] for i in bl[b - 1]]
                block_rows.append([f'S{k:02}', f'{names[x]}-{names[y]}', b, len(genes), orient(d, run),
                                   names[x], sig[x], pos[x][genes[0]][1], pos[x][genes[-1]][2], genes[0], genes[-1],
                                   names[y], sig[y], min(pos[y][order[y][v]][1] for v in dd), max(pos[y][order[y][v]][2] for v in dd),
                                   order[y][dd[0]], order[y][dd[-1]]])
        cons = {pr: len(res[pr][0]) == 1 for pr in pairs}
        bp = {pr: len(res[pr][0]) - 1 for pr in pairs}
        code, text, odd = classify(names, cons)
        est = [(bp[(0, 1)] + bp[(0, 2)] - bp[(1, 2)]) / 2, (bp[(0, 1)] + bp[(1, 2)] - bp[(0, 2)]) / 2,
               (bp[(0, 2)] + bp[(1, 2)] - bp[(0, 1)]) / 2]
        counts[code] += 1
        row = [f'S{k:02}', len(gl)] + list(sig)
        for pr in pairs:
            bl, noise = res[pr]
            o = ''
            if len(bl) == 1:
                r = [rank[pr[1]][g] for g in order[pr[0]]]
                o = '+' if r[bl[0][-1]] >= r[bl[0][0]] else '-'
            row += [len(bl), noise, o]
        row += ['yes' if cons[pr] else 'no' for pr in pairs]
        row += [code, text, odd] + [round(max(e, 0), 1) for e in est]
        seg_rows.append(row)
    head = ['Segment', 'Shared_genes'] + [f'{n}_chr' for n in names]
    for x, y in pairs:
        t = f'{names[x]}-{names[y]}'
        head += [f'{t}_blocks', f'{t}_noise_genes', f'{t}_orient']
    head += [f'{names[x]}-{names[y]}_collinear' for x, y in pairs]
    head += ['Pattern', 'Description', 'Odd_species'] + [f'Est_private_breaks_{n}' for n in names]
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    write(f'{out}_segments.tsv', head, seg_rows)
    write(f'{out}_blocks.tsv', ['Segment', 'Pair', 'Block', 'Genes', 'Orient', 'X_species', 'X_chr', 'X_start', 'X_end',
                                'X_first_gene', 'X_last_gene', 'Y_species', 'Y_chr', 'Y_start', 'Y_end', 'Y_first_gene', 'Y_last_gene'], block_rows)
    tot = sum(counts.values())
    placed = sum(len(groups[s]) for s in segs)
    print(f'species: {", ".join(names)}  shared single-copy genes: {len(shared)}  segments(>={a.min_genes} genes): {tot}  genes in segments: {placed}')
    for code, lab in [('all_consistent', '三者一致'), ('one_differs', '一个不同、两个一致'), ('all_differ', '三者互不一致'), ('ambiguous', '阈值边缘')]:
        print(f'  {lab}: {counts[code]}')
    odd_col = head.index('Odd_species')
    for n in names:
        print(f'    其中 {n} 为不同者: {sum(1 for r in seg_rows if r[odd_col] == n)}')
    print('outputs:', f'{out}_segments.tsv', f'{out}_blocks.tsv')


if __name__ == '__main__':
    main()
