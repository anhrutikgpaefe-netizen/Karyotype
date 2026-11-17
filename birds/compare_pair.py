#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-species gene-order comparison inside corresponding chromosome segments. Standard library only.

Shared genes (one unique locus in both species; LOC<digits>, MT, unplaced dropped) are grouped by the chromosome pair they sit on.
A pair with >= --min-genes shared genes is one segment. Inside a segment the genes are ordered by coordinate in both species;
the rank permutation is split into blocks (rank jump > --window starts a new block, runs shorter than --min-run are noise),
the same rule as compare_orderN.py. One block = collinear (a fully reversed segment is still collinear).
With only two species nothing can be said about which one changed.
Usage: python compare_pair.py A=A_positions.tsv B=B_positions.tsv --out birds/Falco_Balearica
"""
import argparse, csv, html, sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from compare_order3 import read_positions, blocks, nk  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument('species', nargs=2, metavar='NAME=POSITIONS_TSV')
    p.add_argument('--min-genes', type=int, default=30)
    p.add_argument('--window', type=int, default=3)
    p.add_argument('--min-run', type=int, default=5)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    (na, fa), (nb, fb) = [s.split('=', 1) for s in a.species]
    pa, pb = read_positions(fa, False), read_positions(fb, False)
    shared = set(pa) & set(pb)
    groups = defaultdict(list)
    for g in shared:
        groups[(pa[g][0], pb[g][0])].append(g)
    segs = sorted((s for s, gl in groups.items() if len(gl) >= a.min_genes), key=lambda s: (-len(groups[s]), s))
    rows, brows = [], []
    for k, (ca, cb) in enumerate(segs, 1):
        gl = groups[(ca, cb)]
        oa = sorted(gl, key=lambda g: (pa[g][1], pa[g][2]))
        rb = {g: i for i, g in enumerate(sorted(gl, key=lambda g: (pb[g][1], pb[g][2])))}
        r = [rb[g] for g in oa]
        bl, noise, d, runs = blocks(r, a.window, a.min_run)
        sid = f'S{k:02}'
        orient = ''.join('+' if len(run) < 2 or d[run[-1]] >= d[run[0]] else '-' for run in runs)
        rows.append([sid, len(gl), ca, cb, len(bl), noise, orient, 'collinear' if len(bl) == 1 else 'rearranged'])
        for b, run in enumerate(runs, 1):
            genes = [oa[i] for i in bl[b - 1]]
            brows.append([sid, b, len(genes), orient[b - 1], ca, pa[genes[0]][1], pa[genes[-1]][2], genes[0], genes[-1],
                          cb, min(pb[g][1] for g in genes), max(pb[g][2] for g in genes)])
    out = Path(a.out)
    with open(f'{out}_segments.tsv', 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Segment', 'Shared_genes', f'{na}_chr', f'{nb}_chr', 'Blocks', 'Noise_genes', 'Block_orientations', 'Status'])
        w.writerows(rows)
    with open(f'{out}_blocks.tsv', 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Segment', 'Block', 'Genes', 'Orient', f'{na}_chr', f'{na}_start', f'{na}_end', f'{na}_first_gene', f'{na}_last_gene',
                    f'{nb}_chr', f'{nb}_start', f'{nb}_end'])
        w.writerows(brows)
    n_col = sum(r[7] == 'collinear' for r in rows)
    print(f'shared genes {len(shared)}; segments(>={a.min_genes}) {len(rows)}, genes in segments {sum(r[1] for r in rows)}; collinear {n_col}, rearranged {len(rows) - n_col}')
    dist = Counter(r[4] for r in rows)
    print('blocks per segment:', dict(sorted(dist.items())))
    # html
    e = html.escape
    css = ("body{font-family:-apple-system,'Microsoft YaHei',sans-serif;margin:24px auto;max-width:1000px;padding:0 16px;color:#1f2328;background:#fff;line-height:1.6}"
           "h1{font-size:22px}h2{font-size:18px;margin-top:28px;border-bottom:1px solid #d0d7de;padding-bottom:6px}"
           "table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}th,td{border:1px solid #d0d7de;padding:4px 8px;text-align:left}"
           "th{background:#f6f8fa}.num{text-align:right}.note{color:#57606a;font-size:13px}.re{background:#fff1e5}.co{background:#e6f6e6}"
           "@media(prefers-color-scheme:dark){body{background:#0d1117;color:#e6edf3}th{background:#161b22}th,td,h2{border-color:#30363d}.note{color:#8b949e}.re{background:#3a2a14}.co{background:#14301c}}")
    h = [f'<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{e(na)} 与 {e(nb)} 的区段内基因顺序</title><style>{css}</style></head><body>',
         f'<h1>{e(na)} 与 {e(nb)}：对应染色体区段内的基因顺序</h1>',
         f'<p class="note">共有单拷贝基因 {len(shared)} 个；按所在染色体对划分，≥{a.min_genes} 个基因的染色体对为一个区段，共 {len(rows)} 个，覆盖 {sum(r[1] for r in rows)} 个基因。'
         f'区段内按坐标排序后检验共线：名次相差超过 {a.window} 开始新的一块，少于 {a.min_run} 个基因的零星错位视为噪声；只剩一块即共线（整段反向也算）。'
         '只有两个物种，无法判断是哪一个发生了变化。基于基因名对应，非严格直系同源。</p>',
         f'<h2>总览</h2><p>共线 {n_col} 个，有重排 {len(rows) - n_col} 个。</p>',
         '<table><tr><th>块数</th><th class="num">区段数</th></tr>' + ''.join(f'<tr><td>{k}</td><td class="num">{v}</td></tr>' for k, v in sorted(dist.items())) + '</table>',
         '<h2>逐区段结果</h2>',
         f'<table><tr><th>区段</th><th class="num">共有基因</th><th>{e(na)}染色体</th><th>{e(nb)}染色体</th><th class="num">块数</th><th class="num">噪声基因</th><th>各块方向</th><th>结论</th></tr>']
    for r in rows:
        cls = 'co' if r[7] == 'collinear' else 're'
        h.append(f'<tr class="{cls}"><td>{r[0]}</td><td class="num">{r[1]}</td><td>{e(r[2])}</td><td>{e(r[3])}</td><td class="num">{r[4]}</td><td class="num">{r[5]}</td><td>{r[6]}</td><td>{"共线" if r[7] == "collinear" else "有重排"}</td></tr>')
    h.append('</table></body></html>')
    Path(f'{out}_report.html').write_text('\n'.join(h), encoding='utf-8')
    print('outputs:', f'{out}_segments.tsv', f'{out}_blocks.tsv', f'{out}_report.html')


if __name__ == '__main__':
    main()
