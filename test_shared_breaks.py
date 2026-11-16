#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Do species share derived breakpoints more often than chance? Standard library only.

Input: <prefix>_events.tsv from compare_orderN.py (one row = one derived-breakpoint event: the set of species carrying
the derived adjacency state in one segment; polarity-ambiguous junctions are excluded).
Null model: the events x species 0/1 matrix is randomised with fixed row sums (carriers per event) and fixed column sums
(derived events per species) by random 2x2 swaps, so a species with many private changes keeps many but their
co-occurrence with particular partners is shuffled. Statistic per pair = number of events carried by both;
one-sided p = P(null >= observed). Optional contrasts 'A:B,C' count events carried by A and by at least one of B, C.
Usage: python test_shared_breaks.py order5/mammal7j --tree "(elephant,...)" --contrast "Equus:Balae,Hippo" --out order5/mammal7j_test
"""
import argparse, csv, html, random
from itertools import combinations
from pathlib import Path

from compare_orderN import parse_tree, prune, ok_side, tree_str


def main():
    p = argparse.ArgumentParser()
    p.add_argument('prefix')
    p.add_argument('--tree', required=True)
    p.add_argument('--contrast', action='append', default=[], help='"A:B,C" = events with A and at least one of B,C')
    p.add_argument('--n', type=int, default=20000, help='randomisations')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--out', default=None)
    a = p.parse_args()
    out = a.out or f'{a.prefix}_test'
    with open(f'{a.prefix}_segments.tsv', encoding='utf-8-sig', newline='') as f:
        names = [h[:-4] for h in next(csv.reader(f, delimiter='\t')) if h.endswith('_chr')]
    tree = prune(parse_tree(a.tree), set(names))
    ev = []
    with open(f'{a.prefix}_events.tsv', encoding='utf-8-sig', newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            ev.append((r['Segment'], set(r['Carriers'].split('+')), r['Kind']))
    k = len(names)
    idx = {n: i for i, n in enumerate(names)}
    rows = [sum(1 << idx[n] for n in c) for _, c, _ in ev]
    m = len(rows)
    pairs = list(combinations(range(k), 2))
    contrasts = []
    for c in a.contrast:
        left, right = c.split(':')
        contrasts.append((c, idx[left], sum(1 << idx[x] for x in right.split(','))))

    def stats(rs):
        pc = {pr: 0 for pr in pairs}
        for r in rs:
            bits = [i for i in range(k) if r >> i & 1]
            for x, y in combinations(bits, 2):
                pc[(x, y)] += 1
        cc = [sum(1 for r in rs if r >> li & 1 and r & mask) for _, li, mask in contrasts]
        return pc, cc
    obs, obs_c = stats(rows)
    rng = random.Random(a.seed)
    cur = rows[:]
    ge = {pr: 0 for pr in pairs}; sm = {pr: 0 for pr in pairs}
    gec = [0] * len(contrasts); smc = [0] * len(contrasts)
    thin = max(5 * m, 50)
    for _ in range(max(a.n // 10, 2000)):  # burn-in
        pass
    def sweep(rs, steps):
        for _ in range(steps):
            i, j = rng.randrange(m), rng.randrange(m)
            if i == j:
                continue
            x, y = rs[i] & ~rs[j], rs[j] & ~rs[i]
            if x and y:
                bx = 1 << rng.choice([b for b in range(k) if x >> b & 1])
                by = 1 << rng.choice([b for b in range(k) if y >> b & 1])
                rs[i] ^= bx | by; rs[j] ^= bx | by
    sweep(cur, 20 * thin)
    for _ in range(a.n):
        sweep(cur, thin)
        pc, cc = stats(cur)
        for pr in pairs:
            ge[pr] += pc[pr] >= obs[pr]; sm[pr] += pc[pr]
        for i in range(len(contrasts)):
            gec[i] += cc[i] >= obs_c[i]; smc[i] += cc[i]
    res = []
    for pr in pairs:
        pv = (ge[pr] + 1) / (a.n + 1)
        res.append([names[pr[0]], names[pr[1]], obs[pr], sm[pr] / a.n, pv, ok_side(tree, {names[pr[0]], names[pr[1]]})])
    # Benjamini-Hochberg
    order = sorted(range(len(res)), key=lambda i: res[i][4])
    q = [0.0] * len(res); prev = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, res[i][4] * len(res) / rank); q[i] = prev
    for i, r in enumerate(res):
        r.append(q[i])
    with open(f'{out}.tsv', 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Species_1', 'Species_2', 'Observed_events', 'Null_mean', 'P_excess', 'Tree_sisters', 'BH_q'])
        w.writerows([[r[0], r[1], r[2], f'{r[3]:.2f}', f'{r[4]:.4f}', 'yes' if r[5] else 'no', f'{r[6]:.4f}'] for r in res])
    e = html.escape
    css = ("body{font-family:-apple-system,'Microsoft YaHei',sans-serif;margin:24px auto;max-width:1000px;padding:0 16px;color:#1f2328;background:#fff;line-height:1.6}"
           "h1{font-size:22px}h2{font-size:18px;margin-top:30px;border-bottom:1px solid #d0d7de;padding-bottom:6px}"
           "table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}th,td{border:1px solid #d0d7de;padding:4px 8px;text-align:left}"
           "th{background:#f6f8fa}.num{text-align:right}.note{color:#57606a;font-size:13px}.sig{background:#e6f6e6;font-weight:600}.warn{background:#fff1e5;font-weight:600}"
           "@media(prefers-color-scheme:dark){body{background:#0d1117;color:#e6edf3}th{background:#161b22}th,td,h2{border-color:#30363d}.note{color:#8b949e}.sig{background:#14301c}.warn{background:#3a2a14}}")
    h = [f'<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>共有衍生断点显著性检验</title><style>{css}</style></head><body>',
         '<h1>物种两两共有衍生断点的显著性检验</h1>',
         f'<p class="note">数据：{e(a.prefix)}_events.tsv，共 {m} 个衍生断点事件（每个区段内携带物种集合相同的连接合并为一个事件；极性不确定的连接不计）。'
         f'物种树：{e(tree_str(tree))}。零模型：保持每个事件的携带物种数和每个物种的衍生事件总数不变，随机交换 {a.n} 次，'
         '统计两个物种同时携带的事件数，单侧 p = 零模型中不小于观测值的比例；BH q 为多重检验校正。零模型不含物种树，因此树上的姐妹物种预期显著，'
         '真正需要关注的是<b>不是姐妹却显著</b>的物种对（橙色），以及<b>应该显著却不显著</b>的物种对。</p>',
         '<h2>物种两两</h2><table><tr><th>物种1</th><th>物种2</th><th class="num">观测事件数</th><th class="num">零模型均值</th><th class="num">p</th><th class="num">BH q</th><th>树上姐妹</th></tr>']
    for r in sorted(res, key=lambda r: r[4]):
        cls = ('sig' if r[5] else 'warn') if r[4] < 0.05 else ''
        h.append(f'<tr class="{cls}"><td>{e(r[0])}</td><td>{e(r[1])}</td><td class="num">{r[2]}</td><td class="num">{r[3]:.2f}</td>'
                 f'<td class="num">{r[4]:.4f}</td><td class="num">{r[6]:.4f}</td><td>{"是" if r[5] else "否"}</td></tr>')
    h.append('</table>')
    if contrasts:
        h.append('<h2>指定对比</h2><table><tr><th>对比（A:与其中至少一个）</th><th class="num">观测事件数</th><th class="num">零模型均值</th><th class="num">p</th></tr>')
        for i, (lab, _, _) in enumerate(contrasts):
            h.append(f'<tr><td>{e(lab)}</td><td class="num">{obs_c[i]}</td><td class="num">{smc[i] / a.n:.2f}</td><td class="num">{(gec[i] + 1) / (a.n + 1):.4f}</td></tr>')
        h.append('</table>')
    h.append('<p class="note">局限：事件按区段内携带集合合并，仍不能保证彼此独立；基于基因名对应；零模型是一个粗糙的重排参照，不是演化模拟。</p></body></html>')
    Path(f'{out}.html').write_text('\n'.join(h), encoding='utf-8')
    print(f'events={m}  species={k}')
    for r in sorted(res, key=lambda r: r[4])[:12]:
        print(f'  {r[0]:9}{r[1]:9} obs={r[2]:3} null={r[3]:.2f} p={r[4]:.4f} q={r[6]:.4f} sister={r[5]}')
    for i, (lab, _, _) in enumerate(contrasts):
        print(f'  contrast {lab}: obs={obs_c[i]} null={smc[i] / a.n:.2f} p={(gec[i] + 1) / (a.n + 1):.4f}')


if __name__ == '__main__':
    main()
