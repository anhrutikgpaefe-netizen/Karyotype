#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""For 'one species differs' segments (from compare_order3.py): are the two agreeing species close or distant
relatives, and which arrangement do true outgroups (species outside the triple's clade) share?
Usage: python polarize_one_differs.py --dir . --out order3/polarize.html order3/sharks order3/rep_gar ...
Needs <species>_positions.tsv in --dir. Edit TREE if your species relationships differ.
"""
import argparse, csv, html
from collections import Counter, defaultdict
from pathlib import Path
from compare_order3 import read_positions, blocks

# Assumed species tree. Chondrichthyes/Reptilia follow fish8_v2_groups.tsv; the pairing inside the sharks
# (banmasha+jingsha, xiaoyao+zhusha) and bsnake+coriacea is ASSUMED from taxonomy, please check.
TREE = (((('banmasha', 'jingsha'), ('xiaoyao', 'zhusha'))), ('gar', ('Latimeria', ('bsnake', 'coriacea'))))


def leaves(t):
    return {t} if isinstance(t, str) else set().union(*(leaves(x) for x in t))


def lca(t, sp):
    sp = set(sp)
    if not isinstance(t, str):
        for ch in t:
            if sp <= leaves(ch):
                return lca(ch, sp)
    return t


def rd(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def match(genes, x, o, pos, w, min_run, min_genes=20, max_noise=0.3):
    """Is species o collinear with species x over these genes (restricted to o's dominant chromosome)? True/False/None."""
    g = [k for k in genes if k in pos[o]]
    if not g:
        return None
    top, n = Counter(pos[o][k][0] for k in g).most_common(1)[0]
    g = [k for k in g if pos[o][k][0] == top]
    if len(g) < min_genes or len(g) < 0.6 * len(genes):
        return None
    ox = sorted(g, key=lambda k: (pos[x][k][1], pos[x][k][2]))
    oo = {k: i for i, k in enumerate(sorted(g, key=lambda k: (pos[o][k][1], pos[o][k][2])))}
    bl, noise, _, _ = blocks([oo[k] for k in ox], w, min_run)
    if noise > max_noise * len(g):
        return None
    return len(bl) == 1


def junctions(genes, x, y, pos, w, min_run):
    """Gene pairs adjacent in species x but not in y (block boundaries of the x->y comparison)."""
    ox = sorted(genes, key=lambda k: (pos[x][k][1], pos[x][k][2]))
    oy = {k: i for i, k in enumerate(sorted(genes, key=lambda k: (pos[y][k][1], pos[y][k][2])))}
    bl, _, _, _ = blocks([oy[k] for k in ox], w, min_run)
    return [(ox[a_[-1]], ox[b_[0]]) for a_, b_ in zip(bl, bl[1:])]


def outgroup_hits(js, genes, o, pos, w, min_genes=20):
    """(adjacent in outgroup o, testable) among junction pairs, on o's dominant chromosome for these genes."""
    g = [k for k in genes if k in pos[o]]
    if not g:
        return 0, 0
    top = Counter(pos[o][k][0] for k in g).most_common(1)[0][0]
    g = [k for k in g if pos[o][k][0] == top]
    if len(g) < min_genes or len(g) < 0.6 * len(genes):
        return 0, 0
    rk = {k: i for i, k in enumerate(sorted(g, key=lambda k: (pos[o][k][1], pos[o][k][2])))}
    t = [(a_, b_) for a_, b_ in js if a_ in rk and b_ in rk]
    return sum(abs(rk[a_] - rk[b_]) <= w for a_, b_ in t), len(t)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('prefixes', nargs='+')
    p.add_argument('--dir', default='.')
    p.add_argument('--out', default='polarize.html')
    p.add_argument('--window', type=int, default=3)
    p.add_argument('--min-run', type=int, default=5)
    a = p.parse_args()
    allsp = sorted(leaves(TREE))
    pos = {s: read_positions(Path(a.dir) / f'{s}_positions.tsv', False) for s in allsp}
    rows = []
    for pre in a.prefixes:
        segs = rd(f'{pre}_segments.tsv')
        names = [h[:-4] for h in segs[0] if h.endswith('_chr')]
        shared = set.intersection(*(set(pos[n]) for n in names))
        grp = defaultdict(list)
        for g in shared:
            grp[tuple(pos[n][g][0] for n in names)].append(g)
        clade = leaves(lca(TREE, names))
        outs = [s for s in allsp if s not in clade]
        for s in segs:
            if s['Pattern'] != 'one_differs':
                continue
            odd = s['Odd_species']; pair = [n for n in names if n != odd]
            genes = grp[tuple(s[f'{n}_chr'] for n in names)]
            close = odd not in leaves(lca(TREE, pair))
            # junctions private to the odd species (adjacent in D, not in P) vs private to P (adjacent in P, not in D)
            jd = [j for x in pair[:1] for j in junctions(genes, odd, x, pos, a.window, a.min_run)]
            jp = [j for x in pair[:1] for j in junctions(genes, x, odd, pos, a.window, a.min_run)]
            hd = hp = td = tp = 0
            det = []
            for o in outs:
                d1, t1 = outgroup_hits(jd, genes, o, pos, a.window)
                p1, t2 = outgroup_hits(jp, genes, o, pos, a.window)
                hd += d1; td += t1; hp += p1; tp += t2
                det.append(f'{o}: 不同者连接{d1}/{t1}，一致者连接{p1}/{t2}')
            # ancestral = the side whose private junctions the outgroups keep more often (needs >=2 hits and a clear margin)
            fd, fp = (hd / td if td else 0), (hp / tp if tp else 0)
            anc = 'D' if hd >= 2 and fd > 2 * fp else 'P' if hp >= 2 and fp > 2 * fd else ''
            # the triple's own outgroup (first to split off) is the species outside the other two's clade
            outg = next(n for n in names if n not in leaves(lca(TREE, [m for m in names if m != n])))
            if close:  # odd species is the triple's outgroup, the agreeing pair are sisters
                cat = '无法区分：姐妹共有衍生 或 外群私有变化'
                extra = {'D': '；辅助外群检验倾向：姐妹共有衍生', 'P': '；辅助外群检验倾向：外群私有变化', '': ''}[anc]
            else:      # odd species is an ingroup: outgroup + other ingroup agree = ancestral arrangement
                other = next(n for n in pair if n != outg)
                cat = f'祖先排列保留（外群{outg}与{other}一致），{odd}私有衍生'
                extra = '；注意：辅助外群检验有矛盾' if anc == 'D' else ''
            cat += extra
            rows.append(dict(trio='/'.join(names), seg=s['Segment'], genes=len(genes), odd=odd, pair=' + '.join(pair),
                             rel='不同者是外群' if close else '不同者是内群', outs=', '.join(det) or '—', cat=cat, anc=anc,
                             est=s[f'Est_private_breaks_{odd}'], outg=outg))
    with open(Path(a.out).with_suffix('.tsv'), 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['三物种', '区段', '共有基因', '不同者', '一致的两者', '不同者位置', '辅助外群检验', '推断', '不同者私有断点'])
        w.writerows([[r['trio'], r['seg'], r['genes'], r['odd'], r['pair'], r['rel'], r['outs'], r['cat'], r['est']] for r in rows])
    e = html.escape
    css = ("body{font-family:-apple-system,'Microsoft YaHei',sans-serif;margin:24px auto;max-width:1150px;padding:0 16px;color:#1f2328;background:#fff}"
           "h1{font-size:22px}h2{font-size:18px;margin-top:32px;border-bottom:1px solid #d0d7de;padding-bottom:6px}"
           "table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}th,td{border:1px solid #d0d7de;padding:4px 8px;text-align:left;vertical-align:top}"
           "th{background:#f6f8fa}.num{text-align:right}.note{color:#57606a;font-size:13px}"
           "@media(prefers-color-scheme:dark){body{background:#0d1117;color:#e6edf3}th{background:#161b22}th,td,h2{border-color:#30363d}.note{color:#8b949e}}")
    h = [f'<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>一个不同：近缘共有还是远缘祖先</title><style>{css}</style></head><body>',
         '<h1>“一个不同”的区段：姐妹共有突变 还是 祖先排列保留？</h1>',
         '<p class="note">每个三物种按假定物种树确定自己的外群（最早分出者）：((banmasha,jingsha),(xiaoyao,zhusha)) 为鲨类，(gar,(Latimeria,(bsnake,coriacea))) 为硬骨类；鲨类内部配对为按分类学的假定。'
         '不同者是内群之一 → 外群与另一内群一致，该排列为祖先态，不同者私有衍生（确定）；不同者是外群 → 姐妹共有衍生与外群私有变化无法区分。'
         '“辅助外群检验”用三物种之外的其他物种，看双方各自独有的相邻基因对是否仍相邻，功效低，仅作参考。基于基因名对应，非严格直系同源。</p>']
    summ = Counter((r['rel'], r['cat'].split('；')[0] if r['rel'] == '不同者是外群' else '祖先排列保留，不同者私有衍生') for r in rows)
    h.append('<h2>总览</h2><table><tr><th>不同者位置</th><th>推断</th><th class="num">区段数</th></tr>')
    for (rel, cat), n in sorted(summ.items(), key=lambda x: (x[0][0], -x[1])):
        h.append(f'<tr><td>{rel}</td><td>{e(cat)}</td><td class="num">{n}</td></tr>')
    h.append(f'<tr><th colspan="2">合计</th><th class="num">{len(rows)}</th></tr></table>')
    for trio in dict.fromkeys(r['trio'] for r in rows):
        og = next(r['outg'] for r in rows if r['trio'] == trio)
        h.append(f'<h2>{e(trio)}（外群：{e(og)}）</h2><table><tr><th>区段</th><th class="num">基因</th><th>不同者</th><th>一致的两者</th><th>不同者位置</th><th>辅助外群检验</th><th>推断</th><th class="num">私有断点</th></tr>')
        for r in rows:
            if r['trio'] == trio:
                h.append(f'<tr><td>{r["seg"]}</td><td class="num">{r["genes"]}</td><td>{e(r["odd"])}</td><td>{e(r["pair"])}</td><td>{r["rel"]}</td>'
                         f'<td>{e(r["outs"])}</td><td>{e(r["cat"])}</td><td class="num">{r["est"]}</td></tr>')
        h.append('</table>')
    h.append('</body></html>')
    Path(a.out).write_text('\n'.join(h), encoding='utf-8')
    print(f'{len(rows)} segments ->', a.out)
    for (rel, cat), n in sorted(summ.items()):
        print(f'  {rel} | {cat}: {n}')


if __name__ == '__main__':
    main()
