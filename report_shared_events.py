#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML: which derived mutations (breakpoint events) are carried by several species, from compare_orderN.py outputs.
Usage: python report_shared_events.py --out order5/shared_events.html order5/mammal7j order5/mammal8 order5/mammal8bat
"""
import argparse, csv, html
from collections import Counter, defaultdict
from pathlib import Path

CSS = ("body{font-family:-apple-system,'Microsoft YaHei',sans-serif;margin:24px auto;max-width:1200px;padding:0 16px;color:#1f2328;background:#fff;line-height:1.6}"
       "h1{font-size:22px}h2{font-size:19px;margin-top:40px;border-bottom:1px solid #d0d7de;padding-bottom:6px}h3{font-size:15px;margin-top:22px}"
       "table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}th,td{border:1px solid #d0d7de;padding:4px 8px;text-align:left;vertical-align:top}"
       "th{background:#f6f8fa}.num{text-align:right}.note{color:#57606a;font-size:13px}"
       ".c{text-align:center;padding:2px 4px}.clade{background:#1a7f37;color:#fff}.homo{background:#bc4c00;color:#fff}.priv{background:#8250df;color:#fff}"
       ".tag-clade{color:#1a7f37;font-weight:600}.tag-homo{color:#bc4c00;font-weight:600}"
       "@media(prefers-color-scheme:dark){body{background:#0d1117;color:#e6edf3}th{background:#161b22}th,td,h2{border-color:#30363d}.note{color:#8b949e}"
       ".clade{background:#238636}.homo{background:#d9722c}.priv{background:#8957e5}.tag-clade{color:#3fb950}.tag-homo{color:#f0883e}}")


def rd(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('prefixes', nargs='+')
    p.add_argument('--out', default='shared_events.html')
    a = p.parse_args()
    e = html.escape
    h = [f'<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
         f'<title>多物种共有的衍生突变</title><style>{CSS}</style></head><body>',
         '<h1>衍生突变在多个物种中的分布</h1>',
         '<p class="note">一个“突变事件”= 某个区段内、携带同一衍生相邻关系的物种集合（同一区段里携带集合相同的相邻关系变化合并为一个事件；极性不确定的不计）。'
         '“携带”指该物种具有衍生状态（相对祖先状态发生了断点变化）。绿色格 = 携带物种在物种树上构成一个分支（可由一次变化解释的共有衍生）；'
         '橙色格 = 携带物种在树上不相邻（需要多次变化或平行变化）；紫色 = 单个物种私有。基于基因名对应，非严格直系同源。</p>']
    for pre in a.prefixes:
        segs = rd(f'{pre}_segments.tsv')
        names = [k[:-4] for k in segs[0] if k.endswith('_chr')]
        genes = {s['Segment']: int(s['Shared_genes']) for s in segs}
        ev = rd(f'{pre}_events.tsv')
        by = defaultdict(list)
        for r in ev:
            by[(r['Carriers'], r['Kind'])].append(r)
        private = Counter()
        multi = []
        for (car, kind), rs in by.items():
            if kind == 'private':
                private[car] += len(rs)
            else:
                multi.append((car, kind, rs))
        multi.sort(key=lambda x: (-len(x[2]), -len(x[0].split('+')), x[0]))
        nm = sum(len(rs) for _, _, rs in multi)
        ncl = sum(len(rs) for _, k, rs in multi if k == 'shared_clade')
        h.append(f'<h2>{e(Path(pre).name)}（{len(names)} 个物种：{e(", ".join(names))}）</h2>')
        h.append(f'<p>共 {len(ev)} 个突变事件：单物种私有 {sum(private.values())} 个，<b>多物种共有 {nm} 个</b>'
                 f'（其中可由一次变化解释的分支共有 {ncl} 个，树上不相邻的平行变化 {nm - ncl} 个）。</p>')
        sz = Counter(len(c.split('+')) for c, _, rs in multi for _ in rs)
        h.append('<h3>按携带物种数</h3><table><tr><th>携带物种数</th><th class="num">事件数</th></tr>' +
                 ''.join(f'<tr><td>{k}</td><td class="num">{v}</td></tr>' for k, v in sorted(sz.items())) + '</table>')
        h.append('<h3>各物种私有突变数</h3><table><tr>' + ''.join(f'<th>{e(n)}</th>' for n in names) + '</tr><tr>' +
                 ''.join(f'<td class="num">{private[n]}</td>' for n in names) + '</tr></table>')
        h.append('<h3>多物种共有的突变（按事件数排序）</h3><table><tr><th>携带物种</th>' + ''.join(f'<th class="c">{e(n)}</th>' for n in names) +
                 '<th class="num">事件数</th><th class="num">相邻关系变化数</th><th>类型</th><th>区段（基因数）</th></tr>')
        for car, kind, rs in multi:
            cs = set(car.split('+'))
            cls = 'clade' if kind == 'shared_clade' else 'homo'
            cells = ''.join(f'<td class="c {cls if n in cs else ""}">{"●" if n in cs else ""}</td>' for n in names)
            segl = '、'.join(f'{r["Segment"]}({genes[r["Segment"]]})' for r in rs)
            tag = '<span class="tag-clade">分支共有</span>' if kind == 'shared_clade' else '<span class="tag-homo">树上不相邻（平行变化）</span>'
            h.append(f'<tr><td>{e(car)}</td>{cells}<td class="num">{len(rs)}</td><td class="num">{sum(int(r["Junctions"]) for r in rs)}</td><td>{tag}</td><td>{e(segl)}</td></tr>')
        h.append('</table>')
        segm = defaultdict(list)
        for car, kind, rs in multi:
            for r in rs:
                segm[r['Segment']].append((car, kind))
        several = {s: v for s, v in segm.items() if len(v) >= 2}
        rel = Counter()
        rows_html = []
        for s in sorted(several, key=lambda s: -len(several[s])):
            sets = [set(c.split('+')) for c, _ in several[s]]
            kinds = set()
            for i in range(len(sets)):
                for j in range(i + 1, len(sets)):
                    x, y = sets[i], sets[j]
                    kinds.add('disjoint' if not (x & y) else 'nested' if x <= y or y <= x else 'cross')
            for kd in kinds:
                rel[kd] += 1
            lab = '；'.join({'disjoint': '含互不重叠的突变（各自是不同类群的特征）', 'nested': '含嵌套的突变（层级关系）',
                           'cross': '含交叉的突变（相互冲突，无法由同一棵树解释）'}[kd] for kd in ('disjoint', 'nested', 'cross') if kd in kinds)
            chips = ' '.join(f'<span class="tag-{"clade" if k == "shared_clade" else "homo"}">{e(c)}</span>' for c, k in several[s])
            rows_html.append(f'<tr><td>{s}</td><td class="num">{genes[s]}</td><td class="num">{len(several[s])}</td><td>{chips}</td><td>{e(lab)}</td></tr>')
        h.append(f'<h3>同一区段内有多个多物种突变（{len(several)} 个区段）</h3>'
                 f'<p>在 {len(segm)} 个有多物种共有突变的区段里，{len(several)} 个区段同时有 2 个以上不同的多物种突变；'
                 f'其中含互不重叠突变的 {rel["disjoint"]} 个，含嵌套的 {rel["nested"]} 个，含交叉冲突的 {rel["cross"]} 个（一个区段可含多种关系）。'
                 '绿色 = 分支共有，橙色 = 树上不相邻。</p>'
                 '<table><tr><th>区段</th><th class="num">基因</th><th class="num">多物种突变数</th><th>各突变的携带物种</th><th>关系</th></tr>' +
                 ''.join(rows_html) + '</table>')
        pc = Counter()
        for car, _, rs in multi:
            cl = sorted(car.split('+'))
            for i in range(len(cl)):
                for j in range(i + 1, len(cl)):
                    pc[(cl[i], cl[j])] += len(rs)
        h.append('<h3>物种两两共同携带的事件数</h3><table><tr><th></th>' + ''.join(f'<th class="c">{e(n)}</th>' for n in names) + '</tr>')
        mx = max(pc.values()) if pc else 1
        for x in names:
            row = f'<tr><th>{e(x)}</th>'
            for y in names:
                v = pc.get((min(x, y), max(x, y)), 0) if x != y else None
                if v is None:
                    row += '<td class="c" style="background:transparent">—</td>'
                else:
                    alpha = v / mx
                    row += f'<td class="c" style="background:rgba(26,127,55,{alpha:.2f});{"color:#fff;" if alpha > .55 else ""}">{v}</td>'
            h.append(row + '</tr>')
        h.append('</table>')
    h.append('<p class="note">事件数按区段内携带集合合并，不保证彼此独立；两两统计中一个事件会同时计入它所含的所有物种对。</p></body></html>')
    Path(a.out).write_text('\n'.join(h), encoding='utf-8')
    print('written', a.out)


if __name__ == '__main__':
    main()
