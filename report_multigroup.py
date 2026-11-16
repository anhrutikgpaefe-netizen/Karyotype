#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML for segments where >=2 arrangements are each shared by several species (compare_orderN.py output).
Usage: python report_multigroup.py order5/mammal8 --tree "(elephant,(homo,((Hippo,Balae),Equus,((meerkat,cat),Manis))))"
"""
import argparse, csv, html
from pathlib import Path

from compare_orderN import parse_tree, prune, ok_side, tree_str, sankoff

COLORS = ['#ddf4ff', '#fff1e5', '#e6f6e6', '#f3e8ff', '#fff8c5', '#ffe5ea', '#e8ecef']
DARK = ['#12294a', '#3a2a14', '#14301c', '#2d1f4a', '#3b3510', '#43202a', '#262c33']

CSS = ("body{font-family:-apple-system,'Microsoft YaHei',sans-serif;margin:24px auto;max-width:1200px;padding:0 16px;color:#1f2328;background:#fff;line-height:1.6}"
       "h1{font-size:22px}h2{font-size:18px;margin-top:32px;border-bottom:1px solid #d0d7de;padding-bottom:6px}"
       "table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}th,td{border:1px solid #d0d7de;padding:4px 8px;text-align:left;vertical-align:top}"
       "th{background:#f6f8fa}.num{text-align:right}.note{color:#57606a;font-size:13px}"
       ".ok{color:#1a7f37;font-weight:600}.no{color:#bc4c00;font-weight:600}"
       + ''.join(f'.g{i}{{background:{c}}}' for i, c in enumerate(COLORS)) +
       "@media(prefers-color-scheme:dark){body{background:#0d1117;color:#e6edf3}th{background:#161b22}th,td,h2{border-color:#30363d}.note{color:#8b949e}.ok{color:#3fb950}.no{color:#f0883e}"
       + ''.join(f'.g{i}{{background:{c}}}' for i, c in enumerate(DARK)) + '}')


def fmt(s):
    return '+'.join(sorted(s))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('prefix')
    p.add_argument('--tree', required=True)
    p.add_argument('--out', default=None)
    a = p.parse_args()
    out = a.out or f'{a.prefix}_multigroup.html'
    with open(f'{a.prefix}_segments.tsv', encoding='utf-8-sig', newline='') as f:
        rs = list(csv.DictReader(f, delimiter='\t'))
    names = [h[:-4] for h in rs[0] if h.endswith('_chr')]
    tree = prune(parse_tree(a.tree), set(names))
    sel = []
    for r in rs:
        if r['Pattern'] == 'ambiguous':
            continue
        cl = [set(c.split('+')) for c in r['Partition'].split(' | ')]
        if sum(len(c) > 1 for c in cl) >= 2:
            sel.append((r, cl))
    e = html.escape
    det = [x for x in sel if '无法' not in x[0]['Inference']]
    und = [x for x in sel if '无法' in x[0]['Inference']]
    h = [f'<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
         f'<title>两种排列各有多个物种共有的区段</title><style>{CSS}</style></head><body>',
         '<h1>两种排列各自被多个物种共有的区段</h1>',
         f'<p class="note">来源：{e(a.prefix)}_segments.tsv（{len(names)} 个物种：{e(", ".join(names))}）。物种树：{e(tree_str(tree))}（括号内逗号分隔的多个物种为未解析）。'
         f'从全部 {len(rs)} 个区段中筛出“至少有两种排列，每种都被 2 个以上物种共有”的区段，共 <b>{len(sel)}</b> 个，'
         f'其中 {len(det)} 个能判断祖先排列，{len(und)} 个不能。表中同色的染色体格子 = 排列一致的物种；染色体号只是各物种自己的编号，跨物种数字不对应。'
         '基于基因名对应，非严格直系同源。</p>',
         '<h2>区段一览</h2><table><tr><th>区段</th><th class="num">基因</th>' + ''.join(f'<th>{e(n)}</th>' for n in names) +
         '<th>判断</th><th>排列分组与推断</th><th>各组的角色</th></tr>']
    for r, cl in sel:
        cells = ''
        for n in names:
            gi = next(i for i, c in enumerate(cl) if n in c)
            cells += f'<td class="g{gi % len(COLORS)}">{e(r[n + "_chr"])}</td>'
        ok = '无法' not in r['Inference']
        state = {n: i for i, c in enumerate(cl) for n in c}
        _, roots = sankoff(tree, state)
        parts = []
        for i, c in enumerate(cl):
            if i in roots:
                role = '祖先排列（未发生变化的物种）' if len(roots) == 1 else '祖先候选'
            elif len(c) == 1:
                role = '私有衍生（单个物种自己的变化）'
            elif ok_side(tree, c):
                role = '共有衍生（姐妹物种共有一次变化）'
            else:
                role = '衍生，但这些物种在树上不相邻，需多次变化或平行变化'
            parts.append(f'{fmt(c)}：{role}')
        status = '<br>'.join(e(x) for x in parts)
        h.append(f'<tr><td>{r["Segment"]}</td><td class="num">{r["Shared_genes"]}</td>{cells}'
                 f'<td class="{"ok" if ok else "no"}">{"可判断" if ok else "无法唯一确定"}</td><td>{e(r["Partition"])}<br>{e(r["Inference"])}'
                 f'{("<br><b>断点级：</b>" + e(r["Junction_inference"])) if r.get("Junction_inference") else ""}</td><td>{status}</td></tr>')
    h.append('</table>')
    h.append(f'<h2>能判断的区段（{len(det)}）</h2><p>这类区段里，多物种的组中有一个包含外群（elephant 或 homo）和其他物种，它保留祖先排列；'
             '另外的多物种组如果是树上的姐妹物种，就是它们共有的一次衍生变化；单独成一种排列的物种则是各自私有的变化。'
             '祖先排列那一组只是“没有发生变化的物种”，它们在树上不必相邻，所以表中不对它做“是否分支”的检验。</p><ul>')
    for r, cl in det:
        h.append(f'<li>{r["Segment"]}（{r["Shared_genes"]} 基因）：{e(r["Inference"])}</li>')
    h.append(f'</ul><h2>无法判断的区段（{len(und)}）</h2>')
    if und:
        h.append('<p>常见原因有三种：(1) elephant 和 homo 各自成一种排列，祖先候选中总有“某一组 或 elephant 或 homo”，这是树根处的问题，需要更远的外群才能解决；'
                 '(2) 多物种组在树上不是一个分支（例如 Hippo+Manis），要多次变化才能解释，说明有平行变化或判定噪声；(3) 两种排列的差异在树根两侧，三个分支都没法定祖先。</p><ul>')
        for r, cl in und:
            why = []
            single = {n for c in cl if len(c) == 1 for n in c}
            if 'elephant' in single or 'homo' in single:
                why.append(f'{"、".join(n for n in ("elephant", "homo") if n in single)}单独成一种排列')
            bad = [fmt(c) for c in cl if len(c) > 1 and not ok_side(tree, c) and len(c) < len(names) - 2]
            if bad:
                why.append(f'{"、".join(bad)}不是树上的分支')
            if '平行变化' in r['Inference']:
                why.append('含平行变化')
            if '根的两侧' in r['Inference']:
                why.append('差异位于根的两侧')
            h.append(f'<li>{r["Segment"]}（{r["Shared_genes"]} 基因）：{e(r["Partition"])}　→　{e("；".join(why) or "祖先候选不唯一")}</li>')
        h.append('</ul>')
    else:
        h.append('<p>无。</p>')
    h.append('<h2>真正难判的情况</h2><p>如果两种排列各自都有多个物种，而且两边各是一个独立的分支、正好在树根两侧，就和三物种里“姐妹 vs 外群”一样分不清哪边是祖先，'
             '只有更远的外群才能解决。目前树根的一边只有 elephant 一个物种，所以这种情况还没出现；'
             '加入更多非洲兽（Dugong、inunguis）或异关节总目物种后才会出现，报告会标“位于根的两侧”，不会硬给结论。</p></body></html>')
    Path(out).write_text('\n'.join(h), encoding='utf-8')
    print(f'{len(sel)} segments ({len(det)} determinable) ->', out)


if __name__ == '__main__':
    main()
