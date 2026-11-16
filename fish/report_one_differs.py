#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarise 'one species differs' segments from compare_order3.py outputs into one HTML report.
Usage: python report_one_differs.py --out order3/one_differs_report.html order3/sharks order3/rep_gar ...
(each argument is the --out prefix used with compare_order3.py)
"""
import argparse, csv, html
from collections import Counter, defaultdict
from pathlib import Path


def rd(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def kind(blocks):
    """Inversion-only if the blocks keep their order on the other species' chromosome, else rearrangement."""
    ys = [int(b['Y_start']) for b in blocks]
    mono = ys == sorted(ys) or ys == sorted(ys, reverse=True)
    return '倒位型（块顺序不变，仅方向翻转）' if mono else '重排型（块在染色体上位置改变，含易位/转座）'


def analyse(prefix):
    segs, blocks = rd(f'{prefix}_segments.tsv'), rd(f'{prefix}_blocks.tsv')
    by = defaultdict(list)
    for b in blocks:
        by[(b['Segment'], b['Pair'])].append(b)
    names = [h[:-4] for h in segs[0] if h.endswith('_chr')]
    out = []
    for s in segs:
        if s['Pattern'] != 'one_differs':
            continue
        odd = s['Odd_species']
        pairs = []
        for x in names:
            if x == odd:
                continue
            for key in (f'{x}-{odd}', f'{odd}-{x}'):
                if f'{key}_blocks' in s:
                    bl = by[(s['Segment'], key)]
                    pairs.append((key, int(s[f'{key}_blocks']), int(s[f'{key}_noise_genes']), bl))
        ref = [n for n in names if n != odd]
        out.append(dict(seg=s['Segment'], genes=int(s['Shared_genes']), chrs={n: s[f'{n}_chr'] for n in names},
                        odd=odd, ref=ref, est=float(s[f'Est_private_breaks_{odd}']), pairs=pairs,
                        kind=kind(pairs[0][3]) if pairs[0][1] > 1 else ''))
    tot = Counter(s['Pattern'] for s in segs)
    return names, out, tot, len(segs)


CSS = """body{font-family:-apple-system,'Microsoft YaHei',sans-serif;margin:24px auto;max-width:1100px;padding:0 16px;color:#1f2328;background:#fff}
h1{font-size:22px}h2{font-size:18px;margin-top:36px;border-bottom:1px solid #d0d7de;padding-bottom:6px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}th,td{border:1px solid #d0d7de;padding:4px 8px;text-align:left}
th{background:#f6f8fa}.num{text-align:right}.note{color:#57606a;font-size:13px}
details{margin:2px 0}summary{cursor:pointer;color:#0969da}.chip{display:inline-block;padding:1px 8px;border-radius:10px;background:#ddf4ff;margin-right:4px}
.odd{background:#fff1e5;font-weight:600}@media(prefers-color-scheme:dark){body{background:#0d1117;color:#e6edf3}th{background:#161b22}
th,td{border-color:#30363d}.note{color:#8b949e}.odd{background:#3a2a14}.chip{background:#12294a}h2{border-color:#30363d}}"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument('prefixes', nargs='+')
    p.add_argument('--out', default='one_differs_report.html')
    a = p.parse_args()
    e = html.escape
    h = ['<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
         f'<title>三物种区段顺序：一个不同的情况</title><style>{CSS}</style></head><body>',
         '<h1>三物种对应区段内基因顺序：“一个不同、两个一致”</h1>',
         '<p class="note">区段 = 三个物种共有单拷贝基因按“三者所在染色体组合”划分（≥30 基因）。共线 = 两两比较后基因排列只有一块（整段反向也算共线）。'
         '“不同者”= 与另两者都不共线、而另两者互相共线的物种。基于基因名对应，非严格直系同源；私有断点为近似估计。</p>']
    overall = Counter()
    for pre in a.prefixes:
        names, rows, tot, n = analyse(pre)
        title = ' / '.join(names)
        h.append(f'<h2>{e(title)}</h2>')
        h.append(f'<p>共 {n} 个区段：三者一致 {tot["all_consistent"]}，<b>一个不同 {tot["one_differs"]}</b>，三者互不一致 {tot["all_differ"]}，阈值边缘 {tot["ambiguous"]}。</p>')
        cnt = Counter(r['odd'] for r in rows)
        kc = Counter((r['odd'], r['kind'][:3]) for r in rows)
        h.append('<table><tr><th>不同者</th><th class="num">区段数</th><th class="num">涉及基因</th><th class="num">倒位型</th><th class="num">重排型</th><th class="num">平均私有断点</th></tr>')
        for nm in names:
            rs = [r for r in rows if r['odd'] == nm]
            h.append(f'<tr><td>{e(nm)}</td><td class="num">{len(rs)}</td><td class="num">{sum(r["genes"] for r in rs)}</td>'
                     f'<td class="num">{kc[(nm, "倒位型")]}</td><td class="num">{kc[(nm, "重排型")]}</td>'
                     f'<td class="num">{(sum(r["est"] for r in rs) / len(rs)):.1f}</td></tr>' if rs else
                     f'<tr><td>{e(nm)}</td><td class="num">0</td><td class="num">0</td><td class="num">0</td><td class="num">0</td><td class="num">-</td></tr>')
        h.append('</table>')
        for r in rows:
            overall[(title, r['odd'])] += 1
        if rows:
            h.append('<table><tr><th>区段</th><th class="num">共有基因</th>' + ''.join(f'<th>{e(n_)}染色体</th>' for n_ in names) +
                     '<th>不同者</th><th>类型</th><th class="num">私有断点</th><th>与另两者的块数（噪声基因）</th><th>块明细</th></tr>')
            for r in rows:
                cells = ''.join(f'<td class="{"odd" if n_ == r["odd"] else ""}">{e(r["chrs"][n_])}</td>' for n_ in names)
                pr = '；'.join(f'{e(k)}: {b}块({nz})' for k, b, nz, _ in r['pairs'])
                det = ''
                for k, b, nz, bl in r['pairs']:
                    det += f'<details><summary>{e(k)}（{b}块）</summary><table><tr><th>#</th><th>基因数</th><th>方向</th><th>{e(bl[0]["X_species"])}坐标</th><th>{e(bl[0]["Y_species"])}坐标</th><th>起止基因</th></tr>'
                    for x in bl:
                        det += (f'<tr><td>{x["Block"]}</td><td class="num">{x["Genes"]}</td><td>{x["Orient"]}</td>'
                                f'<td>{x["X_chr"]}:{int(x["X_start"]):,}-{int(x["X_end"]):,}</td><td>{x["Y_chr"]}:{int(x["Y_start"]):,}-{int(x["Y_end"]):,}</td>'
                                f'<td>{e(x["X_first_gene"])}…{e(x["X_last_gene"])}</td></tr>')
                    det += '</table></details>'
                h.append(f'<tr><td>{r["seg"]}</td><td class="num">{r["genes"]}</td>{cells}<td>{e(r["odd"])}</td><td>{e(r["kind"])}</td>'
                         f'<td class="num">{r["est"]:.1f}</td><td>{pr}</td><td>{det}</td></tr>')
            h.append('</table>')
    h.append('<p class="note">报告由 report_one_differs.py 根据 compare_order3.py 的 *_segments.tsv / *_blocks.tsv 生成。</p></body></html>')
    Path(a.out).write_text('\n'.join(h), encoding='utf-8')
    print('written', a.out)


if __name__ == '__main__':
    main()
