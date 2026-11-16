#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Map compare8_v4 submodules onto units (gene overlap) and draw it as one HTML.
Usage (from the project folder): python units/map_v4.py [v4_prefix] [units_prefix]
Outputs units/<units_prefix>_v4_map.tsv and units/<units_prefix>_v4_map.html (standard library only)."""
import csv, html, sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
v4p = sys.argv[1] if len(sys.argv) > 1 else 'mammal8_v4'
up = sys.argv[2] if len(sys.argv) > 2 else 'mammal8'
esc = html.escape


def tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def color(uid):
    i = int(uid[1:])
    return f'hsl({(i * 137.508) % 360:.0f},{58 + (i % 3) * 8}%,{44 + (i % 2) * 10}%)'


sub = defaultdict(set)
parent = {}
for r in tsv(ROOT / f'{v4p}_submodule_genes.tsv'):
    if r['Submodule_ID']:
        sub[r['Submodule_ID']].add(r['Gene_symbol'])
        parent[r['Submodule_ID']] = r['Parent_Module_ID']
unit_genes = defaultdict(set)
gu = {}
for r in tsv(HERE / f'{up}_genes.tsv'):
    unit_genes[r['Unit_ID']].add(r['Gene']); gu[r['Gene']] = r['Unit_ID']
units = tsv(HERE / f'{up}_units.tsv')

rows, by_unit, unmatched = [], defaultdict(list), []
for s, gs in sorted(sub.items()):
    c = Counter(gu[g] for g in gs if g in gu)
    if not c:
        unmatched.append((s, len(gs), 0.0, []))
        rows.append([s, parent[s], len(gs), '', 0, 0, ''])
        continue
    u, n = c.most_common(1)[0]
    frac_sub, frac_unit = n / len(gs), n / len(unit_genes[u])
    row = [s, parent[s], len(gs), u, n, round(frac_sub, 3), round(frac_unit, 3)]
    rows.append(row)
    if frac_sub >= 0.9:
        by_unit[u].append((s, len(gs), frac_sub, frac_unit))
    else:
        unmatched.append((s, len(gs), frac_sub, c.most_common(3)))
with open(HERE / f'{up}_v4_map.tsv', 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Submodule', 'Parent_module', 'Submodule_genes', 'Best_unit', 'Shared_genes', 'Frac_of_submodule', 'Frac_of_unit'])
    w.writerows(rows)

n_match = sum(len(v) for v in by_unit.values())
unit_rows = []
for u in units:
    uid = u['Unit_ID']
    chips = ''.join(
        f'<span class="sub" title="{esc(s)} · {n} 基因 · 子模块内 {fs:.0%} 在 {uid} · {uid} 内 {fu:.0%} 被覆盖">{esc(s)}<i>{n}</i></span>'
        for s, n, fs, fu in sorted(by_unit.get(uid, []), key=lambda x: -x[1]))
    cnt = len(by_unit.get(uid, []))
    flag = ' dup' if cnt >= 3 else (' none' if cnt == 0 else '')
    unit_rows.append(f'<div class="urow{flag}"><span class="u" style="background:{color(uid)}">{uid}</span>'
                     f'<span class="g">{u["Genes"]}</span><div class="subs">{chips or "<em>无对应子模块</em>"}</div>'
                     f'<span class="cnt">{cnt}</span></div>')
um_rows = ''.join(
    f'<tr><td>{esc(s)}</td><td>{n}</td><td>{fs:.0%}</td><td>{esc(", ".join(f"{u}:{k}" for u, k in top) or "无单元基因")}</td></tr>'
    for s, n, fs, top in sorted(unmatched, key=lambda x: -x[1]))
css = '''
:root{--bg:#fff;--fg:#1c1c1c;--card:#fafafa;--bd:#ddd}
@media(prefers-color-scheme:dark){:root{--bg:#161616;--fg:#e8e8e8;--card:#1f1f1f;--bd:#333}}
body{background:var(--bg);color:var(--fg);font-family:system-ui,"Microsoft YaHei",sans-serif;margin:0 auto;padding:16px;max-width:1100px}
h1{font-size:20px}h2{font-size:17px;margin-top:26px}.note{font-size:12px;opacity:.7}
.urow{display:flex;align-items:flex-start;gap:8px;border-bottom:1px solid var(--bd);padding:5px 0}
.u{color:#fff;font-size:12px;padding:2px 7px;border-radius:3px;min-width:30px;text-align:center}
.g{width:44px;font-size:12px;opacity:.7;text-align:right;padding-top:2px}.cnt{width:24px;text-align:right;font-size:12px;opacity:.6;padding-top:2px}
.subs{flex:1}.sub{display:inline-block;font-size:11px;border:1px solid var(--bd);background:var(--card);border-radius:3px;padding:1px 6px;margin:2px 3px 2px 0}
.sub i{font-style:normal;opacity:.55;margin-left:5px}.dup .sub{border-color:#e4a347}.none{background:rgba(200,60,60,.1)}
table{border-collapse:collapse;font-size:13px}th,td{border:1px solid var(--bd);padding:3px 9px;text-align:center}
'''
page = f'''<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>v4 子模块与单元对照</title><style>{css}</style></head><body>
<h1>v4 子模块 ↔ 保守连锁单元</h1>
<p class="note">共 {len(sub)} 个 v4 子模块；其中 {n_match} 个有 ≥90% 的基因落在同一个单元内，对应到 {len(by_unit)} / {len(units)} 个单元；
其余 {len(unmatched)} 个见下表。每行是一个单元（颜色与单元组成图一致），右侧是对应它的 v4 子模块（数字 = 子模块基因数），最右列是子模块个数。
橙色边框表示同一单元被 ≥3 个子模块重复描述；红底表示没有子模块对应。悬停可看重叠比例。</p>
{''.join(unit_rows)}
<h2>没有对应到单一单元的子模块（{len(unmatched)} 个）</h2>
<table><tr><th>子模块</th><th>基因数</th><th>最大单元占比</th><th>主要单元:基因数</th></tr>{um_rows}</table>
</body></html>'''
(HERE / f'{up}_v4_map.html').write_text(page, encoding='utf-8')
print('submodules', len(sub), 'matched>=90%:', n_match, 'units hit:', len(by_unit), 'unmatched:', len(unmatched))
print('units with >=3 submodules:', sum(1 for v in by_unit.values() if len(v) >= 3), '| units with none:', sum(1 for u in units if u['Unit_ID'] not in by_unit))
