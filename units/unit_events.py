#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit-level crosswalk and event cards. Standard library only.
Usage (from the project folder): python units/unit_events.py [units_prefix] [top_n]
Reads units/<prefix>_{units,genes}.tsv plus the old module/submodule files of mammal8, v4 and v5.
Writes units/<prefix>_crosswalk.tsv, units/<prefix>_unit_events.tsv, units/<prefix>_unit_events_plot.html."""
import csv, hashlib, html, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
up = sys.argv[1] if len(sys.argv) > 1 else 'mammal8'
top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 80
CN = {'Balae': '鲸', 'Hippo': '河马', 'cat': '猫', 'meerkat': '狐獴', 'chimpanzee': '黑猩猩',
      'elephant': '大象', 'Dugong': '儒艮', 'inunguis': '亚马逊海牛'}
ORDER = ['Balae', 'Hippo', 'cat', 'meerkat', 'chimpanzee', 'elephant', 'Dugong', 'inunguis']
C1, C2, BAR = '#298cb3', '#e4a347', '#111'
esc = html.escape


def tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


units = tsv(HERE / f'{up}_units.tsv')
sig = {u['Unit_ID']: {sp: u[f'{sp}_chr'] for sp in ORDER} for u in units}
size = {u['Unit_ID']: int(u['Genes']) for u in units}
ugenes = defaultdict(set)
for r in tsv(HERE / f'{up}_genes.tsv'):
    ugenes[r['Unit_ID']].add(r['Gene'])


def gene_sets(rows, key, gene):
    d = defaultdict(set)
    for r in rows:
        if r[key]:
            d[r[key]].add(r[gene])
    return d


sources = {
    'M_old(mammal8)': {r['Module_ID']: set(r['Anchor_gene_symbols'].split('|')) for r in tsv(ROOT / 'mammal8_modules.tsv')},
    'v4': gene_sets(tsv(ROOT / 'mammal8_v4_submodule_genes.tsv'), 'Submodule_ID', 'Gene_symbol'),
    'v5': gene_sets(tsv(ROOT / 'mammal8_v5_submodule_genes.tsv'), 'Submodule_ID', 'Gene_symbol'),
}


def cover(genes, pool, thr=0.5):
    """IDs in pool that contain at least thr of the unit's genes; best first."""
    out = []
    for k, s in pool.items():
        f = len(genes & s) / len(genes)
        if f >= thr:
            out.append((f, k))
    return [f'{k}({f:.0%})' for f, k in sorted(out, reverse=True)]


cross = {}
for u in units:
    uid = u['Unit_ID']
    cross[uid] = {name: cover(ugenes[uid], pool) for name, pool in sources.items()}
with open(HERE / f'{up}_crosswalk.tsv', 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Unit_ID', 'Genes'] + [f'{k}_containing_>=50%_of_unit' for k in sources])
    for u in units:
        w.writerow([u['Unit_ID'], u['Genes']] + [';'.join(cross[u['Unit_ID']][k]) for k in sources])

events = []
ids = [u['Unit_ID'] for u in units]
for i, a in enumerate(ids):
    for b in ids[i + 1:]:
        tog = [sp for sp in ORDER if sig[a][sp] and sig[a][sp] == sig[b][sp]]
        apt = [sp for sp in ORDER if sig[a][sp] and sig[b][sp] and sig[a][sp] != sig[b][sp]]
        if not tog or not apt or len(tog) + len(apt) < 5:
            continue
        score = min(size[a], size[b]) * (len(tog) + len(apt)) / len(ORDER)
        events.append((round(score, 2), a, b, tog, apt))
events.sort(key=lambda e: -e[0])
with open(HERE / f'{up}_unit_events.tsv', 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Event', 'Unit_1', 'Unit_2', 'Together_species', 'Apart_species', 'Together_n', 'Apart_n', 'Score'])
    for k, (s, a, b, tog, apt) in enumerate(events, 1):
        w.writerow([f'UE{k:04}', a, b, ','.join(tog), ','.join(apt), len(tog), len(apt), s])


def label(uid):
    parts = []
    for name, short in (('M_old(mammal8)', '旧M'), ('v4', 'v4'), ('v5', 'v5')):
        c = cross[uid][name][:2]
        if c:
            parts.append(f'{short}: {"、".join(c)}')
    return ' ｜ '.join(parts) or '无对应'


exact = {u['Unit_ID']: int(u['Exact_signature_genes']) for u in units}


def ident(uid):
    """Identity string: the same unit ID always means the same frozen gene set (fingerprint proves it)."""
    gs = sorted(ugenes[uid])
    fp = hashlib.sha1('|'.join(gs).encode('utf-8')).hexdigest()[:6]
    extra = size[uid] - exact[uid]
    core = f'核心 {exact[uid]} + 补入 {extra}' if extra else f'{size[uid]} 基因全部签名一致'
    return f'基因指纹 {fp} · {core} · 例：{", ".join(gs[:3])}'


def svg(a, b, tog):
    row_h, gap, x0, w1, w2 = 34, 8, 110, 130, 130
    out = [f'<svg viewBox="0 0 {x0 + w1 + w2 + 20} {len(ORDER) * (row_h + gap)}" width="100%" role="img">']
    for i, sp in enumerate(ORDER):
        y = i * (row_h + gap)
        ca, cb = sig[a][sp] or '—', sig[b][sp] or '—'
        out.append(f'<text x="0" y="{y + 23}" class="sp">{esc(CN[sp])}</text>')
        out.append(f'<rect x="{x0}" y="{y}" width="{w1}" height="{row_h}" rx="4" fill="{C1}"/>')
        out.append(f'<rect x="{x0 + w1}" y="{y}" width="{w2}" height="{row_h}" rx="4" fill="{C2}"/>')
        out.append(f'<text x="{x0 + w1 / 2}" y="{y + 23}" text-anchor="middle" class="lab w">{esc(ca)}</text>')
        out.append(f'<text x="{x0 + w1 + w2 / 2}" y="{y + 23}" text-anchor="middle" class="lab k">{esc(cb)}</text>')
        if sp not in tog:
            out.append(f'<rect x="{x0 + w1 - 4}" y="{y - 1}" width="8" height="{row_h + 2}" rx="1" fill="{BAR}"/>')
    out.append('</svg>')
    return ''.join(out)


css = '''
:root{--bg:#fff;--fg:#1c1c1c;--card:#fafafa;--bd:#ddd}
@media(prefers-color-scheme:dark){:root{--bg:#161616;--fg:#e8e8e8;--card:#1f1f1f;--bd:#333}}
body{background:var(--bg);color:var(--fg);font-family:system-ui,"Microsoft YaHei",sans-serif;margin:0 auto;padding:16px;max-width:1100px}
h1{font-size:20px}.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:12px}
.card h3{margin:0 0 2px;font-size:15px}.meta{font-size:12px;opacity:.7;margin-bottom:6px}
.xw{font-size:11px;opacity:.65;margin-top:6px;line-height:1.5}
.sp{font-size:15px;fill:currentColor}.lab{font-size:15px;font-weight:600}.lab.w{fill:#fff}.lab.k{fill:#2b1a00}
.legend span{display:inline-block;padding:2px 8px;border-radius:4px;margin-right:8px;font-size:12px}.note{font-size:12px;opacity:.7}
'''
cards = []
for k, (s, a, b, tog, apt) in enumerate(events[:top_n], 1):
    cards.append(
        f'<div class="card"><h3>UE{k:04}：{a} + {b}</h3>'
        f'<div class="meta">{size[a]} + {size[b]} 基因 · 同染色体 {len(tog)} 种 / 分开 {len(apt)} 种 · 得分 {s}</div>'
        f'{svg(a, b, set(tog))}'
        f'<div class="xw"><b>{a}</b> {esc(ident(a))}<br>旧编号：{esc(label(a))}<br><b>{b}</b> {esc(ident(b))}<br>旧编号：{esc(label(b))}</div></div>')
page = f'''<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(up)} 单元事件图</title><style>{css}</style></head><body>
<h1>单元层面的事件图：得分最高的 {len(cards)} / {len(events)} 条</h1>
<p class="legend"><span style="background:{C1};color:#fff">单元1所在染色体</span><span style="background:{C2}">单元2所在染色体</span>
<span style="background:{BAR};color:#fff">黑条 = 在该物种分属不同染色体（“—”表示该物种无对应）</span></p>
<p class="note">每个单元是所有物种签名相同的一批基因，不依赖锚点，不重叠，所以同一块基因只出现一次。卡片底部列出这个单元在旧模块（mammal8）、v4、v5 子模块里对应的编号，括号内是该编号包含单元基因的比例。
事件仍是两两配对：一条含 k 个单元的染色体会产生 k(k-1)/2 条，应结合 mammal8_karyotypes.tsv 判断。染色体编号在物种间不可比较；无位置信息。</p>
<div class="grid">{''.join(cards)}</div></body></html>'''
out = HERE / f'{up}_unit_events_plot.html'
out.write_text(page, encoding='utf-8')
print('units', len(units), '| unit events', len(events), '| cards', len(cards))
for uid in ['U07', 'U22', 'U29', 'U58']:
    print(uid, label(uid))
