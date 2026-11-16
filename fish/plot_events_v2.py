#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot compare8_v2 events (two modules per species: same chromosome = joined) and a
module-by-species chromosome matrix as one self-contained HTML. Standard library only.
Usage: python plot_events_v2.py [prefix] [top_n]   (defaults: fish8_v2, 60)"""
import csv, json, sys, html
from pathlib import Path

prefix = sys.argv[1] if len(sys.argv) > 1 else 'fish8_v2'
top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 60
CN = {'banmasha': '斑马鲨', 'bsnake': '蛇', 'coriacea': '棱皮龟', 'gar': '雀鳝', 'jingsha': '鲸鲨',
      'Latimeria': '矛尾鱼', 'xiaoyao': '小鳐', 'zhusha': '竹鲨'}
C1, C2, BAR = '#298cb3', '#e4a347', '#111'
esc = html.escape


def tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


ORDER = json.load(open(f'{prefix}_manifest.json', encoding='utf-8'))['species']
mods = {r['Module_ID']: r for r in tsv(f'{prefix}_modules.tsv')}
all_events = tsv(f'{prefix}_events.tsv')
events = all_events[:top_n]
name = lambda s: CN.get(s, s)


def event_svg(ev):
    m1, m2 = mods[ev['Module_1']], mods[ev['Module_2']]
    together = set(filter(None, ev['Together_species'].split(',')))
    row_h, gap, x0, w1, w2 = 34, 8, 110, 130, 130
    H = len(ORDER) * (row_h + gap)
    out = [f'<svg viewBox="0 0 {x0 + w1 + w2 + 20} {H}" width="100%" role="img">']
    for i, sp in enumerate(ORDER):
        y = i * (row_h + gap)
        a, b = m1[f'{sp}_chr'] or '—', m2[f'{sp}_chr'] or '—'
        out.append(f'<text x="0" y="{y + 23}" class="sp">{esc(name(sp))}</text>')
        out.append(f'<rect x="{x0}" y="{y}" width="{w1}" height="{row_h}" rx="4" fill="{C1}"/>')
        out.append(f'<rect x="{x0 + w1}" y="{y}" width="{w2}" height="{row_h}" rx="4" fill="{C2}"/>')
        out.append(f'<text x="{x0 + w1 / 2}" y="{y + 23}" text-anchor="middle" class="lab w">{esc(a)}</text>')
        out.append(f'<text x="{x0 + w1 + w2 / 2}" y="{y + 23}" text-anchor="middle" class="lab k">{esc(b)}</text>')
        if sp not in together:
            out.append(f'<rect x="{x0 + w1 - 4}" y="{y - 1}" width="8" height="{row_h + 2}" rx="1" fill="{BAR}"/>')
    out.append('</svg>')
    return ''.join(out)


def matrix():
    head = ''.join(f'<th>{esc(name(s))}</th>' for s in ORDER)
    rows = []
    for mid, m in mods.items():
        cells = []
        for s in ORDER:
            c, d = m[f'{s}_chr'], m[f'{s}_dominance']
            try:
                dv = float(d)
            except ValueError:
                dv = 0.0
            a = 0.15 + 0.85 * max(0.0, (dv - 0.5) / 0.5)
            cells.append(f'<td style="background:rgba(41,140,179,{a:.2f})" title="主导率 {esc(d)}">{esc(c)}</td>')
        rows.append(f'<tr><th>{mid}</th><td class="anc">{esc(name(m["Anchor_species"]))} {esc(m["Anchor_chr"])}'
                    f'<br><small>{esc(m["Anchor_gene_count"])} 基因</small></td>{"".join(cells)}</tr>')
    return f'<table><tr><th>模块</th><th>锚点</th>{head}</tr>{"".join(rows)}</table>'


css = '''
:root{--bg:#fff;--fg:#1c1c1c;--card:#fafafa;--bd:#ddd}
@media(prefers-color-scheme:dark){:root{--bg:#161616;--fg:#e8e8e8;--card:#1f1f1f;--bd:#333}}
body{background:var(--bg);color:var(--fg);font-family:system-ui,"Microsoft YaHei",sans-serif;margin:0;padding:16px;max-width:1100px;margin:auto}
h1{font-size:20px}h2{font-size:17px;margin-top:28px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:12px}
.card h3{margin:0 0 2px;font-size:15px}.meta{font-size:12px;opacity:.7;margin-bottom:8px}
.sp{font-size:15px;fill:currentColor}.lab{font-size:15px;font-weight:600}.lab.w{fill:#fff}.lab.k{fill:#2b1a00}
.legend span{display:inline-block;padding:2px 8px;border-radius:4px;margin-right:8px;font-size:12px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{border:1px solid var(--bd);padding:4px 6px;text-align:center}td.anc{text-align:left}
td{color:#fff;font-weight:600}td.anc{color:inherit;font-weight:400}
'''

cards = [
    f'<div class="card"><h3>{ev["Event_ID"]}：{ev["Module_1"]} + {ev["Module_2"]}</h3>'
    f'<div class="meta">同染色体 {ev["Together_n"]} 种 / 分开 {ev["Apart_n"]} 种 · 得分 {ev["Score"]}</div>'
    f'{event_svg(ev)}</div>' for ev in events]

page = f'''<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(prefix)} 事件图</title>
<style>{css}</style></head><body>
<h1>八物种染色体联合结构：得分最高的 {len(events)} / {len(all_events)} 条事件</h1>
<p class="legend"><span style="background:{C1};color:#fff">模块1所在染色体</span>
<span style="background:{C2}">模块2所在染色体</span>
<span style="background:{BAR};color:#fff">黑条 = 两模块在该物种分属不同染色体（“—”表示该物种无对应）</span></p>
<p style="font-size:12px;opacity:.7">各物种染色体编号彼此独立，不能跨物种比较号码。事件为重叠的候选，不是独立的历史重排。
这些物种相距很远（软骨鱼、硬骨鱼、爬行类，中间隔着全基因组复制），蛇染色体数少，“蛇单独合并”类事件需谨慎解读。</p>
<div class="grid">{''.join(cards)}</div>
<h2>模块×物种 染色体对应矩阵（{len(mods)} 个模块；颜色深浅 = 主导率，0.5→1.0）</h2>
{matrix()}
</body></html>'''
out = Path(f'{prefix}_events_plot.html')
out.write_text(page, encoding='utf-8')
print('wrote', out, len(events), 'events,', len(mods), 'modules')
