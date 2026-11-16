#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Event cards at SUBMODULE level (M00003-1, M00019-2 ...) from a compare8_v4 run.
Events whose two submodules have identical chromosome patterns (the same gene block listed
under different parent modules) are collapsed into one card that lists every equivalent ID.
Usage: python plot_submodule_events.py [prefix] [top_n]   (defaults: mammal8_v4, 80)"""
import csv, html, sys
from collections import defaultdict
from pathlib import Path

prefix = sys.argv[1] if len(sys.argv) > 1 else 'mammal8_v4'
top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 80
CN = {'Balae': '鲸', 'Hippo': '河马', 'cat': '猫', 'meerkat': '狐獴', 'chimpanzee': '黑猩猩',
      'elephant': '大象', 'Dugong': '儒艮', 'inunguis': '亚马逊海牛'}
ORDER = ['Balae', 'Hippo', 'cat', 'meerkat', 'chimpanzee', 'elephant', 'Dugong', 'inunguis']
C1, C2, BAR = '#298cb3', '#e4a347', '#111'
esc = html.escape


def tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


chrs = defaultdict(dict)
dom = defaultdict(dict)
for r in tsv(f'{prefix}_submodule_chromosomes.tsv'):
    chrs[r['Submodule_ID']][r['Species']] = r['Main_chromosome'] if r['Pass_module_support'] == '1' else ''
    dom[r['Submodule_ID']][r['Species']] = r['Dominance']
size = {r['Submodule_ID']: int(r['Gene_count']) for r in tsv(f'{prefix}_submodules.tsv')}
pat = lambda s: tuple(chrs[s].get(sp, '') for sp in ORDER)

groups = {}
for e in tsv(f'{prefix}_submodule_events.tsv'):
    a, b = e['Submodule_1'], e['Submodule_2']
    key = frozenset([pat(a), pat(b)])
    g = groups.setdefault(key, {'score': 0.0, 'ev': e, 'members': []})
    g['members'].append((a, b))
    if float(e['Score']) > g['score']:
        g['score'], g['ev'] = float(e['Score']), e
cards_src = sorted(groups.values(), key=lambda g: -g['score'])


def svg(e):
    a, b = e['Submodule_1'], e['Submodule_2']
    together = set(filter(None, e['Together_species'].split(',')))
    row_h, gap, x0, w1, w2 = 34, 8, 110, 130, 130
    out = [f'<svg viewBox="0 0 {x0 + w1 + w2 + 20} {len(ORDER) * (row_h + gap)}" width="100%" role="img">']
    for i, sp in enumerate(ORDER):
        y = i * (row_h + gap)
        ca, cb = chrs[a].get(sp) or '—', chrs[b].get(sp) or '—'
        out.append(f'<text x="0" y="{y + 23}" class="sp">{esc(CN[sp])}</text>')
        out.append(f'<rect x="{x0}" y="{y}" width="{w1}" height="{row_h}" rx="4" fill="{C1}"/>')
        out.append(f'<rect x="{x0 + w1}" y="{y}" width="{w2}" height="{row_h}" rx="4" fill="{C2}"/>')
        out.append(f'<text x="{x0 + w1 / 2}" y="{y + 23}" text-anchor="middle" class="lab w">{esc(ca)}</text>')
        out.append(f'<text x="{x0 + w1 + w2 / 2}" y="{y + 23}" text-anchor="middle" class="lab k">{esc(cb)}</text>')
        if sp not in together:
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
.eq{font-size:11px;opacity:.65;margin-top:6px;line-height:1.5}
.sp{font-size:15px;fill:currentColor}.lab{font-size:15px;font-weight:600}.lab.w{fill:#fff}.lab.k{fill:#2b1a00}
.legend span{display:inline-block;padding:2px 8px;border-radius:4px;margin-right:8px;font-size:12px}
.note{font-size:12px;opacity:.7}
'''
cards = []
for g in cards_src[:top_n]:
    e = g['ev']
    a, b = e['Submodule_1'], e['Submodule_2']
    others = [f'{x}+{y}' for x, y in g['members'] if (x, y) != (a, b)]
    eq = (f'<div class="eq">等价记录（染色体模式相同，来自不同父模块）：{esc("、".join(others[:8]))}'
          f'{" 等 " + str(len(others)) + " 条" if len(others) > 8 else ""}</div>') if others else ''
    cards.append(
        f'<div class="card"><h3>{esc(a)} + {esc(b)}</h3>'
        f'<div class="meta">{size.get(a, "?")} + {size.get(b, "?")} 基因 · 同染色体 {e["Together_n"]} 种 / 分开 {e["Apart_n"]} 种 · 得分 {e["Score"]}</div>'
        f'{svg(e)}{eq}</div>')
page = f'''<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(prefix)} 子模块事件图</title>
<style>{css}</style></head><body>
<h1>子模块层面的事件图：{len(cards)} / {len(groups)} 张（原始 {sum(len(g["members"]) for g in groups.values())} 条子模块事件合并同模式）</h1>
<p class="legend"><span style="background:{C1};color:#fff">子模块1所在染色体</span><span style="background:{C2}">子模块2所在染色体</span>
<span style="background:{BAR};color:#fff">黑条 = 在该物种分属不同染色体（“—”表示该物种无对应）</span></p>
<p class="note">标题里的 M00003-1 这类编号 = 父模块的某一群特定基因。同一批基因出现在多个父模块下时（如 M00032-1、M00033-1、M00034-1）染色体模式完全相同，合并成一张卡片并在卡片底部列出等价记录。染色体编号在物种间不可比较；无位置信息；事件为重叠的候选。</p>
<div class="grid">{''.join(cards)}</div></body></html>'''
out = Path(f'{prefix}_submodule_events_plot.html')
out.write_text(page, encoding='utf-8')
print('wrote', out, '| unique patterns:', len(groups), '| cards:', len(cards))
