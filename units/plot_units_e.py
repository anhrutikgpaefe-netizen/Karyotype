#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Self-contained HTML: per-species chromosome bars split into conserved units, plus the
sensitivity grid. Standard library only. Usage: python plot_units.py [prefix]  (default mammal8)
Click a unit (bar segment or legend chip) to highlight it in every species."""
import csv, json, sys, html, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
prefix = sys.argv[1] if len(sys.argv) > 1 else 'mammal8e'
CN = {'Balae': '鲸', 'Equus': '马', 'cat': '猫', 'meerkat': '狐獴', 'chimpanzee': '黑猩猩',
      'elephant': '大象', 'Dugong': '儒艮', 'inunguis': '亚马逊海牛'}
ORDER = ['Balae', 'Equus', 'cat', 'meerkat', 'chimpanzee', 'elephant', 'Dugong', 'inunguis']
esc = html.escape


def tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def nk(v):
    return tuple((0, int(x)) if x.isdigit() else (1, x.lower()) for x in re.split(r'(\d+)', str(v)))


units = tsv(HERE / f'{prefix}_units.tsv')
kary = tsv(HERE / f'{prefix}_karyotypes.tsv')
size = {u['Unit_ID']: int(u['Genes']) for u in units}


def color(uid):
    i = int(uid[1:])
    h = (i * 137.508) % 360
    return f'hsl({h:.0f},{58 + (i % 3) * 8}%,{44 + (i % 2) * 10}%)'


maxg = max(int(r['Genes_in_units']) for r in kary)
blocks = []
for sp in ORDER:
    rows = sorted((r for r in kary if r['Species'] == sp), key=lambda r: nk(r['Chromosome']))
    bars = []
    for r in rows:
        parts = [p for p in r['Units(gene_count)'].split(';') if p]
        segs = []
        tot = int(r['Genes_in_units'])
        for p in sorted(parts, key=lambda p: int(p[1:3])):
            uid, n = p[:p.index('(')], int(p[p.index('(') + 1:-1])
            w = n / maxg * 100
            lab = uid if w > 4.2 else ''
            segs.append(f'<span class="seg" data-u="{uid}" style="width:{w:.3f}%;background:{color(uid)}" '
                        f'title="{uid} · {n} 基因（该单元共 {size[uid]}）">{lab}</span>')
        bars.append(f'<div class="row"><span class="chr">{esc(r["Chromosome"])}</span>'
                    f'<div class="bar">{"".join(segs)}</div><span class="tot">{tot}</span></div>')
    blocks.append(f'<section class="sp"><h3>{esc(CN[sp])} <small>{len(rows)} 条染色体</small></h3>{"".join(bars)}</section>')

legend = ''.join(f'<span class="chip" data-u="{u["Unit_ID"]}" style="background:{color(u["Unit_ID"])}" '
                 f'title="{u["Genes"]} 基因">{u["Unit_ID"]}</span>' for u in units)

css = '''
:root{--bg:#fff;--fg:#1c1c1c;--card:#fafafa;--bd:#ddd}
@media(prefers-color-scheme:dark){:root{--bg:#161616;--fg:#e8e8e8;--card:#1f1f1f;--bd:#333}}
body{background:var(--bg);color:var(--fg);font-family:system-ui,"Microsoft YaHei",sans-serif;margin:0 auto;padding:16px;max-width:1100px}
h1{font-size:20px}h2{font-size:17px;margin-top:28px}h3{margin:0 0 6px;font-size:15px}h3 small{font-weight:400;opacity:.6}
.sp{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px 12px;margin:12px 0}
.row{display:flex;align-items:center;gap:6px;margin:2px 0}.chr{width:30px;text-align:right;font-size:12px;opacity:.75}
.bar{display:flex;flex:1;height:18px}.tot{width:36px;font-size:11px;opacity:.55}
.seg{display:inline-block;height:100%;box-sizing:border-box;border-right:1px solid rgba(255,255,255,.55);color:#fff;font-size:10px;line-height:18px;text-align:center;overflow:hidden;cursor:pointer;white-space:nowrap}
.chip{display:inline-block;color:#fff;font-size:11px;padding:1px 6px;border-radius:3px;margin:2px;cursor:pointer}
.dim{opacity:.18}
table{border-collapse:collapse;font-size:13px;margin:8px 0}th,td{border:1px solid var(--bd);padding:4px 10px;text-align:center}
td.base{font-weight:700;outline:2px solid #298cb3}.note{font-size:12px;opacity:.7}
'''
js = '''
let cur=null;
function setU(u){cur=(cur===u)?null:u;document.querySelectorAll('[data-u]').forEach(e=>e.classList.toggle('dim',cur!==null&&e.dataset.u!==cur));}
document.addEventListener('click',e=>{const t=e.target.closest('[data-u]');if(t)setU(t.dataset.u);});
'''
page = f'''<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(prefix)} 保守连锁单元</title>
<style>{css}</style></head><body>
<h1>八物种染色体的保守连锁单元组成（{len(units)} 个单元）</h1>
<p class="note">每条横条是一条染色体，按其包含的单元切段，段宽 ∝ 基因数，颜色在所有物种间对应同一个单元。点击任一段或下方色块，可在所有物种里高亮同一个单元。
无位置信息，段在条内的排列顺序只是按编号排，不代表物理顺序。染色体编号在物种间不可比较。</p>
<div>{legend}</div>
{''.join(blocks)}
<script>{js}</script></body></html>'''
out = HERE / f'{prefix}_units_plot.html'
out.write_text(page, encoding='utf-8')
print('wrote', out)
