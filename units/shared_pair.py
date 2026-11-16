#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit pairs where two focal species both differ from a unanimous outgroup (shared fusions / splits).
Usage: python units/shared_pair.py [prefix] [focal1] [focal2]  (default mammal8e Balae Equus)"""
import csv, hashlib, html, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
up = sys.argv[1] if len(sys.argv) > 1 else 'mammal8e'
F1 = sys.argv[2] if len(sys.argv) > 2 else 'Balae'
F2 = sys.argv[3] if len(sys.argv) > 3 else 'Equus'
CN = {'Balae': '鲸', 'Hippo': '河马', 'Equus': '马', 'cat': '猫', 'meerkat': '狐獴', 'chimpanzee': '黑猩猩',
      'elephant': '大象', 'Dugong': '儒艮', 'inunguis': '亚马逊海牛'}
OUT = ['chimpanzee', 'elephant', 'Dugong', 'inunguis']
ORDER = list(dict.fromkeys([F1, F2, 'Balae', 'cat', 'meerkat'] + OUT))
NEED = [F1, F2] + OUT
C1, C2, BAR = '#298cb3', '#e4a347', '#111'
esc = html.escape


def tsv(p):
    with open(p, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


units = tsv(HERE / f'{up}_units.tsv')
sig = {u['Unit_ID']: {s: u[f'{s}_chr'] for s in ORDER} for u in units}
size = {u['Unit_ID']: int(u['Genes']) for u in units}
exact = {u['Unit_ID']: int(u['Exact_signature_genes']) for u in units}
ug = defaultdict(set)
for r in tsv(HERE / f'{up}_genes.tsv'):
    ug[r['Unit_ID']].add(r['Gene'])
cw = {r['Unit_ID']: r for r in tsv(HERE / f'{up}_crosswalk.tsv')}


def old(uid):
    r = cw.get(uid, {})
    return ' ｜ '.join(f'{k}: {r[c].split(";")[0]}' for k, c in
                     (('旧M', 'M_old(mammal8)_containing_>=50%_of_unit'), ('v4', 'v4_containing_>=50%_of_unit'))
                     if r.get(c)) or '无对应'


def ident(uid):
    gs = sorted(ug[uid])
    fp = hashlib.sha1('|'.join(gs).encode('utf-8')).hexdigest()[:6]
    return f'基因指纹 {fp} · {size[uid]} 基因（核心 {exact[uid]}）· 例：{", ".join(gs[:3])}'


ev = []
ids = list(sig)
for i, a in enumerate(ids):
    for b in ids[i + 1:]:
        if not all(sig[a][s] and sig[b][s] for s in NEED):
            continue
        same = {s: (sig[a][s] == sig[b][s]) if sig[a][s] and sig[b][s] else None for s in ORDER}
        o = {same[s] for s in OUT}
        if len(o) != 1:
            continue
        o = o.pop()
        if same[F1] != o and same[F2] != o:
            ev.append(('融合' if same[F1] else '拆分', a, b, same))
ev.sort(key=lambda e: (e[0] != '融合', -min(size[e[1]], size[e[2]])))


def svg(a, b, same):
    rh, gap, x0, w = 34, 8, 110, 130
    out = [f'<svg viewBox="0 0 {x0 + 2 * w + 20} {len(ORDER) * (rh + gap)}" width="100%" role="img">']
    for i, s in enumerate(ORDER):
        y = i * (rh + gap)
        hl = ' font-weight="700"' if s in (F1, F2) else ''
        out.append(f'<text x="0" y="{y + 23}" class="sp"{hl}>{esc(CN[s])}</text>')
        if s in (F1, F2):
            out.append(f'<rect x="{x0 - 6}" y="{y - 3}" width="{2 * w + 12}" height="{rh + 6}" rx="6" fill="none" stroke="#c0392b" stroke-width="2"/>')
        out.append(f'<rect x="{x0}" y="{y}" width="{w}" height="{rh}" rx="4" fill="{C1}"/>')
        out.append(f'<rect x="{x0 + w}" y="{y}" width="{w}" height="{rh}" rx="4" fill="{C2}"/>')
        out.append(f'<text x="{x0 + w / 2}" y="{y + 23}" text-anchor="middle" class="lab w">{esc(sig[a][s])}</text>')
        out.append(f'<text x="{x0 + w * 1.5}" y="{y + 23}" text-anchor="middle" class="lab k">{esc(sig[b][s])}</text>')
        if same[s] is False:
            out.append(f'<rect x="{x0 + w - 4}" y="{y - 1}" width="8" height="{rh + 2}" rx="1" fill="{BAR}"/>')
    return ''.join(out) + '</svg>'


css = '''
:root{--bg:#fff;--fg:#1c1c1c;--card:#fafafa;--bd:#ddd}
@media(prefers-color-scheme:dark){:root{--bg:#161616;--fg:#e8e8e8;--card:#1f1f1f;--bd:#333}}
body{background:var(--bg);color:var(--fg);font-family:system-ui,"Microsoft YaHei",sans-serif;margin:0 auto;padding:16px;max-width:1100px}
h1{font-size:20px}h2{font-size:17px;margin-top:24px}.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:12px}
.card h3{margin:0 0 2px;font-size:15px}.meta{font-size:12px;opacity:.7;margin-bottom:6px}
.xw{font-size:11px;opacity:.65;margin-top:6px;line-height:1.5}
.sp{font-size:15px;fill:currentColor}.lab{font-size:15px;font-weight:600}.lab.w{fill:#fff}.lab.k{fill:#2b1a00}
.legend span{display:inline-block;padding:2px 8px;border-radius:4px;margin-right:8px;font-size:12px}.note{font-size:12px;opacity:.7}
'''
secs = ''
for kind, hint in (('融合', f'外类群（黑猩猩、大象、儒艮、海牛）里两个单元在不同染色体，而{CN[F1]}和{CN[F2]}里在同一条染色体'),
                   ('拆分', f'外类群里两个单元在同一条染色体，而{CN[F1]}和{CN[F2]}里分在不同染色体')):
    cs = [e for e in ev if e[0] == kind]
    cards = ''.join(
        f'<div class="card"><h3>{kind}：{a} + {b}</h3><div class="meta">{size[a]} + {size[b]} 基因</div>{svg(a, b, same)}'
        f'<div class="xw"><b>{a}</b> {esc(ident(a))}<br>旧编号：{esc(old(a))}<br><b>{b}</b> {esc(ident(b))}<br>旧编号：{esc(old(b))}</div></div>'
        for _, a, b, same in cs)
    secs += f'<h2>共有的{kind}（{len(cs)} 个）</h2><p class="note">{hint}。</p><div class="grid">{cards}</div>'
page = f'''<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{CN[F1]}与{CN[F2]}共有的染色体变化</title><style>{css}</style></head><body>
<h1>{CN[F1]}与{CN[F2]}共有的染色体融合与拆分（共 {len(ev)} 个）</h1>
<p class="legend"><span style="background:{C1};color:#fff">单元1所在染色体</span><span style="background:{C2}">单元2所在染色体</span>
<span style="background:{BAR};color:#fff">黑条 = 在该物种分属不同染色体</span><span style="border:2px solid #c0392b">红框 = 比较的两个物种</span></p>
<p class="note">判定条件：4 个外类群在这对单元上状态一致，且两个目标物种都与之不同。注意：只看单元是否同在一条染色体，不含倒位；其余物种不参与判定，仅作参考。
染色体编号在物种间不可比较；“{CN[F1]}和{CN[F2]}都异于外类群”可以是共同祖先继承，也可以是各自独立发生，单靠这两个物种不能区分。</p>
{secs}</body></html>'''
out = HERE / f'{up}_shared_{F1}_{F2}.html'
out.write_text(page, encoding='utf-8')
print('events', len(ev), [(e[0], e[1], e[2]) for e in ev])
