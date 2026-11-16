#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a standalone interactive web page from pairs.json (see chrom_pairs.py).
Usage: python build_chrom_page.py pairs.json out.html [--artifact]
"""
import json, sys
from datetime import date

HTML = r'''<!doctype html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>染色体对应关系表</title>
<style>
:root{--bg:#fff;--fg:#1f2328;--muted:#57606a;--line:#d0d7de;--card:#f6f8fa;--accent:#0969da;--heat:9,105,218;--ok:#1a7f37;--warn:#bc4c00}
@media(prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#0d1117;--fg:#e6edf3;--muted:#8b949e;--line:#30363d;--card:#161b22;--accent:#58a6ff;--heat:88,166,255;--ok:#3fb950;--warn:#f0883e;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#0d1117;--fg:#e6edf3;--muted:#8b949e;--line:#30363d;--card:#161b22;--accent:#58a6ff;--heat:88,166,255;--ok:#3fb950;--warn:#f0883e;color-scheme:dark}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.6 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding-block:24px 48px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
.wrap>*{min-width:0}
td.h{font-variant-numeric:tabular-nums}
#copied{color:var(--ok);font-size:13px}
#fallback{width:100%;height:120px;margin-top:8px;font:12px/1.4 ui-monospace,Consolas,monospace;background:var(--bg);color:var(--fg);border:1px solid var(--line);border-radius:7px;padding:8px}
h1{font-size:24px;margin:0 0 4px}.sub{color:var(--muted);margin:0 0 18px;font-size:14px}
.panel{border:1px solid var(--line);border-radius:10px;background:var(--card);padding:12px 14px;margin-bottom:14px}
.controls{display:flex;flex-wrap:wrap;gap:12px 18px;align-items:center}
.controls label{font-size:13px;color:var(--muted);display:flex;align-items:center;gap:8px}
select,input[type=search],button{font:inherit;font-size:14px;color:var(--fg);background:var(--bg);border:1px solid var(--line);border-radius:7px;padding:5px 10px}
button{cursor:pointer}button.on{background:var(--accent);border-color:var(--accent);color:#fff}
.seg{display:inline-flex}.seg button{border-radius:0;margin-left:-1px}.seg button:first-child{border-radius:7px 0 0 7px;margin-left:0}.seg button:last-child{border-radius:0 7px 7px 0}
input[type=range]{width:150px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:10px;margin-bottom:14px}
.card{border:1px solid var(--line);border-radius:10px;padding:10px 14px;background:var(--bg)}
.card b{display:block;font-size:22px;line-height:1.2}.card span{font-size:12.5px;color:var(--muted)}
.scroll{overflow:auto;border:1px solid var(--line);border-radius:10px;max-height:72vh;background:var(--bg)}
table{border-collapse:separate;border-spacing:0;font-size:13px;width:100%}
th,td{padding:4px 8px;border-bottom:1px solid var(--line);border-right:1px solid var(--line);white-space:nowrap}
thead th{position:sticky;top:0;background:var(--card);z-index:2;text-align:center}
tbody th{position:sticky;left:0;background:var(--card);z-index:1;text-align:right}
thead th:first-child{left:0;z-index:3}
td.h{text-align:center;min-width:46px;padding:3px 4px}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
.t th{cursor:pointer}
.t th.sorted::after{content:attr(data-arrow);margin-left:4px;color:var(--accent)}
.chip{display:inline-block;border:1px solid var(--line);border-radius:999px;padding:1px 9px;margin:2px 4px 2px 0;background:var(--card);font-size:12.5px}
.chip i{font-style:normal;color:var(--muted)}
.badge{font-size:12px;border-radius:5px;padding:1px 7px;border:1px solid var(--line)}
.b1{color:var(--ok);border-color:var(--ok)}.b2{color:var(--warn);border-color:var(--warn)}.b3{color:var(--accent);border-color:var(--accent)}.b0{color:var(--muted)}
em.la{color:var(--muted);font-size:.85em}
.wiki{border-radius:6px}
.wiki table{font-size:14px}
.wiki th,.wiki td{padding:5px 10px;vertical-align:middle}
.wiki thead th{text-align:left;border-bottom:2px solid var(--line)}
.wiki tbody th{font-size:15px;min-width:3.2em}
.wiki>table>tbody>tr:nth-child(even)>td{background:var(--card)}
.wiki>table>tbody>tr:hover>td{background:rgba(var(--heat),.07)}
.simple table{width:auto;min-width:320px}
.simple td{font-variant-numeric:tabular-nums}
.lk.s{display:inline;border:0;background:none;padding:1px 3px;margin:0;border-radius:3px}
.lk.s.hl{box-shadow:0 0 0 2px var(--warn);background:rgba(240,136,62,.16)}
.wiki .mini{width:auto;font-size:13px;border-collapse:collapse;margin:2px 0}
.wiki .mini th,.wiki .mini td{position:static;border:1px solid var(--line);padding:3px 10px;text-align:center;min-width:3em;font-size:13px;background:none;font-variant-numeric:tabular-nums}
.wiki .mini th{background:var(--card)}
.mini th small{display:block;font-weight:400;color:var(--muted);font-size:11px;line-height:1.2}
.mini th.corner{font-weight:400;color:var(--muted);font-size:12px}
.wiki .mini td.hom,.sw.hom{background:rgba(26,127,55,.22)}.mini td.hom{font-weight:600}
.wiki .mini td.part,.sw.part{background:rgba(188,76,0,.2)}
.mini td.no{color:var(--muted)}
.sw{display:inline-block;width:14px;height:14px;border:1px solid var(--line);border-radius:3px}
.legend{flex-wrap:wrap}
.wiki tr.target th,.wiki tr.target td{box-shadow:inset 0 2px 0 var(--accent),inset 0 -2px 0 var(--accent)}
.wiki td.links{white-space:normal;width:100%}
a.anc{color:inherit;text-decoration:none;display:block}a.anc:hover{color:var(--accent)}
.lk{display:inline-flex;align-items:center;gap:6px;padding:2px 8px 2px 6px;margin:2px 6px 2px 0;border:1px solid var(--line);border-radius:6px;background:var(--bg);cursor:pointer;font-variant-numeric:tabular-nums;white-space:nowrap}
.lk b{min-width:1.5em;text-align:right}
.lk .tr{width:56px;height:6px;border-radius:3px;background:rgba(var(--heat),.16);overflow:hidden}
.lk .tr i{display:block;height:100%;background:rgb(var(--heat))}
.lk small{color:var(--muted);font-size:12px}
.lk.weak{border-color:transparent;background:none;padding:0 3px;margin:0 2px 0 0;color:var(--muted);font-size:12px}
.lk.hl{border-color:var(--warn);box-shadow:0 0 0 1px var(--warn);background:rgba(240,136,62,.16);color:var(--fg)}
.note{color:var(--muted);font-size:13px}
.legend{display:flex;align-items:center;gap:8px;font-size:12.5px;color:var(--muted);margin:6px 0 10px}
.bar{width:140px;height:10px;border-radius:5px;background:linear-gradient(90deg,rgba(var(--heat),.1),rgba(var(--heat),1))}
footer{margin-top:28px;border-top:1px solid var(--line);padding-top:14px}
footer h2{font-size:15px;margin:10px 0 4px}footer p,footer li{font-size:13px;color:var(--muted);margin:2px 0}
@media(max-width:640px){h1{font-size:20px}.controls{gap:10px}}
</style>
</head>
<body>
<div class="wrap">
<h1>染色体对应关系表</h1>
<p class="sub">基于两物种间共有基因数的染色体对应 · 数据来源：__SOURCE__ · 生成于 __DATE__</p>

<div class="panel controls">
  <label>行物种 <select id="spa"></select></label>
  <button id="swap" title="交换行列所对应的物种">⇄</button>
  <label>列物种 <select id="spb"></select></label>
  <label>弱对应阈值 <input id="min" type="range" min="1" max="300" value="20"> <b id="minv">20</b> 个基因</label>
  <span class="seg" id="views">
    <button data-v="simple" class="on">简表</button><button data-v="list">条形表</button><button data-v="tri">三物种对照</button><button data-v="heat">热图</button><button data-v="table">明细表</button>
  </span>
</div>

<div class="cards" id="cards"></div>
<div id="view"></div>

<footer>
<h2>如何阅读</h2>
<ul>
<li>简表：和 Excel 原表同样的排法。每行是行物种的一条染色体，右侧每格是“列物种染色体=基因数”，从多到少。</li>
<li>条形表：同样的内容，加上占比条形，低于阈值的弱对应以灰色小字列出。</li>
<li>简表和条形表里，鼠标悬停或点击某条列物种染色体，会高亮它在所有行里的出现位置；点击行首的染色体号可以把这一行固定为链接目标（地址栏会带上位置，可复制引用）。</li>
<li>三物种对照：行物种是参照。它的一条染色体同时对应列物种的某条染色体和第三物种的某条染色体时，小表里的数字是这三条染色体共有的基因数。数字接近两边各自的基因数（绿色），说明那两条染色体对应的是同一批基因，彼此同源；数字接近 0（灰色），说明它们各自对应参照染色体的不同部分，彼此不同源。</li>
<li>热图：每个格子里的数字是行物种的这条染色体与列物种的那条染色体共有的基因数；颜色深浅为该数量占行染色体基因总数的比例。</li>
<li>阈值滑块区分主要对应和零星的小计数（通常是噪声或很小的片段）。低于阈值的对应在对应表里以灰色小字显示，在热图和明细表里隐藏。默认 20。</li>
<li>“一对一”：行染色体只对应一条列染色体，且这条列染色体也只对应它。“一对多”：行染色体分散到多条列染色体上（可能是拆分）。“多对一”：行染色体与其他行染色体共同对应同一条列染色体（可能是融合）。这些是根据数量阈值做的描述性分类，不是对重排事件的推断。</li>
<li>对应关系基于基因名，不是严格的直系同源；染色体编号在不同物种间互相独立，数字相同不代表同源。</li>
</ul>
<h2>数据说明</h2>
<p>由 chrom_pairs.py 根据各物种的“基因名 → 染色体”表计算：只使用在两个物种中都只位于一条染色体上的基因名，线粒体和未定位序列不计。物种名和分组来自 species.tsv。</p>
</footer>
</div>
<script>
const D = __DATA__;
const $ = s => document.querySelector(s);
const nk = v => /^\d+$/.test(v) ? [0, +v, ''] : [1, 0, v];
const nat = (a, b) => { const x = nk(a), y = nk(b); return x[0] - y[0] || x[1] - y[1] || (x[2] < y[2] ? -1 : x[2] > y[2] ? 1 : 0); };
const fmt = n => n.toLocaleString('zh-CN');
const S = { a: D.groups[0].species[0], b: D.groups[0].species[1], c: '', min: 20, mode: 'count', view: 'simple', q: '', sk: 'n', sd: -1, chr: '', pin: '', hover: '' };
const IDX = {}; [...D.alphabet].forEach((ch, i) => IDX[ch] = i);

function spName(id) { return D.species[id][0]; }
function groupOf(id) { return D.groups.find(g => g.species.includes(id)); }

// Count genes by the chromosomes they sit on in the given species (2 or 3 of one group); key = labels joined by tab.
function counts(ids) {
  const gs = ids.map(id => D.genes[id]), cs = ids.map(id => D.chr[id]), m = new Map(), L = gs[0].length;
  for (let i = 0; i < L; i++) {
    let k = '', ok = true;
    for (let j = 0; j < gs.length; j++) { const ch = gs[j][i]; if (ch === '!') { ok = false; break; } k += (j ? '\t' : '') + cs[j][IDX[ch]]; }
    if (ok) m.set(k, (m.get(k) || 0) + 1);
  }
  return m;
}
function pairLinks(a, b) { return [...counts([a, b])].map(([k, n]) => { const [x, y] = k.split('\t'); return { a: x, b: y, n }; }); }

function compute() {
  const A = S.a, B = S.b;
  const links = pairLinks(A, B);
  const totA = {}, totB = {}, allA = {};
  links.forEach(l => { totA[l.a] = (totA[l.a] || 0) + l.n; totB[l.b] = (totB[l.b] || 0) + l.n; (allA[l.a] = allA[l.a] || []).push(l); });
  Object.values(allA).forEach(v => v.sort((x, y) => y.n - x.n || nat(x.b, y.b)));
  const strong = links.filter(l => l.n >= S.min);
  const mapA = {}, mapB = {};
  strong.forEach(l => { (mapA[l.a] = mapA[l.a] || []).push(l); (mapB[l.b] = mapB[l.b] || []).push(l); });
  Object.values(mapA).forEach(v => v.sort((x, y) => y.n - x.n));
  Object.values(mapB).forEach(v => v.sort((x, y) => y.n - x.n));
  const chrA = Object.keys(totA).sort(nat), chrB = Object.keys(totB).sort(nat);
  return { A, B, links, strong, totA, totB, allA, mapA, mapB, chrA, chrB };
}

function copyText(text, done) {
  const fallback = () => {
    $('#fb').innerHTML = '<textarea id="fallback" readonly></textarea>';
    const t = $('#fallback'); t.value = text; t.focus(); t.select();
    $('#copied').textContent = '已选中，请按 Ctrl+C 复制';
  };
  try { navigator.clipboard.writeText(text).then(() => { $('#copied').textContent = done; }, fallback); } catch (e) { fallback(); }
}

function saveHash() {
  try { history.replaceState(null, '', '#a=' + S.a + '&b=' + S.b + (S.chr ? '&chr=' + encodeURIComponent(S.chr) : '')); } catch (e) {}
}
function loadHash() {
  let h = {};
  try { h = Object.fromEntries(new URLSearchParams(location.hash.slice(1))); } catch (e) {}
  if (D.species[h.a] && D.species[h.b] && h.a !== h.b && groupOf(h.a) === groupOf(h.b)) { S.a = h.a; S.b = h.b; S.chr = h.chr || ''; }
}

function kind(c, a) {
  const m = c.mapA[a];
  if (!m) return ['无足够对应', 'b0'];
  if (m.length > 1) return ['一对多', 'b2'];
  return c.mapB[m[0].b].length > 1 ? ['多对一', 'b3'] : ['一对一', 'b1'];
}

function renderCards(c) {
  let one = 0, split = 0, merge = 0, none = 0;
  c.chrA.forEach(a => { const k = kind(c, a)[0]; if (k === '一对一') one++; else if (k === '一对多') split++; else if (k === '多对一') merge++; else none++; });
  const cards = [
    [fmt(c.links.reduce((s, l) => s + l.n, 0)), '两物种共有的基因数'],
    [c.chrA.length, `${spName(c.A)}染色体数`], [c.chrB.length, `${spName(c.B)}染色体数`],
    [one, '一对一'], [split, '一对多（分散到多条）'], [merge, '多对一（共用同一条）'], [none, '无足够对应'],
  ];
  $('#cards').innerHTML = cards.map(([n, t]) => `<div class="card"><b>${n}</b><span>${t}</span></div>`).join('');
}

function renderHeat(c) {
  const idxA = Object.fromEntries(c.chrA.map((a, i) => [a, i]));
  const cols = c.chrB.filter(b => c.mapB[b]).sort((x, y) => {
    const dx = idxA[c.mapB[x][0].a], dy = idxA[c.mapB[y][0].a];
    return dx - dy || c.mapB[y][0].n - c.mapB[x][0].n;
  });
  const cell = {};
  c.strong.forEach(l => cell[l.a + '\t' + l.b] = l.n);
  let h = `<div class="legend">颜色 = 占${spName(c.A)}该染色体基因总数的比例 <span class="bar"></span> 0 → 100%　<span>${S.mode === 'count' ? '格内为基因数' : '格内为占比（%）'}</span>
    <button id="mode">${S.mode === 'count' ? '改为显示占比' : '改为显示基因数'}</button></div>`;
  h += '<div class="scroll"><table><thead><tr><th>' + spName(c.A) + ' \\ ' + spName(c.B) + '</th>' + cols.map(b => `<th>${b}</th>`).join('') + '</tr></thead><tbody>';
  c.chrA.forEach(a => {
    h += `<tr><th title="共 ${fmt(c.totA[a])} 个基因">${a} <em class="la">${fmt(c.totA[a])}</em></th>`;
    cols.forEach(b => {
      const n = cell[a + '\t' + b];
      if (!n) { h += '<td class="h"></td>'; return; }
      const f = n / c.totA[a], al = 0.1 + 0.9 * f;
      h += `<td class="h" style="background:rgba(var(--heat),${al.toFixed(2)});color:${f > 0.5 ? '#fff' : 'inherit'}" title="${spName(c.A)} ${a} ↔ ${spName(c.B)} ${b}：${fmt(n)} 个基因，占 ${(100 * f).toFixed(1)}%">${S.mode === 'count' ? fmt(n) : Math.round(100 * f)}</td>`;
    });
    h += '</tr>';
  });
  $('#view').innerHTML = h + '</tbody></table></div>';
  $('#mode').onclick = () => { S.mode = S.mode === 'count' ? 'pct' : 'count'; render(); };
}

function applyHl() {
  const b = S.hover || S.pin;
  document.querySelectorAll('#view .lk').forEach(el => el.classList.toggle('hl', !!b && el.dataset.b === b));
}

function bindLinks(c) {
  const v = $('#view');
  v.onmouseover = e => { const t = e.target.closest('.lk'); const b = t ? t.dataset.b : ''; if (b !== S.hover) { S.hover = b; applyHl(); } };
  v.onmouseleave = () => { S.hover = ''; applyHl(); };
  v.onclick = e => {
    const an = e.target.closest('a.anc');
    if (an) {
      e.preventDefault(); S.chr = S.chr === an.dataset.chr ? '' : an.dataset.chr; saveHash();
      v.querySelectorAll('tr[data-row]').forEach(tr => tr.classList.toggle('target', tr.dataset.row === S.chr));
      return;
    }
    const t = e.target.closest('.lk');
    if (t) { S.pin = S.pin === t.dataset.b ? '' : t.dataset.b; applyHl(); }
  };
  applyHl();
  $('#xl').onclick = () => {
    const rows = c.chrA.map(a => [a].concat(c.allA[a].filter(l => l.n >= S.min).map(l => `${l.b}=${l.n}`)).join('\t'));
    copyText([spName(c.A) + '\t' + spName(c.B)].concat(rows).join('\n'), `已复制 ${rows.length} 行（≥ ${S.min} 个基因的对应）`);
  };
}
function rowHead(c, a) {
  return `<tr class="${S.chr === a ? 'target' : ''}" data-row="${a}"><th><a class="anc" data-chr="${a}" href="#a=${c.A}&b=${c.B}&chr=${encodeURIComponent(a)}" title="把这一行设为链接目标">${a}</a></th>`;
}

// The plain table: same layout as the Excel sheet, one 'chromosome=gene count' cell per correspondence.
function renderSimple(c) {
  const A = spName(c.A), B = spName(c.B);
  const rows = c.chrA.map(a => c.allA[a].filter(l => l.n >= S.min));
  const w = Math.max(1, ...rows.map(r => r.length));
  let h = `<div class="controls" style="margin-bottom:10px"><button id="xl">复制为 Excel 表</button><span id="copied"></span>
    <span class="note">每格为“${B}染色体=基因数”，从多到少 · 只列出不少于 ${S.min} 个基因的对应，把阈值调到 1 可显示全部</span></div><div id="fb"></div>`;
  h += `<div class="scroll wiki simple"><table><thead><tr><th>${A}</th><th colspan="${w}">${B}</th></tr></thead><tbody>`;
  c.chrA.forEach((a, i) => {
    h += rowHead(c, a) + rows[i].map(l => `<td><span class="lk s" data-b="${l.b}" title="${A} ${a} ↔ ${B} ${l.b}：${fmt(l.n)} 个基因">${l.b}=${l.n}</span></td>`).join('') +
      '<td></td>'.repeat(w - rows[i].length) + '</tr>';
  });
  $('#view').innerHTML = h + '</tbody></table></div>';
  bindLinks(c);
}

// Three species: A is the reference; for each chromosome of A, are the B and C chromosomes it corresponds to the same genes?
function renderTri(c) {
  const A = spName(c.A), B = spName(c.B);
  const others = groupOf(c.A).species.filter(id => id !== c.A && id !== c.B);
  if (!others.includes(S.c)) S.c = others[0] || '';
  if (!S.c) { $('#view').innerHTML = '<p class="note">这个分组里没有第三个物种。</p>'; return; }
  const C = spName(S.c), mapC = {};
  pairLinks(c.A, S.c).filter(l => l.n >= S.min).forEach(l => (mapC[l.a] = mapC[l.a] || []).push(l));
  Object.values(mapC).forEach(v => v.sort((x, y) => y.n - x.n));
  const tri = counts([c.A, c.B, S.c]);
  let h = `<div class="controls" style="margin-bottom:10px"><label>第三物种 <select id="spc">${others.map(id => `<option value="${id}">${spName(id)} ${D.species[id][1]}</option>`).join('')}</select></label>
    <span class="note">以${A}为参照：每行是${A}的一条染色体，小表的行是它对应的${B}染色体，列是它对应的${C}染色体，格内是三个物种里同时落在这三条染色体上的基因数</span></div>
    <div class="legend"><span class="sw hom"></span>同一批基因：这段在${B}和${C}之间同源　<span class="sw part"></span>部分重叠　<span class="sw no"></span>几乎没有共同基因：两者对应的是${A}染色体的不同部分</div>`;
  h += `<div class="scroll wiki"><table><thead><tr><th>${A}</th><th style="white-space:normal">${B}（行）× ${C}（列）共同基因数</th></tr></thead><tbody>`;
  c.chrA.forEach(a => {
    const bs = c.mapA[a] || [], cs = mapC[a] || [];
    h += rowHead(c, a) + '<td class="links">';
    if (!bs.length || !cs.length) h += `<span class="note">${!bs.length ? B : C}没有不少于 ${S.min} 个基因的对应</span>`;
    else {
      h += `<table class="mini"><tr><th class="corner">${B} ＼ ${C}</th>` + cs.map(l => `<th>${l.b}<small>${fmt(l.n)}</small></th>`).join('') + '</tr>';
      bs.forEach(lb => {
        h += `<tr><th>${lb.b}<small>${fmt(lb.n)}</small></th>` + cs.map(lc => {
          const t = tri.get(a + '\t' + lb.b + '\t' + lc.b) || 0, f = t / Math.min(lb.n, lc.n);
          const cls = f >= 0.5 ? 'hom' : f >= 0.1 ? 'part' : 'no';
          return `<td class="${cls}" title="${A} ${a} 上有 ${fmt(lb.n)} 个基因在${B} ${lb.b}、${fmt(lc.n)} 个基因在${C} ${lc.b}，其中 ${fmt(t)} 个是同一批基因（占较少一方的 ${(100 * f).toFixed(0)}%）">${t ? fmt(t) : '–'}</td>`;
        }).join('') + '</tr>';
      });
      h += '</table>';
    }
    h += '</td></tr>';
  });
  const v = $('#view');
  v.innerHTML = h + '</tbody></table></div>';
  $('#spc').value = S.c;
  $('#spc').onchange = e => { S.c = e.target.value; render(); };
  v.onclick = e => {
    const an = e.target.closest('a.anc');
    if (!an) return;
    e.preventDefault(); S.chr = S.chr === an.dataset.chr ? '' : an.dataset.chr; saveHash();
    v.querySelectorAll('tr[data-row]').forEach(tr => tr.classList.toggle('target', tr.dataset.row === S.chr));
  };
}

function renderList(c) {
  const A = spName(c.A), B = spName(c.B);
  let h = `<div class="controls" style="margin-bottom:10px"><button id="xl">复制为 Excel 表（染色体=基因数）</button><span id="copied"></span>
    <span class="note">条形 = 占该行染色体基因数的比例 · 灰色小字 = 少于 ${S.min} 个基因的弱对应</span></div><div id="fb"></div>`;
  h += `<div class="scroll wiki"><table><thead><tr><th>${A}</th><th class="n">基因数</th><th>类型</th><th>对应的${B}染色体（基因数从多到少）</th></tr></thead><tbody>`;
  c.chrA.forEach(a => {
    const [k, cls] = kind(c, a);
    const items = c.allA[a].map(l => {
      const f = l.n / c.totA[a], tip = `${A} ${a} ↔ ${B} ${l.b}：${fmt(l.n)} 个基因，占${A} ${a} 的 ${(100 * f).toFixed(1)}%，占${B} ${l.b} 的 ${(100 * l.n / c.totB[l.b]).toFixed(1)}%`;
      return l.n >= S.min
        ? `<span class="lk" data-b="${l.b}" title="${tip}"><b>${l.b}</b><span class="tr"><i style="width:${(100 * f).toFixed(1)}%"></i></span>${fmt(l.n)}<small>${(100 * f).toFixed(0)}%</small></span>`
        : `<span class="lk weak" data-b="${l.b}" title="${tip}">${l.b}=${l.n}</span>`;
    }).join('');
    h += rowHead(c, a) + `<td class="n">${fmt(c.totA[a])}</td><td><span class="badge ${cls}">${k}</span></td><td class="links">${items}</td></tr>`;
  });
  $('#view').innerHTML = h + '</tbody></table></div>';
  bindLinks(c);
}

function tableRows(c) {
  const q = S.q.trim();
  let r = c.strong.map(l => ({ a: l.a, b: l.b, n: l.n, pa: l.n / c.totA[l.a], pb: l.n / c.totB[l.b] }));
  if (q) r = r.filter(x => x.a === q || x.b === q);
  r.sort((x, y) => {
    let v = S.sk === 'a' || S.sk === 'b' ? nat(x[S.sk], y[S.sk]) : x[S.sk] - y[S.sk];
    return (v || x.n - y.n) * S.sd;
  });
  return r;
}

function renderTable(c) {
  const r = tableRows(c);
  const cols = [['a', `${spName(c.A)}染色体`], ['b', `${spName(c.B)}染色体`], ['n', '基因数'], ['pa', `占${spName(c.A)}染色体`], ['pb', `占${spName(c.B)}染色体`]];
  let h = `<div class="controls" style="margin-bottom:10px"><input id="q" type="search" placeholder="按染色体编号筛选，如 3 或 X" value="${S.q}">
    <button id="csv">复制表格（可粘贴到 Excel）</button><span id="copied"></span><span class="note">共 ${r.length} 条对应（≥ ${S.min} 个基因）</span></div><div id="fb"></div>`;
  h += '<div class="scroll"><table class="t"><thead><tr>' + cols.map(([k, t]) => `<th class="${k === S.sk ? 'sorted' : ''}" data-k="${k}" data-arrow="${S.sd < 0 ? '↓' : '↑'}">${t}</th>`).join('') + '</tr></thead><tbody>';
  r.forEach(x => h += `<tr><td>${x.a}</td><td>${x.b}</td><td class="n">${fmt(x.n)}</td><td class="n">${(100 * x.pa).toFixed(1)}%</td><td class="n">${(100 * x.pb).toFixed(1)}%</td></tr>`);
  $('#view').innerHTML = h + '</tbody></table></div>';
  $('#q').oninput = e => { S.q = e.target.value; const p = e.target.selectionStart; renderTable(c); const el = $('#q'); el.focus(); el.setSelectionRange(p, p); };
  document.querySelectorAll('th[data-k]').forEach(th => th.onclick = () => { const k = th.dataset.k; S.sd = S.sk === k ? -S.sd : -1; S.sk = k; renderTable(c); });
  $('#csv').onclick = () => {
    const head = [`${spName(c.A)}染色体`, `${spName(c.B)}染色体`, '基因数', `占${spName(c.A)}染色体比例`, `占${spName(c.B)}染色体比例`];
    const text = [head.join('\t')].concat(r.map(x => [x.a, x.b, x.n, x.pa.toFixed(4), x.pb.toFixed(4)].join('\t'))).join('\n');
    copyText(text, `已复制 ${r.length} 行`);
  };
}

function render() {
  const c = compute();
  const v = $('#view'); v.onmouseover = v.onmouseleave = v.onclick = null;
  renderCards(c);
  ({ simple: renderSimple, tri: renderTri, heat: renderHeat, list: renderList, table: renderTable })[S.view](c);
}

function fillSelects() {
  const opt = id => `<option value="${id}">${spName(id)} ${D.species[id][1]}</option>`;
  $('#spa').innerHTML = D.groups.map(g => `<optgroup label="${g.label}">${g.species.map(opt).join('')}</optgroup>`).join('');
  $('#spa').value = S.a;
  $('#spb').innerHTML = groupOf(S.a).species.filter(id => id !== S.a).map(opt).join('');
  $('#spb').value = S.b;
}

function pick(a, b) {
  const others = groupOf(a).species.filter(id => id !== a);
  S.a = a; S.b = others.includes(b) ? b : others[0];
  S.q = S.chr = S.pin = S.hover = '';
  fillSelects(); saveHash(); render();
}

function init() {
  loadHash();
  fillSelects();
  $('#spa').onchange = e => pick(e.target.value, S.b);
  $('#spb').onchange = e => pick(S.a, e.target.value);
  $('#swap').onclick = () => pick(S.b, S.a);
  $('#min').oninput = e => { S.min = +e.target.value; $('#minv').textContent = S.min; render(); };
  document.querySelectorAll('#views button').forEach(b => b.onclick = () => {
    S.view = b.dataset.v; document.querySelectorAll('#views button').forEach(x => x.classList.toggle('on', x === b)); render();
  });
  render();
  const t = document.querySelector('#view tr.target');
  if (t) t.scrollIntoView({ block: 'center' });
}
init();
</script>
</body>
</html>
'''


def main():
    data = json.load(open(sys.argv[1], encoding='utf-8'))
    pairs = [1] * sum(len(g['species']) * (len(g['species']) - 1) // 2 for g in data['groups'])
    html = (HTML.replace('__DATA__', json.dumps(data, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/'))
            .replace('__SOURCE__', f'{len(data["species"])} 个物种的基因表（chrom_pairs.py 计算）').replace('__DATE__', date.today().isoformat()))
    if len(sys.argv) > 3 and sys.argv[3] == '--artifact':  # the viewer wraps the page in its own skeleton
        html = html.replace('<!doctype html>\n<html lang="zh">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1">\n',
                            '<meta name="claude:viewport-fit-sides" content="cover">\n')
        html = html.replace('</head>\n<body>\n', '').replace('</body>\n</html>\n', '')
    open(sys.argv[2], 'w', encoding='utf-8').write(html)
    print('written', sys.argv[2], len(html), 'bytes;', len(pairs), 'pairs')


if __name__ == '__main__':
    main()
