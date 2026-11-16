#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N-species (N>=3) gene-order comparison inside corresponding chromosome segments. Standard library only.

1. Genes with one unique locus in every species form the shared set (LOC<digits>, MT/Un dropped).
   A species given as a legacy 'SYMBOL<tab>CHR' file (no coordinates) still defines segments but is not
   order-compared.
2. Shared genes are grouped by the tuple of chromosomes they sit on in all species; a tuple with >= --min-genes
   genes is one corresponding segment.
3. In each segment every pair of coordinate species is tested for collinearity (one block after noise removal;
   a fully reversed segment is still collinear). Species are then grouped into mutually collinear clusters.
4. Pattern: all collinear / one species differs / two clusters (split) / >=3 arrangements.
   For two clusters the species tree (--tree, polytomies allowed) says which side can be the derived one:
   a side that is a clade (or a union of children of one node) can carry a single change; the other side,
   which then contains the basal lineages, is taken as the ancestral arrangement.
Usage: python compare_orderN.py --species A=A_positions.tsv --species B=B_positions.tsv --species C=C.txt ...
       --tree "(elephant,(hxx,(meerkat,Equus,Balae)))" --out order5/mammal5
"""
import argparse, csv, html, re
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

from compare_order3 import read_positions, blocks, JUNK_CHR, LOC_RE

DEFAULT_TREE = '(elephant,(hxx,(meerkat,Equus,Balae)))'


def read_legacy(path, keep_loc):
    loc = defaultdict(set)
    with open(path, encoding='utf-8-sig', newline='') as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 2:
                continue
            g, c = row[0].strip().upper(), row[1].strip()
            if not g or not c or c.upper() in JUNK_CHR or (not keep_loc and LOC_RE.match(g)):
                continue
            loc[g].add(c)
    return {g: (next(iter(v)), None, None) for g, v in loc.items() if len(v) == 1}


def read_any(path, keep_loc):
    with open(path, encoding='utf-8-sig') as f:
        head = f.readline()
    if 'Chromosome' in head and 'Start' in head:
        return read_positions(path, keep_loc), True
    return read_legacy(path, keep_loc), False


def parse_tree(s):
    toks = re.findall(r'\(|\)|,|[^(),\s]+', s)
    pos = 0

    def node():
        nonlocal pos
        if toks[pos] == '(':
            pos += 1; ch = [node()]
            while toks[pos] == ',':
                pos += 1; ch.append(node())
            pos += 1
            return tuple(ch)
        pos += 1
        return toks[pos - 1]
    return node()


def leaves(t):
    return {t} if isinstance(t, str) else set().union(*(leaves(x) for x in t))


def prune(t, keep):
    if isinstance(t, str):
        return t if t in keep else None
    ch = [c for c in (prune(x, keep) for x in t) if c is not None]
    return None if not ch else ch[0] if len(ch) == 1 else tuple(ch)


def ok_side(t, side):
    """side is a union of whole children of a single node (clade, or part of a polytomy)."""
    if isinstance(t, str):
        return False
    sets = [leaves(c) for c in t]
    inside = [s for s in sets if s <= side]
    if inside and set().union(*inside) == side:
        return True
    return any(ok_side(c, side) for c in t)


def write(path, head, rows):
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(head); w.writerows(rows)


def fmt(s):
    return '+'.join(sorted(s))


def infer(s, t, tree):
    ok_s, ok_t = ok_side(tree, set(s)), ok_side(tree, set(t))
    if ok_s and not ok_t:
        der, anc = s, t
    elif ok_t and not ok_s:
        der, anc = t, s
    elif ok_s and ok_t:
        return 'undetermined', f'树上单次变化可解释，但无法判断哪侧为祖先（{fmt(s)} | {fmt(t)} 位于根的两侧）'
    else:
        return 'incompatible', '与物种树不相容：至少需要两次独立变化（平行变化，或排列判定有误）'
    if len(der) == 1:
        return 'private', f'祖先排列由{fmt(anc)}保留，{fmt(der)}私有衍生'
    return 'shared', f'{fmt(der)}共有衍生，祖先排列由{fmt(anc)}保留'


def sankoff(tree, state):
    """Minimum number of changes on the (n-ary) tree and the set of root states achieving it. state: species -> label."""
    labels = sorted(set(state.values()))
    inf = 10 ** 6

    def cost(t):
        if isinstance(t, str):
            return {s: 0 if s == state[t] else inf for s in labels}
        cs = [cost(c) for c in t]
        return {s: sum(min(c[s], min(c.values()) + 1) for c in cs) for s in labels}
    root = cost(tree)
    best = min(root.values())
    return best, [s for s in labels if root[s] == best]


def infer_multi(clusters, tree):
    state = {n: i for i, c in enumerate(clusters) for n in c}
    best, roots = sankoff(tree, state)
    anc = ' 或 '.join(fmt(clusters[r]) for r in roots)
    homo = '（含平行变化）' if best > len(clusters) - 1 else ''
    if len(roots) == 1:
        derived = '、'.join(fmt(c) for i, c in enumerate(clusters) if i != roots[0])
        return f'最少{best}次变化{homo}；祖先排列={anc}，其余各自衍生：{derived}'
    return f'最少{best}次变化{homo}；祖先排列候选：{anc}（无法唯一确定）'


def sankoff_bin(tree, state):
    """Binary-character Sankoff; leaves with state None (unknown) cost nothing in either state."""
    inf = 10 ** 6

    def cost(t):
        if isinstance(t, str):
            return {s: 0 if state[t] is None or state[t] == s else inf for s in (False, True)}
        cs = [cost(c) for c in t]
        return {s: sum(min(c[s], min(c.values()) + 1) for c in cs) for s in (False, True)}
    root = cost(tree)
    best = min(root.values())
    return best, [s for s in (False, True) if root[s] == best]


def junction_analysis(ordered, od, flag, tree, window):
    """Per-adjacency (junction) inference. A junction is a pair of genes adjacent in at least one species; it is present in a
    species when their rank distance there is <= window. Genes moved by short runs in a species (flagged in most of its
    pairwise comparisons) are unknown there. Each variable junction is polarised on the tree; carriers of the derived state
    are private (1 species), shared (a clade) or homoplastic."""
    k = len(ordered)
    noise = {s: {g for g, pr in flag.items() if sum(1 for q in pr if s in q) > (k - 1) / 2} for s in ordered}
    order = {s: [g for g in od[s] if g not in noise[s]] for s in ordered}
    rk = {s: {g: i for i, g in enumerate(order[s])} for s in ordered}
    cand = set()
    for s in ordered:
        for a, b in zip(order[s], order[s][1:]):
            cand.add((a, b) if a < b else (b, a))
    gains, derived, shared, homo, amb, amb_root, nvar = Counter(), Counter(), Counter(), 0, 0, 0, 0
    events = Counter()  # carrier set of the derived state -> number of junctions (one event per carrier set per segment)
    kids = [leaves(c) for c in tree] if isinstance(tree, tuple) else []
    for a, b in cand:
        st = {}
        for s in ordered:
            ra, rb = rk[s].get(a), rk[s].get(b)
            st[s] = None if ra is None or rb is None else abs(ra - rb) <= window
        known = [v for v in st.values() if v is not None]
        if all(known) or not any(known):
            continue
        nvar += 1
        _, roots = sankoff_bin(tree, st)
        if len(roots) != 1:
            amb += 1
            for ch in kids:  # the states differ exactly between the two sides of the root
                ins = {st[s] for s in ch if st[s] is not None}
                outs = {st[s] for s in ordered if s not in ch and st[s] is not None}
                if len(ins) == 1 and len(outs) == 1 and ins != outs:
                    amb_root += 1
                    break
            continue
        d = not roots[0]
        carriers = {s for s in ordered if st[s] == d}
        events[frozenset(carriers)] += 1
        for s in carriers:
            derived[s] += 1
            gains[s] += d
        if len(carriers) > 1:
            if ok_side(tree, carriers):
                shared[fmt(carriers)] += d
            else:
                homo += 1
    anc = [s for s in ordered if not derived[s]]
    parts = []
    if nvar == 0:
        parts.append('所有物种的相邻关系一致')
    else:
        parts.append(f'变化的连接 {nvar} 个')
        parts.append('祖先排列保留者=' + (fmt(anc) if anc else '无') + ('' if not amb else '（暂定，不含极性不确定连接涉及的物种差异）'))
        own = [f'{s}{gains[s]}' for s in ordered if gains[s]]
        if own:
            parts.append('衍生断点（新增连接）：' + '、'.join(own))
        sh = [f'{g}:{n}' for g, n in shared.items() if n]
        if sh:
            parts.append('共有衍生断点：' + '、'.join(sh))
        if homo:
            parts.append(f'不相邻物种共有的连接 {homo} 个（平行变化）')
        if amb:
            parts.append(f'极性不确定 {amb} 个（其中 {amb_root} 个是树根两侧的差异，需要更远的外群）')
    status = 'determined' if (nvar > 0 and not amb and anc) else 'root_only' if (amb and amb_root == amb) else 'undetermined'
    if nvar == 0:
        status = 'determined'
    return dict(nvar=nvar, anc=anc, amb=amb, amb_root=amb_root, homo=homo, gains=dict(gains), shared=dict(shared), events=events,
                determined=(status == 'determined'), status=status, text='；'.join(parts))


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--species', action='append', required=True, metavar='NAME=FILE', help='positions tsv (ordered) or legacy txt (chromosome only)')
    p.add_argument('--tree', default=DEFAULT_TREE, help='rooted tree, e.g. "(elephant,(hxx,(meerkat,Equus,Balae)))"; polytomies allowed')
    p.add_argument('--min-genes', type=int, default=30)
    p.add_argument('--window', type=int, default=3)
    p.add_argument('--min-run', type=int, default=5)
    p.add_argument('--keep-loc', action='store_true')
    p.add_argument('--out', default='orderN')
    a = p.parse_args()
    specs = [tuple(s.split('=', 1)) for s in a.species]
    names = [n for n, _ in specs]
    pos, has = {}, {}
    for n, f in specs:
        pos[n], has[n] = read_any(f, a.keep_loc)
    ordered = [n for n in names if has[n]]
    if len(ordered) < 3:
        p.error('need at least three species with coordinates')
    tree = prune(parse_tree(a.tree), set(ordered))
    shared = set.intersection(*(set(pos[n]) for n in names))
    groups = defaultdict(list)
    for g in shared:
        groups[tuple(pos[n][g][0] for n in names)].append(g)
    segs = sorted((s for s, gl in groups.items() if len(gl) >= a.min_genes), key=lambda s: (-len(groups[s]), s))
    pairs = list(combinations(ordered, 2))
    rows, brows, recs = [], [], []
    for k, sig in enumerate(segs, 1):
        sid = f'S{k:02}'
        gl = groups[sig]
        od = {n: sorted(gl, key=lambda g, n=n: (pos[n][g][1], pos[n][g][2])) for n in ordered}
        rk = {n: {g: j for j, g in enumerate(od[n])} for n in ordered}
        col, info = {}, {}
        flag = defaultdict(set)
        for x, y in pairs:
            r = [rk[y][g] for g in od[x]]
            bl, noise, d, runs = blocks(r, a.window, a.min_run)
            kept = {i for b_ in bl for i in b_}
            for i in range(len(r)):
                if i not in kept:
                    flag[od[x][i]].add((x, y))
            col[(x, y)] = len(bl) == 1
            info[(x, y)] = (len(bl), noise)
            for b, run in enumerate(runs, 1):
                genes = [od[x][i] for i in bl[b - 1]]
                dd = [r[i] for i in bl[b - 1]]
                brows.append([sid, f'{x}-{y}', b, len(genes), '+' if len(run) < 2 or d[run[-1]] >= d[run[0]] else '-',
                              x, pos[x][genes[0]][0], pos[x][genes[0]][1], pos[x][genes[-1]][2], genes[0], genes[-1],
                              y, pos[y][od[y][dd[0]]][0], min(pos[y][od[y][v]][1] for v in dd), max(pos[y][od[y][v]][2] for v in dd),
                              od[y][dd[0]], od[y][dd[-1]]])
        # clusters of mutually collinear species
        par = {n: n for n in ordered}

        def find(n):
            while par[n] != n:
                par[n] = par[par[n]]; n = par[n]
            return n
        for (x, y), c in col.items():
            if c:
                par[find(x)] = find(y)
        cl = defaultdict(set)
        for n in ordered:
            cl[find(n)].add(n)
        clusters = sorted(cl.values(), key=lambda s: (-len(s), sorted(s)))
        clique = all(col[tuple(sorted(pr, key=ordered.index))] for c in clusters for pr in combinations(sorted(c, key=ordered.index), 2))
        part = ' | '.join(fmt(c) for c in clusters)
        odd = ''
        if not clique:
            code, text = 'ambiguous', '阈值边缘：同簇内并非两两共线'
        elif len(clusters) == 1:
            code, text = 'all_consistent', '所有物种顺序一致'
        elif len(clusters) == 2:
            small, big = clusters[1], clusters[0]
            if len(small) == 1:
                code = 'one_differs'; odd = next(iter(small))
            else:
                code = 'split'
            _, text = infer(small, big, tree)
        else:
            code = 'all_differ' if len(clusters) == len(ordered) else 'multi'
            text = infer_multi(clusters, tree)
        ja = junction_analysis(ordered, od, flag, tree, a.window)
        rows.append([sid, len(gl)] + list(sig) + [f'{info[pr][0]}({info[pr][1]})' for pr in pairs] +
                    [len(clusters), part, code, odd, text, ja['nvar'], fmt(ja['anc']), ja['amb'], ja['homo'],
                     ja['status'], ja['text']])
        recs.append(dict(sid=sid, genes=len(gl), sig=sig, part=part, code=code, odd=odd, text=text, ja=ja,
                         clusters=clusters, info={f'{x}-{y}': info[(x, y)] for x, y in pairs}))
    head = ['Segment', 'Shared_genes'] + [f'{n}_chr' for n in names] + [f'{x}-{y}_blocks(noise)' for x, y in pairs] + \
           ['Arrangements', 'Partition', 'Pattern', 'Odd_species', 'Inference', 'Junction_variable', 'Junction_ancestral_like',
            'Junction_ambiguous', 'Junction_homoplasy', 'Junction_status', 'Junction_inference']
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    write(f'{out}_segments.tsv', head, rows)
    write(f'{out}_blocks.tsv', ['Segment', 'Pair', 'Block', 'Genes', 'Orient', 'X_species', 'X_chr', 'X_start', 'X_end', 'X_first_gene',
                                'X_last_gene', 'Y_species', 'Y_chr', 'Y_start', 'Y_end', 'Y_first_gene', 'Y_last_gene'], brows)
    erows = []
    for r in recs:
        for car, n in sorted(r['ja']['events'].items(), key=lambda x: (-len(x[0]), sorted(x[0]))):
            kind = 'private' if len(car) == 1 else 'shared_clade' if ok_side(tree, car) else 'homoplasy'
            erows.append([r['sid'], fmt(car), len(car), n, kind])
    write(f'{out}_events.tsv', ['Segment', 'Carriers', 'Size', 'Junctions', 'Kind'], erows)
    cnt = Counter(r['code'] for r in recs)
    two = len(ordered) == 4
    label = {'all_consistent': '全部一致', 'one_differs': '一个不同', 'split': '分成两组(各≥2)',
             'multi': '两个不同（一对共排列，另两个各自独有）' if two else '多种排列（部分物种共排列）',
             'all_differ': '三个及以上不同（每个物种排列互异）' if two else '全部互异', 'ambiguous': '阈值边缘'}
    print(f'species: {", ".join(names)}  coordinate species: {", ".join(ordered)}  shared genes: {len(shared)}  '
          f'segments(>={a.min_genes}): {len(recs)}  genes in segments: {sum(len(groups[s]) for s in segs)}')
    print('pruned tree:', tree)
    for c, l in label.items():
        print(f'  {l}: {cnt[c]}')
    write_html(f'{out}_report.html', names, ordered, tree, recs, cnt, label, a, shared, has)
    print('outputs:', f'{out}_segments.tsv', f'{out}_blocks.tsv', f'{out}_report.html')


CSS = ("body{font-family:-apple-system,'Microsoft YaHei',sans-serif;margin:24px auto;max-width:1150px;padding:0 16px;color:#1f2328;background:#fff}"
       "h1{font-size:22px}h2{font-size:18px;margin-top:32px;border-bottom:1px solid #d0d7de;padding-bottom:6px}"
       "table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}th,td{border:1px solid #d0d7de;padding:4px 8px;text-align:left;vertical-align:top}"
       "th{background:#f6f8fa}.num{text-align:right}.note{color:#57606a;font-size:13px}.odd{background:#fff1e5;font-weight:600}.pair{background:#ddf4ff}"
       "@media(prefers-color-scheme:dark){body{background:#0d1117;color:#e6edf3}th{background:#161b22}th,td,h2{border-color:#30363d}.note{color:#8b949e}.odd{background:#3a2a14}.pair{background:#12294a}}")


def tree_str(t):
    return t if isinstance(t, str) else '(' + ','.join(tree_str(x) for x in t) + ')'


def write_html(path, names, ordered, tree, recs, cnt, label, a, shared, has):
    e = html.escape
    h = [f'<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
         f'<title>{len(names)}物种区段内基因顺序</title><style>{CSS}</style></head><body>',
         f'<h1>{len(names)}个物种对应染色体区段内的基因顺序（{e(", ".join(names))}）</h1>']
    nocoord = [n for n in names if not has[n]]
    h.append(f'<p class="note">共有单拷贝基因 {len(shared)} 个；按各物种所在染色体组合划分，≥{a.min_genes} 个基因的组合为一个对应区段，共 {len(recs)} 个。'
             f'顺序比较的物种：{e(", ".join(ordered))}。'
             + (f'<b>{e(", ".join(nocoord))} 没有坐标文件</b>，只参与划分区段，不参与顺序比较。' if nocoord else '')
             + f' 采用的物种树（已删去无坐标物种）：{e(tree_str(tree))}；括号内逗号分隔的多个物种表示未解析（多分支）。'
             '共线 = 两物种间基因排列只有一块（整段反向仍算共线），断点用窗口 '
             f'{a.window}、小于 {a.min_run} 个基因的零星错位视为噪声。基于基因名对应，非严格直系同源。</p>')
    h.append('<h2>总览</h2><table><tr><th>类型</th><th class="num">区段数</th></tr>')
    for c, l in label.items():
        h.append(f'<tr><td>{l}</td><td class="num">{cnt[c]}</td></tr>')
    h.append(f'<tr><th>合计</th><th class="num">{len(recs)}</th></tr></table>')
    bad = [r for r in recs if r['code'] not in ('all_consistent', 'ambiguous')]
    old_ok = lambda r: '无法' not in r['text']
    h.append('<h2>断点连接级别推断与整段排列推断的对比</h2>'
             f'<p class="note">在 {len(bad)} 个有分歧的区段里，比较两种方法能否确定祖先排列。'
             '断点级：把每对相邻基因当作一个连接，判断它在每个物种里是否保留，并在物种树上分别推断哪个状态是祖先的；'
             '没有任何衍生连接的物种即“保留祖先排列者”。所有变化连接的极性都能确定才算“可确定”；'
             '“仅根处不确定”表示不确定的连接全部是树根两侧的差异（例如只有一个外群与其余物种不同），除此之外的连接都已确定，'
             '这类不确定只有加入更远的外群才能解决。</p>'
             '<table><tr><th></th><th class="num">断点级可确定</th><th class="num">断点级仅根处不确定</th><th class="num">断点级其他不确定</th></tr>')
    for lab, f in (('整段排列可确定', old_ok), ('整段排列无法确定', lambda r: not old_ok(r))):
        c = [sum(f(r) and r['ja']['status'] == s for r in bad) for s in ('determined', 'root_only', 'undetermined')]
        h.append(f'<tr><th>{lab}</th>' + ''.join(f'<td class="num">{x}</td>' for x in c) + '</tr>')
    h.append('</table>')

    def seg_table(rs, cols):
        t = '<table><tr><th>区段</th><th class="num">基因</th>' + ''.join(f'<th>{e(n)}染色体</th>' for n in names) + \
            '<th>排列分组</th><th>整段排列推断</th><th>断点级推断</th></tr>'
        for r in rs:
            single = {n for c in r['clusters'] if len(c) == 1 for n in c} if r['code'] == 'multi' else {r['odd']}
            paired = {n for c in r['clusters'] if len(c) > 1 for n in c} if r['code'] == 'multi' else set()
            cells = ''.join(f'<td class="{"odd" if n in single else "pair" if n in paired else ""}">{e(c)}</td>' for n, c in zip(names, r['sig']))
            t += (f'<tr><td>{r["sid"]}</td><td class="num">{r["genes"]}</td>{cells}<td>{e(r["part"])}</td><td>{e(r["text"])}</td>'
                  f'<td>{e(r["ja"]["text"])}</td></tr>')
        return t + '</table>'
    one = [r for r in recs if r['code'] == 'one_differs']
    h.append(f'<h2>一个不同（{len(one)}）</h2>')
    if one:
        h.append('<table><tr><th>不同者</th><th class="num">区段数</th><th class="num">涉及基因</th><th>推断</th></tr>')
        grp = Counter((r['odd'], r['text']) for r in one)
        for (o, tx), n in sorted(grp.items(), key=lambda x: (-x[1], x[0])):
            h.append(f'<tr><td>{e(o)}</td><td class="num">{n}</td><td class="num">{sum(r["genes"] for r in one if r["odd"] == o and r["text"] == tx)}</td><td>{e(tx)}</td></tr>')
        h.append('</table>')
        h.append(seg_table(one, names))
    for code in ('split', 'multi', 'all_differ', 'ambiguous'):
        rs = [r for r in recs if r['code'] == code]
        if not rs:
            continue
        h.append(f'<h2>{label[code]}（{len(rs)}）</h2>')
        grp = Counter((r['part'], r['text']) for r in rs)
        h.append('<table><tr><th>排列分组</th><th>推断</th><th class="num">区段数</th><th class="num">涉及基因</th></tr>')
        for (pt, tx), n in sorted(grp.items(), key=lambda x: -x[1])[:40]:
            h.append(f'<tr><td>{e(pt)}</td><td>{e(tx)}</td><td class="num">{n}</td>'
                     f'<td class="num">{sum(r["genes"] for r in rs if r["part"] == pt and r["text"] == tx)}</td></tr>')
        h.append('</table>')
        h.append(seg_table(rs, names))
    h.append('<p class="note">表中染色体格子的颜色：橙色加粗 = 排列与其他物种都不同的物种（一个不同的不同者；两个不同里各自独有的两个），蓝色 = 两个不同里共排列的一对；四种排列互异的区段不着色。</p>')
    h.append('<p class="note">完整逐区段结果见 *_segments.tsv，共线块明细见 *_blocks.tsv。</p></body></html>')
    Path(path).write_text('\n'.join(h), encoding='utf-8')


if __name__ == '__main__':
    main()
