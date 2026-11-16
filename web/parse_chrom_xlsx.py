#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse the pairwise chromosome-correspondence list workbook (cells 'partner_chr=gene_count') into long format.
Usage: python parse_chrom_xlsx.py 20260611矛尾鱼.xlsx out.json
Blocks are located automatically: two adjacent header cells (species A over its chromosome-label column, species B over the
partner cells), optionally followed by a Latin-name header row; data run down until the label column is empty.
"""
import json, re, sys
import openpyxl

DATA = re.compile(r'^(\w+)=(\d+)$')
SPECIES = {  # any header text -> (id, Chinese, Latin)
    '海牛': ('manatee', '海牛', 'Trichechus inunguis'), 'Amazon manatee': ('manatee', '海牛', 'Trichechus inunguis'),
    '儒艮': ('dugong', '儒艮', 'Dugong dugon'), 'Dugong dugon': ('dugong', '儒艮', 'Dugong dugon'),
    '大象': ('elephant', '大象', 'Loxodonta africana'), 'Loxodonta africana': ('elephant', '大象', 'Loxodonta africana'),
    '矛尾鱼': ('latimeria', '矛尾鱼', 'Latimeria chalumnae'), 'Latimeria': ('latimeria', '矛尾鱼', 'Latimeria chalumnae'),
    '棱皮龟': ('dermochelys', '棱皮龟', 'Dermochelys coriacea'), 'Dermochelys': ('dermochelys', '棱皮龟', 'Dermochelys coriacea'),
    '雀鳝': ('lepisosteus', '雀鳝', 'Lepisosteus'), 'Lepisosteus': ('lepisosteus', '雀鳝', 'Lepisosteus'),
}


def is_label(v):
    return v is not None and re.fullmatch(r'\d+|[XYxy]', str(v).strip()) is not None


def is_header(v):
    return isinstance(v, str) and v.strip() in SPECIES


def parse_sheet(ws):
    blocks = []
    for r in range(1, ws.max_row + 1):
        for c in range(1, ws.max_column):
            a, b = ws.cell(r, c).value, ws.cell(r, c + 1).value
            if not (is_header(a) and is_header(b)):
                continue
            nxt = ws.cell(r + 1, c).value
            if is_header(nxt):          # this row holds Chinese names, the next one Latin names: handle at the Latin row
                continue
            first = r + 1
            rows = []
            rr = first
            while ws.cell(rr, c).value is not None and is_label(ws.cell(rr, c).value):
                a_chr = str(ws.cell(rr, c).value).strip()
                cc = c + 1
                while ws.cell(rr, cc).value is not None:
                    m = DATA.match(str(ws.cell(rr, cc).value).strip())
                    if not m:   # next block's chromosome-label column starts here
                        break
                    rows.append((a_chr, m.group(1), int(m.group(2))))
                    cc += 1
                rr += 1
            blocks.append(dict(a=SPECIES[a.strip()][0], b=SPECIES[b.strip()][0], rows=rows, where=f'{ws.title}!R{r}C{c}'))
    return blocks


def main():
    wb = openpyxl.load_workbook(sys.argv[1])
    allb = []
    for ws in wb:
        allb += [dict(b, sheet=ws.title) for b in parse_sheet(ws)]
    for b in allb:
        print(b['sheet'], b['where'], b['a'], b['b'], len(b['rows']), 'links', sum(n for _, _, n in b['rows']), 'genes')
    # Sheet3 should be a thresholded subset of the matching Sheet2 block
    key = lambda b: (b['a'], b['b'])
    full = {key(b): b for b in allb if b['sheet'] == 'Sheet2'}
    for b in allb:
        if b['sheet'] == 'Sheet3' and key(b) in full:
            fs = {(x, y): n for x, y, n in full[key(b)]['rows']}
            miss = [(x, y, n) for x, y, n in b['rows'] if fs.get((x, y)) != n]
            print('Sheet3 vs Sheet2', key(b), 'subset check:', 'OK' if not miss else f'{len(miss)} differ {miss[:3]}')
    json.dump(allb, open(sys.argv[2], 'w', encoding='utf-8'), ensure_ascii=False)


if __name__ == '__main__':
    main()
