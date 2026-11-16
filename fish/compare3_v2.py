#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Three-species Gene Symbol chromosome comparison (standard library only).

Inputs are parse_ncbi.py merged TSV files: SYMBOL<TAB>CHROMOSOME.
A/B/C chromosome pair counts use set intersections (no arbitrary overwrite).
Triple counts use the *same* gene symbols in all three chromosomes.
"""
import argparse
import csv
import re
import zipfile
from xml.sax.saxutils import escape
from collections import defaultdict
from pathlib import Path


def natural_key(value):
    return tuple((0, int(p)) if p.isdigit() else (1, p.casefold())
                 for p in re.split(r'(\d+)', value))


def read_chromosomes(path, exclude_trna=False):
    result = defaultdict(set)
    with open(path, encoding='utf-8-sig', newline='') as f:
        for lineno, row in enumerate(csv.reader(f, delimiter='\t'), 1):
            if len(row) < 2:
                continue
            symbol, chrom = row[0].strip().upper(), row[1].strip()
            if not symbol or not chrom or (symbol in {'GENE_SYMBOL', 'SYMBOL'} and chrom.lower() in {'chromosome', 'chrom'}):
                continue
            if exclude_trna and (symbol.startswith('TRNA') or symbol.startswith('MT-TRN')):
                continue
            result[chrom].add(symbol)
    return dict(result)


def pair_counts(first, second):
    # Index all occurrences instead of overwriting cross-chromosome duplicates.
    loc = defaultdict(set)
    for chrom, symbols in second.items():
        for gene in symbols:
            loc[gene].add(chrom)
    out = defaultdict(int)
    for c1, symbols in first.items():
        for gene in symbols:
            for c2 in loc.get(gene, ()):
                out[(c1, c2)] += 1
    return out


def write_tsv(path, headers, rows):
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        writer.writerows(rows)



def _col(n):
    value = ''
    while n:
        n, digit = divmod(n - 1, 26)
        value = chr(65 + digit) + value
    return value


def write_xlsx(path, sheets):
    """Write real Excel XLSX with stdlib only. Sheets: [(name, headers, rows)]."""
    def sheet_xml(headers, rows):
        width = [max(12, min(22, len(h) + 3)) for h in headers]
        cols = ''.join(f'<col min="{i}" max="{i}" width="{w}" customWidth="1"/>'
                       for i, w in enumerate(width, 1))
        xml_rows = []
        for ri, row in enumerate([headers, *rows], 1):
            cells = []
            for ci, value in enumerate(row, 1):
                ref = f'{_col(ci)}{ri}'
                style = ' s="1"' if ri == 1 else ''
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    cells.append(f'<c r="{ref}"{style}><v>{value}</v></c>')
                else:
                    cells.append(f'<c r="{ref}" t="inlineStr"{style}><is><t>{escape(str(value))}</t></is></c>')
            xml_rows.append(f'<row r="{ri}">{"".join(cells)}</row>')
        last = f'{_col(len(headers))}{len(rows)+1}'
        return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                f'<dimension ref="A1:{last}"/><sheetViews><sheetView workbookViewId="0">'
                '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
                '</sheetView></sheetViews><sheetFormatPr defaultRowHeight="16"/>'
                f'<cols>{cols}</cols><sheetData>{"".join(xml_rows)}</sheetData>'
                f'<autoFilter ref="A1:{last}"/></worksheet>').encode('utf-8')
    names = [name for name, _, _ in sheets]
    book = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            '<sheets>' + ''.join(f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>'
                               for i, name in enumerate(names, 1)) + '</sheets></workbook>')
    rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + ''.join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1, len(sheets)+1))
            + '</Relationships>')
    content = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
               '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
               '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
               '<Default Extension="xml" ContentType="application/xml"/>'
               '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
               '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
               + ''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,len(sheets)+1))
               + '</Types>')
    styles = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
              '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
              '<font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font></fonts>'
              '<fills count="3"><fill><patternFill patternType="none"/></fill>'
              '<fill><patternFill patternType="gray125"/></fill>'
              '<fill><patternFill patternType="solid"><fgColor rgb="FF234F70"/><bgColor indexed="64"/></patternFill></fill></fills>'
              '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
              '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
              '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
              '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFill="1" applyFont="1"/></cellXfs>'
              '</styleSheet>')
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml', content)
        z.writestr('_rels/.rels', '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        z.writestr('xl/workbook.xml', book)
        z.writestr('xl/_rels/workbook.xml.rels', rels)
        z.writestr('xl/styles.xml', styles)
        for i, (_, headers, rows) in enumerate(sheets, 1):
            z.writestr(f'xl/worksheets/sheet{i}.xml', sheet_xml(headers, rows))



def detect_shared_blocks(A, B, C, min_triple=20, min_fraction=0.15):
    """Find X/Y chromosome groups sharing one anchor chromosome in Z.

    ABC count is the *same gene set* in all three species. All three
    arrangements are tested, regardless of command-line input order.
    Returns rows and small pair-level summaries. They are candidates only.
    """
    import itertools
    species = {'A': A, 'B': B, 'C': C}
    output = []
    for anchor_label in ('A', 'B', 'C'):
        others = [name for name in ('A', 'B', 'C') if name != anchor_label]
        left_label, right_label = others
        X, L, R = species[anchor_label], species[left_label], species[right_label]
        for xchr in sorted(X, key=natural_key):
            triple = []
            for lchr, lg in L.items():
                xy = X[xchr] & lg
                if len(xy) < min_triple:
                    continue
                for rchr, rg in R.items():
                    common = xy & rg
                    count = len(common)
                    if count >= min_triple:
                        triple.append((lchr, rchr, count, len(xy), len(X[xchr] & rg)))
            # Distinct left + distinct right: diagonal supports ancestral blocks;
            # using disjoint members prevents pairing the same chromosome twice.
            for t1, t2 in itertools.combinations(triple, 2):
                l1, r1, n1, ab1, ac1 = t1
                l2, r2, n2, ab2, ac2 = t2
                if l1 == l2 or r1 == r2:
                    continue
                # deterministic one-row-per-diagonal (avoid mirror duplicates)
                if (natural_key(l2), natural_key(r2)) < (natural_key(l1), natural_key(r1)):
                    l1, r1, n1, ab1, ac1, l2, r2, n2, ab2, ac2 = l2, r2, n2, ab2, ac2, l1, r1, n1, ab1, ac1
                cross1 = len(X[xchr] & L[l1] & R[r2])
                cross2 = len(X[xchr] & L[l2] & R[r1])
                purity = (n1+n2)/(n1+n2+cross1+cross2) if n1+n2+cross1+cross2 else 0
                min_strength = min(n1, n2)
                if min_strength / max(n1, n2) < min_fraction:
                    continue
                output.append((anchor_label, xchr, left_label, l1, l2, right_label, r1, r2,
                               n1, n2, cross1, cross2, round(purity, 4),
                               min_strength, n1+n2, 'candidate_only'))
    # Rank within each anchor, prioritizing balanced, high-confidence diagonal blocks
    output.sort(key=lambda z: (-z[13], -z[12], -z[14], z[0], natural_key(z[1])))
    return output



def detect_conserved_pairs(A, B, C, min_triple=20, min_fraction=0.15):
    """Two species retain one chromosome each; third has >=2 strong chromosomes."""
    import itertools
    species = {'A': A, 'B': B, 'C': C}
    out = []
    for split_label in ('A', 'B', 'C'):
        kept = [v for v in ('A', 'B', 'C') if v != split_label]
        k1, k2 = kept
        X, Y, Z = species[k1], species[k2], species[split_label]
        for xchr, xgenes in X.items():
            for ychr, ygenes in Y.items():
                common = xgenes & ygenes
                if len(common) < 2*min_triple:
                    continue
                matches = [(zchr, len(common & zgenes)) for zchr, zgenes in Z.items()]
                matches = [(z, n) for z, n in matches if n >= min_triple]
                if len(matches) < 2:
                    continue
                matches.sort(key=lambda x: (-x[1], natural_key(x[0])))
                top = matches[0][1]
                for (z1, n1), (z2, n2) in itertools.combinations(matches, 2):
                    if min(n1,n2)/max(n1,n2) < min_fraction:
                        continue
                    out.append((k1,xchr,k2,ychr,split_label,z1,z2,n1,n2,
                                n1+n2,len(common),round((n1+n2)/len(common),4),
                                min(n1,n2),'candidate_only'))
    out.sort(key=lambda x: (-x[12],-x[11],-x[9],x[0],natural_key(x[1])))
    return out


PAIR_HEADERS = ['Shared_species_1','Shared_chr_1','Shared_species_2',
                'Shared_chr_2','Split_species','Split_chr_1','Split_chr_2',
                'ABC_count_1','ABC_count_2','Combined_ABC','Shared_AB_total',
                'Fraction_in_two_split_chr','Min_ABC','Interpretation']


BLOCK_HEADERS = ['Anchor_species','Anchor_chr','Other1_species','Other1_chr_1',
                 'Other1_chr_2','Other2_species','Other2_chr_1','Other2_chr_2',
                 'ABC_diagonal_1','ABC_diagonal_2','ABC_cross_1','ABC_cross_2',
                 'Diagonal_fraction','Min_diagonal','Total_diagonal','Interpretation']


def load_positions(path, exclude_trna=False):
    """Read parse_ncbi_v2 _positions.tsv, with one locus per symbol per species.
    Duplicate symbols are conservatively excluded across all chromosomes.
    """
    import csv
    from collections import Counter
    rows=[]
    with open(path,encoding='utf-8-sig',newline='') as f:
        for line in csv.DictReader(f,delimiter='\t'):
            try:
                ch=line['Chromosome'].strip(); symbol=line['Symbol'].strip().upper()
                start=int(line['Start']); end=int(line['End'])
                if symbol and ch and (not exclude_trna or not (symbol.startswith('TRNA') or symbol.startswith('MT-TRN'))):
                    rows.append((ch,symbol,start,end))
            except (ValueError,KeyError,TypeError): continue
    counts=Counter(r[1] for r in rows)
    filtered=[r for r in rows if counts[r[1]]==1]
    result=defaultdict(list)
    for ch,symbol,start,end in filtered: result[ch].append((start,end,symbol))
    for ch in result: result[ch].sort()
    return dict(result),len(rows)-len(filtered)


def position_blocks(first, second, first_label, second_label, window=30, min_anchors=12, min_rho=0.65):
    """Sliding local gene-order correlation; conservative exploratory collinearity candidates.
    Windows require matched unique symbols, a consistent target chromosome and
    predominantly increasing or decreasing target order. Neighbour windows merge.
    """
    from collections import defaultdict
    from math import sqrt
    target={g:(ch,i,start,end) for ch,seq in second.items() for i,(start,end,g) in enumerate(seq)}
    results=[]
    for chr1,seq in first.items():
        anchors=defaultdict(list)
        for i,(start,end,g) in enumerate(seq):
            if g in target:
                chr2,j,tstart,tend=target[g]
                anchors[chr2].append((i,j,start,end,tstart,tend,g))
        for chr2,match in anchors.items():
            if len(match)<min_anchors: continue
            # Each window contains up to 'window' shared unique anchors.
            accepted=[]
            step=max(1,window//3)
            for k in range(0,len(match),step):
                sub=match[k:k+window]
                if len(sub)<min_anchors: continue
                vals=[r[1] for r in sub]; n=len(vals)
                # Spearman's rank: first coordinate is already sorted by anchor order.
                ranks={value:i for i,value in enumerate(sorted(vals))}
                diff=sum((i-ranks[v])**2 for i,v in enumerate(vals))
                rho=1-6*diff/(n*(n*n-1))
                if abs(rho)>=min_rho:
                    accepted.append((k,min(len(match),k+window),rho))
            # Merge overlapping, same-direction candidate windows.
            merged=[]
            for lo,hi,rho in accepted:
                direction='+' if rho>=0 else '-'
                if merged and direction==merged[-1][2] and lo<=merged[-1][1]:
                    merged[-1]=(merged[-1][0],max(merged[-1][1],hi),direction)
                else: merged.append((lo,hi,direction))
            for lo,hi,direction in merged:
                segment=match[lo:hi]
                if len(segment)<min_anchors:continue
                results.append((first_label,chr1,second_label,chr2,direction,len(segment),
                     min(x[2] for x in segment),max(x[3] for x in segment),
                     min(x[4] for x in segment),max(x[5] for x in segment),
                     min(x[0] for x in segment)+1,max(x[0] for x in segment)+1,
                     'local_order_candidate'))
    return sorted(results,key=lambda r:(-r[5],r[0],natural_key(r[1]),r[2],natural_key(r[3])))

POSITION_HEADERS=['Species_1','Chr_1','Species_2','Chr_2','Order_direction','Anchor_count',
 'Start_1','End_1','Start_2','End_2','Gene_index_start_1','Gene_index_end_1','Interpretation']

def main():
    parser = argparse.ArgumentParser(description='Three-species chromosome correspondence and exact triple gene intersections.')
    parser.add_argument('a', help='species A merged SYMBOL\\tCHROMOSOME file')
    parser.add_argument('b', help='species B merged SYMBOL\\tCHROMOSOME file')
    parser.add_argument('c', help='species C merged SYMBOL\\tCHROMOSOME file')
    parser.add_argument('--out', default='compare3', help='output path prefix (default: compare3)')
    parser.add_argument('--min-pair', type=int, default=10, help='minimum A-B and A-C counts to inspect (default 10)')
    parser.add_argument('--min-triple', type=int, default=0, help='minimum ABC count to output a row (default 0)')
    parser.add_argument('--exclude-trna', action='store_true', help='exclude genes beginning with TRNA or MT-TRN')
    parser.add_argument('--min-block', type=int, default=20, help='minimum triple-shared genes per block (default 20)')
    parser.add_argument('--all', action='store_true', help='include all A-B/A-C pairs with count >=1 regardless of --min-pair')
    parser.add_argument('--positions-a', help='A coordinate TSV from parse_ncbi_v2')
    parser.add_argument('--positions-b', help='B coordinate TSV')
    parser.add_argument('--positions-c', help='C coordinate TSV')
    parser.add_argument('--order-window', type=int, default=30, help='matching anchors per local window')
    parser.add_argument('--min-order-anchors', type=int, default=12)
    parser.add_argument('--min-order-rho', type=float, default=0.65)
    args = parser.parse_args()
    if args.min_pair < 1 or args.min_triple < 0 or args.min_block < 1:
        parser.error('thresholds must be nonnegative; --min-pair must be >= 1')
    for path in (args.a, args.b, args.c):
        if not Path(path).is_file():
            parser.error(f'file not found: {path}')
    A, B, C = (read_chromosomes(p, args.exclude_trna) for p in (args.a, args.b, args.c))
    ab, ac, bc = pair_counts(A, B), pair_counts(A, C), pair_counts(B, C)
    threshold = 1 if args.all else args.min_pair
    b_by_a, c_by_a = defaultdict(list), defaultdict(list)
    for (a, b), n in ab.items():
        if n >= threshold:
            b_by_a[a].append(b)
    for (a, c), n in ac.items():
        if n >= threshold:
            c_by_a[a].append(c)

    summary, details = [], []
    for a in sorted(A, key=natural_key):
        for b in sorted(b_by_a.get(a, []), key=natural_key):
            ab_genes = A[a] & B[b]
            for c in sorted(c_by_a.get(a, []), key=natural_key):
                triple = sorted(ab_genes & C[c])
                abc = len(triple)
                if abc < args.min_triple:
                    continue
                n_ab, n_ac, n_bc = ab[(a, b)], ac[(a, c)], bc.get((b, c), 0)
                summary.append((a, b, c, n_ab, n_ac, n_bc, abc,
                                f'{abc/n_ab:.3f}', f'{abc/n_ac:.3f}',
                                f'{abc/n_bc:.3f}' if n_bc else '0.000'))
                details.extend((a, b, c, gene) for gene in triple)

    prefix = Path(args.out)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    summary_path = str(prefix) + '_summary.tsv'
    genes_path = str(prefix) + '_genes.tsv'
    sum_headers = ['A_chr','B_chr','C_chr','AB_count','AC_count','BC_count',
                   'ABC_count','ABC_over_AB','ABC_over_AC','ABC_over_BC']
    gene_headers = ['A_chr','B_chr','C_chr','Gene_symbol']
    write_tsv(summary_path, sum_headers, summary)
    write_tsv(genes_path, gene_headers, details)
    blocks = detect_shared_blocks(A, B, C, min_triple=args.min_block)
    block_path = str(prefix) + '_shared_blocks.tsv'
    write_tsv(block_path, BLOCK_HEADERS, blocks)
    conserved = detect_conserved_pairs(A, B, C, min_triple=args.min_block)
    conserved_path = str(prefix) + '_conserved_pairs.tsv'
    write_tsv(conserved_path, PAIR_HEADERS, conserved)
    excel_summary = [tuple([*row[:7], *[float(x) for x in row[7:]]]) for row in summary]
    excel_path = str(prefix) + '.xlsx'
    write_xlsx(excel_path, [('Summary', sum_headers, excel_summary), ('Genes', gene_headers, details), ('Shared_blocks', BLOCK_HEADERS, blocks), ('Conserved_pairs', PAIR_HEADERS, conserved)])
    positional_paths=[args.positions_a,args.positions_b,args.positions_c]
    if any(positional_paths):
        if not all(positional_paths): parser.error('--positions-a/b/c must all be supplied')
        if args.order_window<3 or args.min_order_anchors<3 or not 0<=args.min_order_rho<=1:
            parser.error('invalid gene-order window/rho thresholds')
        for p in positional_paths:
            if not Path(p).is_file(): parser.error(f'position file not found: {p}')
        positional=[]
        for label,p in zip('ABC',positional_paths):
            data, removed=load_positions(p,args.exclude_trna)
            print(f'{label}坐标记录: {sum(map(len,data.values()))}, 过滤重复Symbol基因座: {removed}')
            positional.append(data)
        orders=[]
        for i,j in [(0,1),(0,2),(1,2)]:
            orders.extend(position_blocks(positional[i],positional[j],'ABC'[i],'ABC'[j],
                              args.order_window,args.min_order_anchors,args.min_order_rho))
        position_path=str(prefix)+'_ordered_blocks.tsv'
        write_tsv(position_path,POSITION_HEADERS,orders)
        print('基因顺序共线性候选:',len(orders),'文件:',position_path)
        print('注意: ordered_blocks是局部基因顺序相关候选，尚不能自动断定祖先状态或结构变异。')
    print('染色体数:', f'A={len(A)}', f'B={len(B)}', f'C={len(C)}')
    print('阈值: AB和AC至少', threshold, '个共享Symbol; ABC至少', args.min_triple)
    print('候选三元组:', len(summary), '三元组中基因记录:', len(details))
    print('Excel文件:', excel_path, '(Summary、Genes、Shared_blocks、Conserved_pairs四个工作表)')
    print('兼容TSV:', summary_path, '和', genes_path)
    print('共有连锁单元候选:', len(blocks), '文件:', block_path)
    print('双物种共有结构候选:', len(conserved), '文件:', conserved_path)
    print('输入顺序不会遗漏模式：程序轮流以A、B、C为锚点。')
    print('提示: 三元交集是相同Symbol的交集，不等同于严格直系同源或物理共线性。')


if __name__ == '__main__':
    main()
