#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Four-species chromosome linkage screening using shared gene Symbols.

Input: four merged SYMBOL<TAB>CHROMOSOME files made by parse_ncbi.py.
No sequence, coordinates, or phylogenetic inference is performed.

Outputs: *_modules4.tsv, *_modules3.tsv, *_events.tsv,
         *_best_events.tsv, *_genes4.tsv, *_qc.tsv and .xlsx.
Default excludes Symbols found on >1 chromosome *within a species*.
Pure Python standard library; no third-party packages needed.
"""
import argparse
import csv
import itertools
import re
from collections import defaultdict, Counter
from pathlib import Path


def nkey(v):
    return tuple((0,int(x)) if x.isdigit() else (1,x.lower()) for x in re.split(r'(\d+)',str(v)))


def read_species(path, exclude_trna=False):
    loc=defaultdict(set)
    with open(path,encoding='utf-8-sig',newline='') as f:
        for row in csv.reader(f,delimiter='\t'):
            if len(row)<2: continue
            sym,chrom=row[0].strip().upper(),row[1].strip()
            if not sym or not chrom or (sym in ('GENE_SYMBOL','SYMBOL') and chrom.lower() in ('chromosome','chrom')): continue
            if exclude_trna and (sym.startswith('TRNA') or sym.startswith('MT-TRN')): continue
            loc[sym].add(chrom)
    return dict(loc)


def put_tsv(path,cols,rows):
    with open(path,'w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow(cols);writer.writerows(rows)


def xlsx_write(path,sheets):
    # Borrowed standalone OOXML writer compatible with current compare3 script.
    import zipfile
    from xml.sax.saxutils import escape
    def col(n):
        t=''
        while n:
            n,r=divmod(n-1,26);t=chr(65+r)+t
        return t
    def sheet_data(headers, rows):
        data=[]
        for j,row in enumerate(itertools.chain([headers],rows),1):
            cells=[]
            for k,v in enumerate(row,1):
                style=' s="1"' if j==1 else ''
                pos=f'{col(k)}{j}'
                if isinstance(v,(int,float)) and not isinstance(v,bool):
                    cells.append(f'<c r="{pos}"{style}><v>{v}</v></c>')
                else:
                    cells.append(f'<c r="{pos}" t="inlineStr"{style}><is><t>{escape(str(v))}</t></is></c>')
            data.append(f'<row r="{j}">{"".join(cells)}</row>')
        last=f'{col(len(headers))}{len(rows)+1}'
        widths=''.join(f'<col min="{k}" max="{k}" width="{min(36,max(14,len(h)+3))}" customWidth="1"/>' for k,h in enumerate(headers,1))
        return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                f'<dimension ref="A1:{last}"/><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
                f'<cols>{widths}</cols><sheetData>{"".join(data)}</sheetData><autoFilter ref="A1:{last}"/></worksheet>')
    names=[name for name,_,_ in sheets]
    book=('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'+
          ''.join(f'<sheet name="{escape(s)}" sheetId="{i}" r:id="rId{i}"/>' for i,s in enumerate(names,1))+'</sheets></workbook>')
    rels=('<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'+
          ''.join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,len(sheets)+1))+'</Relationships>')
    types=('<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'+
           ''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,len(sheets)+1))+'</Types>')
    styles=('<?xml version="1.0" encoding="UTF-8"?><styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font></fonts><fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF24506F"/></patternFill></fill></fills><borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/></cellXfs></styleSheet>')
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml',types)
        z.writestr('_rels/.rels','<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        z.writestr('xl/workbook.xml',book);z.writestr('xl/_rels/workbook.xml.rels',rels);z.writestr('xl/styles.xml',styles)
        for i,(_,header,rows) in enumerate(sheets,1):z.writestr(f'xl/worksheets/sheet{i}.xml',sheet_data(header,rows))


def classify(mask,labels):
    together=[labels[i] for i,v in enumerate(mask) if v]
    apart=[labels[i] for i,v in enumerate(mask) if not v]
    if len(together)==1:
        return 'One_species_combined', ','.join(together), ','.join(apart)
    if len(together)==2:
        return 'Two_species_combined', ','.join(together), ','.join(apart)
    if len(together)==3:
        return 'Three_species_combined', ','.join(together), ','.join(apart)
    return 'Not_informative','',''


def analyze(paths,labels,mini=30,min_triple=10,exclude_trna=False,keep_dups=False):
    original=[read_species(p,exclude_trna) for p in paths]
    multiple=[{s for s,loc in d.items() if len(loc)>1} for d in original]
    all_names=set().union(*(set(d) for d in original))
    duplicate_symbols=set().union(*multiple)
    # This prevents a multi-copy gene from appearing in several different modules.
    eligible=all_names if keep_dups else all_names-duplicate_symbols
    chrom_sets=[defaultdict(set) for _ in range(4)]
    for k,genes in enumerate(original):
        for sym in eligible & genes.keys():
            for chrname in genes[sym]:chrom_sets[k][chrname].add(sym)

    # The exact shared-Symbol quartets and triples. Quartet counts avoid products of separate pairwise hit numbers.
    qgenes=defaultdict(set)
    tgenes=defaultdict(set)
    for sym in sorted(eligible):
        present=[i for i in range(4) if sym in original[i]]
        if len(present)<3:continue
        for combo in itertools.combinations(present,3):
            chrom_options=[original[i][sym] for i in combo]
            for chromosomes in itertools.product(*chrom_options):
                tgenes[(combo,chromosomes)].add(sym)
        if len(present)==4:
            for four_chrom in itertools.product(*(original[i][sym] for i in range(4))):
                qgenes[four_chrom].add(sym)

    modules={k:v for k,v in qgenes.items() if len(v)>=mini}
    modkeys=sorted(modules, key=lambda k:(-len(modules[k]),tuple(nkey(x) for x in k)))
    # Make event pairs where at least one species shares a chromosome AND at least one differs.
    # The pair's exact 4-way gene supports are disjoint under default unique-symbol filtering.
    events=[]
    for idx,a in enumerate(modkeys):
        for b in modkeys[idx+1:]:
            mask=tuple(a[i]==b[i] for i in range(4))
            if not 0<sum(mask)<4:continue
            etype,shared,split=classify(mask,labels)
            if not etype.startswith(('One_','Two_','Three_')):continue
            n1,n2=len(modules[a]),len(modules[b])
            event={'type':etype,'shared':shared,'split':split,'pattern':''.join('S' if x else 'D' for x in mask),
                   'a':a,'b':b,'n1':n1,'n2':n2,'total':n1+n2,'balance':round(min(n1,n2)/max(n1,n2),4),
                   'event_id':f'E{len(events)+1:05d}'}
            events.append(event)
    events.sort(key=lambda e:(-min(e['n1'],e['n2']),-e['total'],-e['balance'],e['type'],e['a']))
    for i,e in enumerate(events,1):e['event_id']=f'E{i:05d}'

    # Nonoverlapping summary: one best combination per pair of ancestral modules;
    # no module may be used twice in this *greedy presentation only*, not a phylogenetic optimization.
    best=[];used=set()
    for e in events:
        if e['a'] in used or e['b'] in used:continue
        best.append(e);used.update((e['a'],e['b']))

    quartet_rows=[(*c,len(modules[c]),'|'.join(sorted(modules[c]))) for c in modkeys]
    triple_rows=[]
    for (indices,chromosomes),genes in tgenes.items():
        if len(genes)<min_triple:continue
        key=[str('-')]*4
        for i,c in zip(indices,chromosomes):key[i]=c
        # Quartet subset support: number of exact four-way genes in this triple (may be split across fourth chr).
        fourth=next(j for j in range(4) if j not in indices)
        with_four=sum(len(qv) for qkey,qv in qgenes.items() if all(qkey[i]==c for i,c in zip(indices,chromosomes)))
        triple_rows.append((*key,len(genes),with_four,labels[fourth]))
    triple_rows.sort(key=lambda r:(-r[4],r[:4]))

    gene_rows=[]
    for key in modkeys:
        for sym in sorted(modules[key]):gene_rows.append((*key,sym))

    def event_row(e):
        return (e['event_id'],e['type'],e['shared'],e['split'],e['pattern'],*e['a'],*e['b'],e['n1'],e['n2'],e['total'],e['balance'],'Candidate_only')
    event_rows=[event_row(e) for e in events]
    best_rows=[event_row(e) for e in best]

    qc=[]
    for i,label in enumerate(labels):
        numchrom=len(set().union(*(loc for loc in original[i].values()))) if original[i] else 0
        qc.append((label,str(paths[i]),numchrom,len(original[i]),len(multiple[i]),len(set(original[i])&eligible)))
    qc.append(('ALL','unique-or-multicopy', '',len(all_names),len(duplicate_symbols),len(eligible)))
    return {'modules':quartet_rows,'triples':triple_rows,'events':event_rows,'best':best_rows,'genes':gene_rows,'qc':qc,
            'eligible':len(eligible),'quartet_total':len(qgenes),'quadpass':len(modules),'event_total':len(events),'best_total':len(best)}


def main():
    ap=argparse.ArgumentParser(description='Four-species conserved linkage and candidate structural differences, based on Gene Symbol.')
    for x in ('a','b','c','d'):ap.add_argument(x,help='Merged SYMBOL<TAB>CHROMOSOME file')
    ap.add_argument('--names',nargs=4,metavar=('A','B','C','D'),help='Species display names; defaults to filenames')
    ap.add_argument('--out',default='fish4')
    ap.add_argument('--min-block',type=int,default=30,help='Minimum 4-way symbols per module (default:30)')
    ap.add_argument('--min-triple',type=int,default=10,help='Minimum 3-way symbols for supplementary triple table (default:10)')
    ap.add_argument('--min-pair',type=int,default=None,help='Compatibility alias for --min-triple')
    ap.add_argument('--exclude-trna',action='store_true')
    ap.add_argument('--include-duplicates',action='store_true',help='Allow symbols repeated on multiple chromosomes (can inflate counts)')
    args=ap.parse_args()
    paths=[Path(x) for x in (args.a,args.b,args.c,args.d)]
    for p in paths:
        if not p.is_file():ap.error(f'File not found: {p}')
    labels=args.names or [p.stem for p in paths]
    if len(set(labels))<4:ap.error('Species display names must be different (pass --names)')
    if args.min_block<1 or args.min_triple<1:ap.error('Minimum counts must be positive')
    triple_cutoff=args.min_pair if args.min_pair is not None else args.min_triple
    if triple_cutoff<1:ap.error('--min-pair must be positive')
    res=analyze(paths,labels,args.min_block,triple_cutoff,args.exclude_trna,args.include_duplicates)
    col4=[f'{name}_chr' for name in labels]
    sheets=[
        ('Modules4',col4+['ABC D_shared_count','Gene_symbols'],res['modules']),
        ('Modules3',col4+['Triple_count','With_fourth_count','Missing_species'],res['triples']),
        ('Events',['ID','Pattern_type','Combined_species','Separated_species','S_equal_D_different',
                    *('Module1_'+n for n in col4),*('Module2_'+n for n in col4),
                    'Support1','Support2','Sum_support','Balance','Status'],res['events']),
        ('Best_events',['ID','Pattern_type','Combined_species','Separated_species','S_equal_D_different',
                    *('Module1_'+n for n in col4),*('Module2_'+n for n in col4),
                    'Support1','Support2','Sum_support','Balance','Status'],res['best']),
        ('Genes4',col4+['Gene_symbol'],res['genes']),
        ('QC',['Species','Input','Chromosome_count','Total_symbols','Multichrom_symbols','Symbols_after_filter'],res['qc'])
    ]
    prefix=Path(args.out)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    paths_out={}
    for name,columns,rows in sheets:
        namefile=name.lower()
        p=Path(f'{prefix}_{namefile}.tsv');put_tsv(p,columns,rows);paths_out[name]=p
    book=Path(f'{prefix}.xlsx');xlsx_write(book,sheets)
    print('物种顺序: '+' → '.join(labels))
    print(f'有效Gene Symbol：{res["eligible"]:,}；达到阈值的四物种模块：{res["quadpass"]}（未经筛选共 {res["quartet_total"]} 组）')
    print(f'候选模块组合: {res["event_total"]}；去重展示（贪心、不代表独立事件）: {res["best_total"]}')
    print(f'Excel: {book}（'+', '.join(n for n,_,_ in sheets)+'）')
    print('TSV: '+', '.join(str(p) for p in paths_out.values()))
    print('说明：同名Symbol不是严格直系同源；模块组合不是已确定的融合/裂变，方向需外群与基因位置。')

if __name__=='__main__':main()
