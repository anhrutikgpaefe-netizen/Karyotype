#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NCBI gene TSV -> legacy SYMBOL/CHROMOSOME plus optional coordinate table."""
import argparse
import csv
import re
from pathlib import Path
from collections import defaultdict


def natural_key(s):
    return tuple((0,int(p)) if p.isdigit() else (1,p.casefold()) for p in re.split(r'(\d+)',s))


def locate(header, candidates, fallback=None):
    lookup={re.sub(r'[^a-z0-9]','',h.lower()): i for i,h in enumerate(header)}
    for cand in candidates:
        key=re.sub(r'[^a-z0-9]','',cand.lower())
        if key in lookup: return lookup[key]
    return fallback


def main():
    ap=argparse.ArgumentParser(description='Parse NCBI gene TSV; preserve chromosome positions where provided')
    ap.add_argument('input'); ap.add_argument('species',nargs='?'); ap.add_argument('output_dir',nargs='?')
    args=ap.parse_args()
    src=Path(args.input); species=args.species or src.stem; dest=Path(args.output_dir) if args.output_dir else src.parent
    dest.mkdir(parents=True,exist_ok=True)
    with src.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f, delimiter='\t')
        header=next(reader)
        c=locate(header,['Chromosome'],3)
        s=locate(header,['Symbol','Gene Symbol'],6)
        st=locate(header,['Start','Start Position','Start_position','Chromosome Start','Genomic Start','Begin'])
        en=locate(header,['End','End Position','End_position','Chromosome End','Genomic End','Stop'])
        strand=locate(header,['Strand','Orientation'])
        gid=locate(header,['GeneID','Gene ID','NCBI Gene ID'])
        assembly=locate(header,['Assembly','Assembly Accession','Assembly accession'])
        print('检测到列:', ', '.join(f'{name}={header[i] if i is not None and i<len(header) else "无"}' for name,i in [('chrom',c),('symbol',s),('start',st),('end',en),('strand',strand),('gene_id',gid),('assembly',assembly)]))
        positions=[]; legacy=defaultdict(set); skipped=defaultdict(int)
        for row in reader:
            if len(row)<=max(c,s): skipped['invalid_row']+=1; continue
            chrom=row[c].strip(); symbol=row[s].strip().upper()
            if not chrom or not symbol or symbol=='-': skipped['missing_symbol_or_chrom']+=1; continue
            legacy[chrom].add(symbol)
            def val(idx): return row[idx].strip() if idx is not None and idx<len(row) else ''
            x,y=val(st),val(en)
            try:
                a=int(x.replace(',','')); b=int(y.replace(',',''))
                if a>b: a,b=b,a
                if a<0: raise ValueError()
            except ValueError:
                skipped['missing_or_invalid_coordinates']+=1
                continue
            positions.append([chrom,symbol,a,b,val(strand),val(gid),val(assembly)])
    old=dest/(species+'.txt')
    with old.open('w',encoding='utf-8',newline='') as f:
        for chrom in sorted(legacy,key=natural_key):
            for symbol in sorted(legacy[chrom]): f.write(f'{symbol}\t{chrom}\n')
    print(f'兼容旧文件: {old}, {sum(map(len,legacy.values()))} 个染色体-Symbol记录')
    if st is None or en is None:
        print('警告: 原始TSV没有检测到 Start/End 坐标列；不能生成物理共线性文件。请从NCBI获取带基因坐标的数据。')
    elif positions:
        ordered=dest/(species+'_positions.tsv')
        positions.sort(key=lambda r:(natural_key(r[0]),r[2],r[3],r[1]))
        with ordered.open('w',encoding='utf-8-sig',newline='') as f:
            wr=csv.writer(f,delimiter='\t');wr.writerow(['Chromosome','Symbol','Start','End','Strand','GeneID','Assembly','Gene_order'])
            counters=defaultdict(int)
            for row in positions:
                counters[row[0]]+=1; wr.writerow(row+[counters[row[0]]])
        print(f'坐标文件: {ordered}, {len(positions)} 条有坐标记录')
        print('提示: 同名基因可能有多处基因座，坐标文件保留重复记录；compare 会过滤不唯一 Symbol。')
    else:
        print('警告: 没有有效坐标，未生成 positions.tsv')
    print('跳过统计:',dict(skipped))

if __name__=='__main__': main()
