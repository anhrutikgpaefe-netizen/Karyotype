#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Multi-species chromosome linkage and coordinate-aware anchor-window analysis.

Standard library only. Inputs are legacy parse_ncbi_v2 SYMBOL<tab>CHROMOSOME
and optional *_positions.tsv. Never asserts evolutionary direction or strict orthology.
"""
import argparse, csv, json, re, sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path


def nk(v):
    return tuple((0,int(x)) if x.isdigit() else (1,x.lower()) for x in re.split(r'(\d+)',str(v)))


LOC_RE=re.compile(r'^LOC\d+$')
JUNK_CHR={'MT','MITO','MITOCHONDRION','UN','UNPLACED','UNLOCALIZED'}


def skip_symbol(g,exclude_trna,drop_loc):
    # LOC<digits> are NCBI placeholder names, species-specific and not comparable across species.
    if drop_loc and LOC_RE.match(g):return True
    return exclude_trna and (g.startswith('TRNA') or g.startswith('MT-TRN'))


def read_legacy(path, exclude_trna=False, drop_loc=True):
    loc=defaultdict(set)
    with open(path,encoding='utf-8-sig',newline='') as f:
        for row in csv.reader(f,delimiter='\t'):
            if len(row)<2: continue
            g,c=row[0].strip().upper(),row[1].strip()
            if not g or not c or (g in ('SYMBOL','GENE_SYMBOL') and c.lower() in ('chromosome','chrom')): continue
            if skip_symbol(g,exclude_trna,drop_loc):continue
            loc[g].add(c)
    return {g:next(iter(v)) for g,v in loc.items() if len(v)==1},len(loc),sum(len(v)>1 for v in loc.values())


def read_positions(path,exclude_trna=False,drop_loc=True):
    loci=defaultdict(list)
    with open(path,encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        required={'Chromosome','Symbol','Start','End'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: missing columns {required-set(reader.fieldnames or [])}')
        for r in reader:
            g=r['Symbol'].strip().upper(); c=r['Chromosome'].strip()
            if not g or not c or skip_symbol(g,exclude_trna,drop_loc):continue
            try: start=int(r['Start']); end=int(r['End'])
            except (ValueError,TypeError):continue
            loci[g].append((c,start,end))
    # Do not arbitrarily pick among duplicated loci even if on one chromosome.
    return {g:v[0] for g,v in loci.items() if len(v)==1},len(loci),sum(len(v)>1 for v in loci.values())


def write(path,headers,rows):
    with open(path,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(headers);w.writerows(rows)


def read_args():
    p=argparse.ArgumentParser(description='Multi-species (recommended 8) chromosome linkage / ordered module explorer')
    p.add_argument('files',nargs='*',help='eight legacy files in order (alternative to --species)')
    p.add_argument('--species',action='append',default=[],metavar='NAME=FILE',help='repeat 8 times; also supports 3+ species')
    p.add_argument('--positions',action='append',default=[],metavar='NAME=FILE',help='optional position input for each species')
    p.add_argument('--group',action='append',default=[],metavar='GROUP=NAME,NAME',help='optional descriptive clade grouping (not automatic phylogeny)')
    p.add_argument('--min-block',type=int,default=30,help='minimum supported anchor genes per module')
    p.add_argument('--min-species',type=int,default=4,help='minimum species with >= min-support genes (including anchor)')
    p.add_argument('--min-support',type=int,default=8,help='minimum genes for mapping a module to another species chromosome')
    p.add_argument('--window',type=int,default=150,help='anchor coordinate window size in genes')
    p.add_argument('--purity',type=float,default=.65,help='minimum chromosome assignment dominance')
    p.add_argument('--min-genes-present',type=int,default=3,help='required independent species presence for anchor marker')
    p.add_argument('--max-modules',type=int,default=4000)
    p.add_argument('--out',default='mammal8')
    p.add_argument('--exclude-trna',action='store_true')
    p.add_argument('--keep-loc',action='store_true',help='keep LOC<digits> placeholder symbols (dropped by default)')
    p.add_argument('--min-chr-genes',type=int,default=20,help='drop chromosomes/scaffolds with fewer usable genes (default 20); MT/Un are always dropped')
    p.add_argument('--merge-overlap',type=float,default=.80,help='modules from different anchors sharing >= this fraction of genes (vs the smaller one) with no conflicting chromosome are merged as duplicates (default 0.80)')
    a=p.parse_args()
    if a.species and a.files:p.error('use positional files OR --species, not both')
    if a.species:
        specs=[]
        for item in a.species:
            if '=' not in item:p.error('species needs NAME=FILE')
            name,path=item.split('=',1);specs.append((name.strip(),path.strip()))
    else:specs=[(Path(f).stem,f) for f in a.files]
    if len(specs)<3:p.error('at least 3 species required')
    if len(set(x for x,_ in specs))!=len(specs):p.error('species names must be unique')
    if a.min_species>len(specs):p.error('min-species exceeds number of species')
    if a.window<10 or a.min_block<2 or a.min_support<1 or not 0<a.purity<=1:p.error('invalid thresholds')
    if not 0<a.merge_overlap<=1 or a.min_chr_genes<0:p.error('invalid --merge-overlap / --min-chr-genes')
    pos={}
    for item in a.positions:
        if '=' not in item:p.error('--positions needs NAME=FILE')
        n,f=item.split('=',1);pos[n.strip()]=f.strip()
    for n,f in specs:
        if not Path(f).is_file():p.error('missing '+f)
    for n,f in pos.items():
        if n not in dict(specs):p.error('position species not in inputs: '+n)
        if not Path(f).is_file():p.error('missing '+f)
    return a,specs,pos


def main():
    args,specs,pospaths=read_args()
    names=[n for n,_ in specs];loc={}; coords={}; qc=[]
    for name,path in specs:
        drop_loc=not args.keep_loc
        locus,total,dups=read_legacy(path,args.exclude_trna,drop_loc)
        positions={}; ptot=pdup=0
        if name in pospaths:
            positions,ptot,pdup=read_positions(pospaths[name],args.exclude_trna,drop_loc)
            # Coordinate-unique subset only; chromosome assignment from coordinate input.
            locus={g:c for g,(c,s,e) in positions.items()}
        nchr_raw=len(set(locus.values()))
        size=Counter(locus.values())
        bad={c for c,n in size.items() if n<args.min_chr_genes or c.upper() in JUNK_CHR}
        locus={g:c for g,c in locus.items() if c not in bad}
        positions={g:v for g,v in positions.items() if v[0] not in bad}
        loc[name]=locus;coords[name]=positions
        qc.append([name,path,len(locus),total,dups,len(positions),ptot,pdup,len(set(locus.values())), 'coordinate' if positions else 'chromosome_only',nchr_raw,len(bad)])
    presence=Counter(g for n in names for g in loc[n])
    modules=[]
    for anchor in names:
        by_chr=defaultdict(list)
        for g,c in loc[anchor].items():
            if presence[g]>=args.min_genes_present:by_chr[c].append(g)
        for chrom in sorted(by_chr,key=nk):
            genes=by_chr[chrom]
            if coords[anchor]:genes.sort(key=lambda g:(coords[anchor][g][1],g))
            else:genes.sort()
            # In chromosome-only mode this is one unsplit chromosome, NOT alphabetical windows.
            segments=[genes[i:i+args.window] for i in range(0,len(genes),args.window)] if coords[anchor] else [genes]
            window_rows=[]
            for segment in segments:
                signature=[]
                for other in names:
                    if other==anchor:signature.append(chrom);continue
                    ct=Counter(loc[other][g] for g in segment if g in loc[other])
                    if not ct:signature.append('');continue
                    best,count=sorted(ct.items(),key=lambda x:(-x[1],nk(x[0])))[0]
                    signature.append(best if count>=args.min_support and count/sum(ct.values())>=args.purity else '')
                window_rows.append((tuple(signature),segment))
            # merge adjacent physical windows if their dominant chromosome assignments agree
            groups=[]
            for sig,geneseg in window_rows:
                if groups and groups[-1][0]==sig:groups[-1][1].extend(geneseg)
                else:groups.append((sig,list(geneseg)))
            for sig,gset in groups:
                if len(gset)<args.min_block:continue
                assignment=[];counts=[];purities=[]
                for other in names:
                    if other==anchor:
                        assignment.append(chrom);counts.append(len(gset));purities.append(1.0);continue
                    ct=Counter(loc[other][g] for g in gset if g in loc[other])
                    if ct:
                        best,n=sorted(ct.items(),key=lambda x:(-x[1],nk(x[0])))[0]
                        purity=n/sum(ct.values())
                        if n>=args.min_support and purity>=args.purity:
                            assignment.append(best);counts.append(n);purities.append(round(purity,4));continue
                    assignment.append('');counts.append(0);purities.append(0)
                supported=sum(c>=args.min_support for c in counts)
                if supported<args.min_species:continue
                start=min(coords[anchor][g][1] for g in gset) if coords[anchor] else ''
                end=max(coords[anchor][g][2] for g in gset) if coords[anchor] else ''
                modules.append({'anchor':anchor,'anchor_chr':chrom,'start':start,'end':end,'genes':set(gset),'map':assignment,'counts':counts,'purity':purities,'support':supported,'type':'ordered_window' if coords[anchor] else 'chromosome_level'})
    # Dedup across anchor species when overlapping genes and mapping signature agree.
    modules.sort(key=lambda m:(-m['support'],-min(x for x in m['counts'] if x),-len(m['genes'])))
    selected=[];index=defaultdict(list)
    for m in modules:
        # avoid quadratic comparisons using an informative anchor symbol
        if not m['genes']:continue
        witnesses=sorted(m['genes'])[:min(8,len(m['genes']))]
        candidates=set(i for g in witnesses for i in index[g])
        duplicate=False
        for i in candidates:
            old=selected[i]
            shared=len(m['genes'] & old['genes']);overlap=shared/min(len(m['genes']),len(old['genes']))
            conflict=any(a and b and a!=b for a,b in zip(m['map'],old['map']))
            if overlap>=args.merge_overlap and not conflict:duplicate=True;break
        if duplicate:continue
        if len(selected)>=args.max_modules:break
        i=len(selected);selected.append(m)
        for g in m['genes']:index[g].append(i)
    # Pair of modules in at least 2 shared species; chromosomes together/apart.
    events=[]
    for i,m in enumerate(selected):
        if i>=1200:break # bounded candidate exploration
        for j in range(i+1,min(len(selected),1200)):
            q=selected[j]
            if len(m['genes']&q['genes'])>min(len(m['genes']),len(q['genes']))*.20:continue
            together=[names[k] for k in range(len(names)) if m['map'][k] and m['map'][k]==q['map'][k]]
            apart=[names[k] for k in range(len(names)) if m['map'][k] and q['map'][k] and m['map'][k]!=q['map'][k]]
            if not together or not apart or len(together)+len(apart)<args.min_species:continue
            score=min(sum(m['counts']),sum(q['counts'])) * (len(together)+len(apart))/len(names)
            events.append((round(score,2),i,j,','.join(together),','.join(apart),len(together),len(apart)))
    events.sort(reverse=True)
    prefix=Path(args.out);prefix.parent.mkdir(parents=True,exist_ok=True)
    head=['Module_ID','Anchor_species','Anchor_chr','Anchor_start','Anchor_end','Anchor_gene_count','Supported_species','Evidence_level']
    for n in names:head.extend([n+'_chr',n+'_anchor_shared',n+'_dominance'])
    head+=['Anchor_gene_symbols']
    mr=[]
    for i,m in enumerate(selected,1):
        row=[f'M{i:05}',m['anchor'],m['anchor_chr'],m['start'],m['end'],len(m['genes']),m['support'],m['type']]
        for c,num,p in zip(m['map'],m['counts'],m['purity']):row.extend((c,num,p))
        row+=['|'.join(sorted(m['genes']))];mr.append(row)
    er=[]
    for k,(score,i,j,to,apart,tc,ac) in enumerate(events,1):
        er.append([f'E{k:05}',f'M{i+1:05}',f'M{j+1:05}',to,apart,tc,ac,score,'candidate_only'])
    write(str(prefix)+'_modules.tsv',head,mr)
    write(str(prefix)+'_events.tsv',['Event_ID','Module_1','Module_2','Together_species','Apart_species','Together_n','Apart_n','Score','Interpretation'],er)
    write(str(prefix)+'_qc.tsv',['Species','Input','Usable_unique_symbols','Raw_symbols','Duplicate_symbols','Unique_position_loci','Raw_position_symbols','Duplicate_position_symbols','Chromosomes','Analysis_mode','Chromosomes_before_filter','Chromosomes_dropped'],qc)
    groups=[]
    for item in args.group:
        if '=' not in item:raise ValueError('--group expects GROUP=NAME,NAME')
        label,s=item.split('=',1);members=s.split(',')
        if any(m not in names for m in members):raise ValueError('unknown group member '+item)
        groups.append((label,members))
    grows=[]
    for i,m in enumerate(selected,1):
        for label,members in groups:
            matched=[m['map'][names.index(n)] for n in members]
            agrees=bool(all(matched) and all(x and x==matched[0] for x in matched))
            grows.append([f'M{i:05}',label,','.join(members),int(agrees),matched[0] if agrees else '',len([x for x in matched if x])])
    if groups:write(str(prefix)+'_groups.tsv',['Module_ID','Group','Members','Same_chromosome_label_within_group','Label_if_same','Represented_species'],grows)
    cfg={'species':names,'files':dict(specs),'positions':pospaths,'thresholds':{'min_block':args.min_block,'min_species':args.min_species,'min_support':args.min_support,'window':args.window,'purity':args.purity,'merge_overlap':args.merge_overlap,'min_chr_genes':args.min_chr_genes,'drop_loc_symbols':not args.keep_loc},'module_count':len(selected),'event_count':len(events),'warnings':['Comparisons based on Gene Symbols, not validated orthology.','Modules are heuristic; cannot prove physical collinearity or event direction.','Without positions, chromosome-level mapping only.','Across species chromosome labels are independent; group labels are not homologous identifiers.','Events are overlapping candidate pairs, not independent historical rearrangements.']}
    with open(str(prefix)+'_manifest.json','w',encoding='utf-8') as f:json.dump(cfg,f,ensure_ascii=False,indent=2)
    print('物种:', ' → '.join(names))
    print('基因位置输入:',len(pospaths),'/',len(names))
    print('候选模块:',len(selected),'联合/分离组合:',len(events))
    print('结果:',str(prefix)+'_modules.tsv,',str(prefix)+'_events.tsv,',str(prefix)+'_qc.tsv')
    print('提示: chromosome_only 模式不会用字母排序冒充基因物理顺序。结果都是候选，不等于独立演化事件。')

if __name__=='__main__':main()
