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


def read_legacy(path, exclude_trna=False):
    loc=defaultdict(set)
    with open(path,encoding='utf-8-sig',newline='') as f:
        for row in csv.reader(f,delimiter='\t'):
            if len(row)<2: continue
            g,c=row[0].strip().upper(),row[1].strip()
            if not g or not c or (g in ('SYMBOL','GENE_SYMBOL') and c.lower() in ('chromosome','chrom')): continue
            if exclude_trna and (g.startswith('TRNA') or g.startswith('MT-TRN')):continue
            loc[g].add(c)
    return {g:next(iter(v)) for g,v in loc.items() if len(v)==1},len(loc),sum(len(v)>1 for v in loc.values())


def read_positions(path,exclude_trna=False):
    loci=defaultdict(list)
    with open(path,encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        required={'Chromosome','Symbol','Start','End'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: missing columns {required-set(reader.fieldnames or [])}')
        for r in reader:
            g=r['Symbol'].strip().upper(); c=r['Chromosome'].strip()
            if not g or not c or (exclude_trna and (g.startswith('TRNA') or g.startswith('MT-TRN'))):continue
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
    p.add_argument('--merge-overlap',type=float,default=.60,help='minimum smaller-set overlap for same chromosome signature merging')
    p.add_argument('--out',default='mammal8')
    p.add_argument('--exclude-trna',action='store_true')
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
    if a.window<10 or a.min_block<2 or a.min_support<1 or not 0<a.purity<=1 or not 0<a.merge_overlap<=1:p.error('invalid thresholds')
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


def module_signature(m):
    return tuple(m['map'])


def recalc_module(m, names, loc, args):
    """Recompute support after union, never reuse original per-species counts."""
    gset=m['genes']
    assignment=[];counts=[];purities=[]
    for name in names:
        ct=Counter(loc[name][g] for g in gset if g in loc[name])
        if not ct:
            assignment.append('');counts.append(0);purities.append(0.0);continue
        best,n=sorted(ct.items(),key=lambda x:(-x[1],nk(x[0])))[0]
        purity=n/sum(ct.values())
        assignment.append(best if n>=args.min_support and purity>=args.purity else '')
        counts.append(n if assignment[-1] else 0)
        purities.append(round(purity,4) if assignment[-1] else 0.0)
    m['map']=assignment;m['counts']=counts;m['purity']=purities
    m['support']=sum(n>=args.min_support for n in counts)
    return m


def consolidate_modules(modules, names, loc, args):
    """Merge redundant chromosome-only modules (with gene union); dedup ordered windows.

    Same across-species signature + >= overlap of smaller gene set required. A
    union must retain the original chromosome signature and evidence thresholds.
    Coordinate windows only deduplicate if they overlap physically on the SAME
    anchor chromosome: this prevents merging disjoint conserved segments.
    """
    modules.sort(key=lambda m:(-m['support'],-min(x for x in m['counts'] if x),-len(m['genes'])))
    selected=[];by_signature=defaultdict(list);ordered_index=defaultdict(set);report=[]
    for m in modules:
        if not m['genes']:continue
        sig=module_signature(m)
        if m['type']=='chromosome_level':
            possible=by_signature[('chromosome_level',sig)]
        else:
            # Each distinct position-supported interval remains separate.
            # Full inverted index avoids the old first-eight-witness bug.
            possible=set()
            for g in m['genes']:possible.update(ordered_index.get(g,()))
        merged=False
        for i in possible:
            old=selected[i]
            if old['type']!=m['type']:continue
            shared=len(m['genes'] & old['genes'])
            containment=shared/min(len(m['genes']),len(old['genes']))
            if m['type']=='chromosome_level':
                if containment<args.merge_overlap:continue
                union=old['genes'] | m['genes']
                candidate=dict(old);candidate['genes']=union
                recalc_module(candidate,names,loc,args)
                if module_signature(candidate)!=sig or candidate['support']<args.min_species:
                    continue
                before=len(old['genes']);old.update(candidate)
                old['merged_sources']=old.get('merged_sources',1)+m.get('merged_sources',1)
                report.append(['merged',m['anchor'],m['anchor_chr'],old['anchor'],old['anchor_chr'],shared,round(containment,4),len(union),len(union)-before])
                merged=True;break
            else:
                if containment<.8 or old['anchor']!=m['anchor'] or old['anchor_chr']!=m['anchor_chr']:continue
                if not all(not a or not b or a==b for a,b in zip(old['map'],m['map'])):continue
                if old['start']=='' or m['start']=='':continue
                intersection_bp=min(old['end'],m['end'])-max(old['start'],m['start'])
                if intersection_bp<0:continue
                report.append(['duplicate_window',m['anchor'],m['anchor_chr'],old['anchor'],old['anchor_chr'],shared,round(containment,4),len(old['genes']),0])
                merged=True;break
        if merged:continue
        if len(selected)>=args.max_modules:break
        i=len(selected);m['merged_sources']=m.get('merged_sources',1)
        selected.append(m)
        by_signature[(m['type'],sig)].append(i)
        if m['type']=='ordered_window':
            for g in m['genes']:ordered_index[g].add(i)
    # A merged cluster might now match another accepted cluster; iterate.
    # No module may be silently deleted: merge all genes in chromosome-only mode.
    again=True
    while again:
        again=False
        for i in range(len(selected)):
            a=selected[i]
            if a['type']!='chromosome_level':continue
            for j in range(i+1,len(selected)):
                b=selected[j]
                if b['type']!='chromosome_level' or module_signature(a)!=module_signature(b):continue
                shared=len(a['genes']&b['genes'])
                containment=shared/min(len(a['genes']),len(b['genes']))
                if containment<args.merge_overlap:continue
                cand=dict(a);cand['genes']=a['genes']|b['genes']
                recalc_module(cand,names,loc,args)
                if module_signature(cand)!=module_signature(a):continue
                cand['merged_sources']=a['merged_sources']+b['merged_sources']
                report.append(['merge_clusters',b['anchor'],b['anchor_chr'],a['anchor'],a['anchor_chr'],shared,round(containment,4),len(cand['genes']),len(cand['genes'])-len(a['genes'])])
                a.update(cand);selected.pop(j);again=True;break
            if again:break
    selected.sort(key=lambda m:(-m['support'],-min(x for x in m['counts'] if x),-len(m['genes'])))
    return selected,report


def main():
    args,specs,pospaths=read_args()
    names=[n for n,_ in specs];loc={}; coords={}; qc=[]
    for name,path in specs:
        locus,total,dups=read_legacy(path,args.exclude_trna)
        positions={}; ptot=pdup=0
        if name in pospaths:
            positions,ptot,pdup=read_positions(pospaths[name],args.exclude_trna)
            # Coordinate-unique subset only; chromosome assignment from coordinate input.
            locus={g:c for g,(c,s,e) in positions.items()}
        loc[name]=locus;coords[name]=positions
        qc.append([name,path,len(locus),total,dups,len(positions),ptot,pdup,len(set(locus.values())), 'coordinate' if positions else 'chromosome_only'])
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
    # Full-set duplicate detection and chromosome-only consensus merging.
    raw_module_count=len(modules)
    selected,merge_report=consolidate_modules(modules,names,loc,args)
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
    write(str(prefix)+'_merges.tsv',['Action','Source_anchor','Source_chr','Kept_anchor','Kept_chr','Shared_genes','Smaller_set_overlap','Result_gene_count','New_genes_added'],merge_report)
    write(str(prefix)+'_events.tsv',['Event_ID','Module_1','Module_2','Together_species','Apart_species','Together_n','Apart_n','Score','Interpretation'],er)
    write(str(prefix)+'_qc.tsv',['Species','Input','Usable_unique_symbols','Raw_symbols','Duplicate_symbols','Unique_position_loci','Raw_position_symbols','Duplicate_position_symbols','Chromosomes','Analysis_mode'],qc)
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
    cfg={'species':names,'files':dict(specs),'positions':pospaths,'thresholds':{'min_block':args.min_block,'min_species':args.min_species,'min_support':args.min_support,'window':args.window,'purity':args.purity,'merge_overlap':args.merge_overlap},'module_count':len(selected),'raw_module_count':raw_module_count,'merge_actions':len(merge_report),'event_count':len(events),'warnings':['Comparisons based on Gene Symbols, not validated orthology.','Modules are heuristic; cannot prove physical collinearity or event direction.','Without positions, chromosome-level mapping only.','Across species chromosome labels are independent; group labels are not homologous identifiers.','Events are overlapping candidate pairs, not independent historical rearrangements.']}
    with open(str(prefix)+'_manifest.json','w',encoding='utf-8') as f:json.dump(cfg,f,ensure_ascii=False,indent=2)
    print('物种:', ' → '.join(names))
    print('基因位置输入:',len(pospaths),'/',len(names))
    print('候选模块: 原始',raw_module_count,'→ 合并后',len(selected),'合并操作:',len(merge_report),'联合/分离组合:',len(events))
    print('结果:',str(prefix)+'_modules.tsv,',str(prefix)+'_events.tsv,',str(prefix)+'_qc.tsv')
    print('提示: chromosome_only 模式不会用字母排序冒充基因物理顺序。结果都是候选，不等于独立演化事件。')

if __name__=='__main__':main()
