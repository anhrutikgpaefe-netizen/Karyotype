#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Multi-species linkage with frozen per-parent gene submodules and event-specific evidence.

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
    p.add_argument('--merge-added-purity',type=float,default=.65,help='genes a merge would add must have at least this fraction on the target chromosome in every species where it is defined (default 0.65)')
    p.add_argument('--merge-purity',type=float,default=.50,help='purity floor used only when re-evaluating a merged gene union (default 0.50; never above --purity)')
    p.add_argument('--min-submodule',type=int,default=15,help='minimum size of a fixed chromosome-signature submodule')
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
    if a.min_submodule<2 or a.window<10 or a.min_block<2 or a.min_support<1 or not 0<a.purity<=1 or not 0<a.merge_overlap<=1:p.error('invalid thresholds')
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


def recalc_module(m, names, loc, args, purity=None):
    """Recompute support after union, never reuse original per-species counts.

    purity overrides args.purity (v5: merged unions of near-threshold modules are judged
    with a relaxed floor so a purity dip alone does not block an otherwise conflict-free merge).
    """
    floor=args.purity if purity is None else purity
    gset=m['genes']
    assignment=[];counts=[];purities=[]
    for name in names:
        ct=Counter(loc[name][g] for g in gset if g in loc[name])
        if not ct:
            assignment.append('');counts.append(0);purities.append(0.0);continue
        best,n=sorted(ct.items(),key=lambda x:(-x[1],nk(x[0])))[0]
        purity=n/sum(ct.values())
        assignment.append(best if n>=args.min_support and purity>=floor else '')
        counts.append(n if assignment[-1] else 0)
        purities.append(round(purity,4) if assignment[-1] else 0.0)
    m['map']=assignment;m['counts']=counts;m['purity']=purities
    m['support']=sum(n>=args.min_support for n in counts)
    return m


def compatible_signature(candidate, sig):
    """No species may be assigned different chromosomes; a blank on either side is not a conflict."""
    return all(not a or not b or a==b for a,b in zip(candidate,sig))


def merge_union(base, other_genes, sig, names, loc, args):
    """Union of gene sets re-evaluated with the relaxed merge purity.

    Returns (candidate, relaxed) or None. relaxed is True when the strict --purity would have
    changed the signature, i.e. the merge only passes thanks to --merge-purity.
    """
    # Genes that would be ADDED must themselves sit on the target chromosomes, otherwise a
    # mixed module (e.g. one chromosome carrying two conserved blocks) would drag a foreign
    # block into this one. Judged per species where the target signature is defined.
    added=other_genes-base['genes']
    for name,chrom in zip(names,sig):
        if not chrom:continue
        ct=Counter(loc[name][g] for g in added if g in loc[name])
        n=sum(ct.values())
        if n>=args.min_support and ct.get(chrom,0)/n<args.merge_added_purity:return None
    candidate=dict(base);candidate['genes']=base['genes']|other_genes
    floor=min(args.merge_purity,args.purity)
    recalc_module(candidate,names,loc,args,purity=floor)
    if not compatible_signature(candidate['map'],sig) or candidate['support']<args.min_species:
        return None
    strict=dict(candidate);recalc_module(strict,names,loc,args)
    return candidate,module_signature(strict)!=sig


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
                res=merge_union(old,m['genes'],sig,names,loc,args)
                if res is None:continue
                candidate,relaxed=res;union=candidate['genes']
                before=len(old['genes']);old.update(candidate)
                old['merged_sources']=old.get('merged_sources',1)+m.get('merged_sources',1)
                old['relaxed_merge']=old.get('relaxed_merge',False) or relaxed
                report.append(['merged',m['anchor'],m['anchor_chr'],old['anchor'],old['anchor_chr'],shared,round(containment,4),len(union),len(union)-before,int(relaxed)])
                merged=True;break
            else:
                if containment<.8 or old['anchor']!=m['anchor'] or old['anchor_chr']!=m['anchor_chr']:continue
                if not all(not a or not b or a==b for a,b in zip(old['map'],m['map'])):continue
                if old['start']=='' or m['start']=='':continue
                intersection_bp=min(old['end'],m['end'])-max(old['start'],m['start'])
                if intersection_bp<0:continue
                report.append(['duplicate_window',m['anchor'],m['anchor_chr'],old['anchor'],old['anchor_chr'],shared,round(containment,4),len(old['genes']),0,0])
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
                res=merge_union(a,b['genes'],module_signature(a),names,loc,args)
                if res is None:continue
                cand,relaxed=res
                cand['merged_sources']=a['merged_sources']+b['merged_sources']
                cand['relaxed_merge']=a.get('relaxed_merge',False) or b.get('relaxed_merge',False) or relaxed
                report.append(['merge_clusters',b['anchor'],b['anchor_chr'],a['anchor'],a['anchor_chr'],shared,round(containment,4),len(cand['genes']),len(cand['genes'])-len(a['genes']),int(relaxed)])
                a.update(cand);selected.pop(j);again=True;break
            if again:break
    selected.sort(key=lambda m:(-m['support'],-min(x for x in m['counts'] if x),-len(m['genes'])))
    return selected,report


def write_network_reports(prefix, selected, names, loc, args, events):
    """Separate within-module dispersion from between-module co-assignment.

    These are *chromosome membership* relationships, not physical adjacency.
    """
    member_heads=['Module_ID','Gene_count','Mapped_species','Species','Distinct_chromosomes',
                  'Main_chr','Main_count','Mapped_gene_count','Main_fraction',
                  'Other_chr_counts','Absent_gene_count','Status']
    integrity=[]
    for i,m in enumerate(selected,1):
        mid=f'M{i:05}'
        for name in names:
            ct=Counter(loc[name][g] for g in m['genes'] if g in loc[name])
            present=sum(ct.values())
            ordered=sorted(ct.items(),key=lambda z:(-z[1],nk(z[0])))
            main,count=ordered[0] if ordered else ('',0)
            fraction=round(count/present,4) if present else 0.0
            if not present:status='NO_SYMBOL_MATCH'
            elif count<args.min_support:status='LOW_SUPPORT'
            elif fraction<args.purity:status='DISPERSED'
            elif len(ordered)>1:status='DOMINANT_WITH_MINOR_LOCATIONS'
            else:status='SINGLE_CHROMOSOME'
            integrity.append([mid,len(m['genes']),m['support'],name,len(ordered),main,count,present,
                              fraction,'|'.join(f'{c}:{n}' for c,n in ordered[1:]),
                              len(m['genes'])-present,status])
    write(str(prefix)+'_module_integrity.tsv',member_heads,integrity)

    # Every module pair is kept; shared hub membership is not a duplicate.
    edges=[]; per_module=defaultdict(list); chrom_members=defaultdict(list)
    for i,m in enumerate(selected):
        for k,name in enumerate(names):
            if m['map'][k]:chrom_members[(name,m['map'][k])].append((i,m['counts'][k]))
    for i,m in enumerate(selected):
        for j in range(i+1,len(selected)):
            q=selected[j]
            together=[];apart=[];missing=[]
            for k,name in enumerate(names):
                a,b=m['map'][k],q['map'][k]
                if not a or not b:missing.append(name)
                elif a==b:together.append(name)
                else:apart.append(name)
            if not together:continue
            mid1,mid2=f'M{i+1:05}',f'M{j+1:05}'
            # Cross-species co-assignment, not proven physical neighborhood.
            edges.append([mid1,mid2,','.join(together),','.join(apart),
                          ','.join(missing),len(together),len(apart),
                          round(len(together)/(len(together)+len(apart)),4) if together or apart else 0,
                          int(bool(apart)),
                          round(len(m['genes']&q['genes'])/min(len(m['genes']),len(q['genes'])),4)])
            per_module[i].append((j,together,apart))
            per_module[j].append((i,together,apart))
    write(str(prefix)+'_module_links.tsv',
          ['Module_1','Module_2','Same_chr_species','Different_chr_species','Missing_species',
           'Same_n','Different_n','Same_fraction_observed','Variable_link','Gene_set_overlap_smaller'],edges)
    neighbors=[]
    for i,m in enumerate(selected):
        ls=per_module[i]
        associates=sorted(ls,key=lambda z:(-len(z[1]),-len(z[2]),z[0]))
        neighbors.append([f'M{i+1:05}',len(m['genes']),len(ls),
                          sum(bool(t) and bool(a) for _,t,a in ls),
                          '|'.join(f'M{j+1:05}' for j,_,_ in associates),
                          '|'.join(f'M{j+1:05}:{len(t)}same/{len(a)}different' for j,t,a in associates)])
    write(str(prefix)+'_module_neighbors.tsv',
          ['Module_ID','Gene_count','Linked_module_count','Variable_link_count',
           'Linked_module_IDs','Link_profiles'],neighbors)
    composition=[]
    for (name,chrom),v in sorted(chrom_members.items(),key=lambda x:(names.index(x[0][0]),nk(x[0][1]))):
        v.sort(key=lambda z:(-z[1],z[0]))
        composition.append([name,chrom,len(v),'|'.join(f'M{i+1:05}' for i,n in v),
                            '|'.join(str(n) for i,n in v),
                            'chromosome_membership_not_physical_order'])
    write(str(prefix)+'_chromosome_composition.tsv',
          ['Species','Chromosome','Module_count','Modules','Anchor_symbol_counts','Interpretation'],composition)
    # Group original pairwise events by their exact tested same/apart pattern,
    # without deleting any member of the group or assuming identical history.
    grouped=defaultdict(list)
    for score,i,j,to,apart,tc,ac in events:
        grouped[(to,apart)].append((i,j,score))
    erows=[]
    for n,((to,apart),rows) in enumerate(sorted(grouped.items(),key=lambda t:(-len(t[1]),t[0])),1):
        erows.append([f'G{n:05}',to,apart,len(rows),
                      '|'.join(f'M{i+1:05}+M{j+1:05}' for i,j,_ in rows),
                      max(round(sc,2) for _,_,sc in rows)])
    write(str(prefix)+'_event_groups.tsv',
          ['Group_ID','Same_chr_species','Different_chr_species','Pair_count','Module_pairs','Max_score'],erows)
    return len(edges),len(grouped)


def build_submodules(selected, names, loc, args):
    """Partition each parental module into disjoint, stable chromosome-signature gene sets.

    Missing symbols are treated as unknown, not a distinct chromosome. A gene may
    belong to only one child of a given parent. The same child ID always denotes
    the same frozen Symbol set across all downstream events.
    """
    result=[]; metadata=[]; assignments=[]
    for i,parent in enumerate(selected,1):
        pid=f'M{i:05}'
        # Seed patterns with all available species, including missing. Group exact signatures first.
        buckets=defaultdict(set)
        for gene in parent['genes']:
            signature=tuple(loc[n].get(gene,'') for n in names)
            buckets[signature].add(gene)
        # A greedy consensus grouping: merge compatible buckets when no species contradicts.
        # Use larger, well-observed patterns first; never merge conflicting chromosomes.
        seeds=sorted(buckets.items(),key=lambda kv:(-len(kv[1]),-sum(bool(x) for x in kv[0]),kv[0]))
        clusters=[]
        for signature,genes in seeds:
            candidates=[]
            for ci,(consensus,gset) in enumerate(clusters):
                conflict=sum(bool(a and b and a!=b) for a,b in zip(consensus,signature))
                matches=sum(bool(a and b and a==b) for a,b in zip(consensus,signature))
                if conflict==0 and matches>=max(2,min(args.min_species-1,3)):
                    candidates.append((matches,len(gset),-ci,ci))
            if candidates:
                ci=max(candidates)[-1];consensus,gset=clusters[ci]
                clusters[ci]=(tuple(a or b for a,b in zip(consensus,signature)),gset|genes)
            else:clusters.append((signature,set(genes)))
        # Absorb undersized clusters only if compatible, otherwise retain as remainder.
        min_size=args.min_submodule
        large=[(sig,set(gs)) for sig,gs in clusters if len(gs)>=min_size]
        small=[(sig,set(gs)) for sig,gs in clusters if len(gs)<min_size]
        remainder=set()
        for signature,genes in small:
            candidates=[]
            for ci,(consensus,gset) in enumerate(large):
                conflict=sum(bool(a and b and a!=b) for a,b in zip(consensus,signature))
                matches=sum(bool(a and b and a==b) for a,b in zip(consensus,signature))
                if conflict==0 and matches>=2:candidates.append((matches,len(gset),-ci,ci))
            if candidates:
                ci=max(candidates)[-1];consensus,gset=large[ci]
                large[ci]=(tuple(a or b for a,b in zip(consensus,signature)),gset|genes)
            else:remainder.update(genes)
        # If the parent is globally stable, preserving a single complete child is desirable.
        if not large:
            large=[(tuple(parent['map']),set(parent['genes']))]
            remainder=set()
        elif remainder and len(remainder)>=min_size:
            large.append((tuple('' for _ in names),remainder))
            remainder=set()
        # Freeze disjoint child sets; unresolved remainder is explicitly exported but not in events.
        large.sort(key=lambda z:(-len(z[1]),z[0]))
        for j,(sig,genes) in enumerate(large,1):
            child=dict(parent)
            child['genes']=frozenset(genes)
            recalc_module(child,names,loc,args)
            sid=f'{pid}-{j}'
            child['id']=sid;child['parent_id']=pid
            child['eligible']=len(genes)>=min_size and child['support']>=args.min_species
            result.append(child)
            metadata.append([sid,pid,j,len(genes),child['support'],int(child['eligible']),
                             'fixed_symbol_set_chromosome_pattern','|'.join(sorted(genes))])
            for gene in sorted(genes):assignments.append([sid,pid,gene])
        if remainder:
            for gene in sorted(remainder):assignments.append(['',pid,gene])
    return result,metadata,assignments


def write_submodule_reports(prefix, parents, children, metadata, assignments, names, loc, args):
    write(str(prefix)+'_submodules.tsv',
          ['Submodule_ID','Parent_Module_ID','Index','Gene_count','Supported_species',
           'Event_eligible','Definition','Gene_symbols'],metadata)
    write(str(prefix)+'_submodule_genes.tsv',['Submodule_ID','Parent_Module_ID','Gene_symbol'],assignments)
    chromrows=[]
    for ch in children:
        for n in names:
            ct=Counter(loc[n][g] for g in ch['genes'] if g in loc[n])
            ordered=sorted(ct.items(),key=lambda z:(-z[1],nk(z[0])))
            total=sum(ct.values()); dominant=ordered[0][0] if ordered else ''
            dominant_count=ordered[0][1] if ordered else 0
            chromrows.append([ch['id'],ch['parent_id'],n,dominant,dominant_count,total,
                             round(dominant_count/total,4) if total else 0,
                             '|'.join(f'{c}:{v}' for c,v in ordered),
                             int(bool(ch['map'][names.index(n)]))])
    write(str(prefix)+'_submodule_chromosomes.tsv',
          ['Submodule_ID','Parent_Module_ID','Species','Main_chromosome',
           'Main_gene_count','Mapped_gene_count','Dominance','All_chromosome_counts',
           'Pass_module_support'],chromrows)
    eligible=[m for m in children if m['eligible']]
    events=[]
    for i,a in enumerate(eligible):
        for b in eligible[i+1:]:
            if a['parent_id']==b['parent_id']:continue
            if len(a['genes']&b['genes'])>min(len(a['genes']),len(b['genes']))*.20:continue
            together=[n for k,n in enumerate(names) if a['map'][k] and a['map'][k]==b['map'][k]]
            apart=[n for k,n in enumerate(names) if a['map'][k] and b['map'][k] and a['map'][k]!=b['map'][k]]
            if not together or not apart or len(together)+len(apart)<args.min_species:continue
            score=round(min(sum(a['counts']),sum(b['counts']))*(len(together)+len(apart))/len(names),2)
            events.append((score,a,b,together,apart))
    events.sort(key=lambda e:(-e[0],e[1]['id'],e[2]['id']))
    er=[];eg=[]
    for idx,(score,a,b,together,apart) in enumerate(events,1):
        eid=f'E{idx:05}'
        er.append([eid,a['id'],b['id'],a['parent_id'],b['parent_id'],
                   ','.join(together),','.join(apart),len(together),len(apart),score,'candidate_only'])
        # Event evidence is per child, not union of both (their overlap is not presumed).
        for child in (a,b):
            for gene in sorted(child['genes']):eg.append([eid,child['id'],gene])
    write(str(prefix)+'_submodule_events.tsv',
          ['Event_ID','Submodule_1','Submodule_2','Parent_1','Parent_2',
           'Together_species','Apart_species','Together_n','Apart_n','Score','Interpretation'],er)
    write(str(prefix)+'_event_genes.tsv',['Event_ID','Submodule_ID','Gene_symbol'],eg)
    return len(eligible),len(events)


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
    head+=['Anchor_gene_symbols','Merged_sources','Relaxed_merge']
    mr=[]
    for i,m in enumerate(selected,1):
        row=[f'M{i:05}',m['anchor'],m['anchor_chr'],m['start'],m['end'],len(m['genes']),m['support'],m['type']]
        for c,num,p in zip(m['map'],m['counts'],m['purity']):row.extend((c,num,p))
        row+=['|'.join(sorted(m['genes'])),m.get('merged_sources',1),int(m.get('relaxed_merge',False))];mr.append(row)
    er=[]
    for k,(score,i,j,to,apart,tc,ac) in enumerate(events,1):
        er.append([f'E{k:05}',f'M{i+1:05}',f'M{j+1:05}',to,apart,tc,ac,score,'candidate_only'])
    write(str(prefix)+'_modules.tsv',head,mr)
    write(str(prefix)+'_merges.tsv',['Action','Source_anchor','Source_chr','Kept_anchor','Kept_chr','Shared_genes','Smaller_set_overlap','Result_gene_count','New_genes_added','Relaxed_purity'],merge_report)
    write(str(prefix)+'_events.tsv',['Event_ID','Module_1','Module_2','Together_species','Apart_species','Together_n','Apart_n','Score','Interpretation'],er)
    link_count,pattern_count=write_network_reports(prefix,selected,names,loc,args,events)
    children,submeta,subassign=build_submodules(selected,names,loc,args)
    eligible_children,subevents=write_submodule_reports(prefix,selected,children,submeta,subassign,names,loc,args)
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
    cfg={'species':names,'files':dict(specs),'positions':pospaths,'thresholds':{'min_block':args.min_block,'min_species':args.min_species,'min_support':args.min_support,'window':args.window,'purity':args.purity,'merge_overlap':args.merge_overlap,'merge_purity':args.merge_purity,'merge_added_purity':args.merge_added_purity},'module_count':len(selected),'raw_module_count':raw_module_count,'merge_actions':len(merge_report),'event_count':len(events),'module_link_count':link_count,'event_pattern_count':pattern_count,'submodule_count':len(children),'eligible_submodule_count':eligible_children,'submodule_event_count':subevents,'warnings':['Comparisons based on Gene Symbols, not validated orthology.','Modules are heuristic; cannot prove physical collinearity or event direction.','Without positions, chromosome-level mapping only.','Across species chromosome labels are independent; group labels are not homologous identifiers.','Events are overlapping candidate pairs, not independent historical rearrangements.']}
    with open(str(prefix)+'_manifest.json','w',encoding='utf-8') as f:json.dump(cfg,f,ensure_ascii=False,indent=2)
    print('物种:', ' → '.join(names))
    print('基因位置输入:',len(pospaths),'/',len(names))
    print('候选模块: 原始',raw_module_count,'→ 合并后',len(selected),'合并操作:',len(merge_report),'联合/分离组合:',len(events))
    print('结果:',str(prefix)+'_modules.tsv,',str(prefix)+'_events.tsv,',str(prefix)+'_qc.tsv')
    print('模块完整性:',str(prefix)+'_module_integrity.tsv','模块连接:',str(prefix)+'_module_links.tsv')
    print('模块邻居:',str(prefix)+'_module_neighbors.tsv','染色体组成:',str(prefix)+'_chromosome_composition.tsv')
    print('事件模式:',str(prefix)+'_event_groups.tsv', '模块连接数:',link_count,'事件模式数:',pattern_count)
    print('固定子模块:',len(children),'符合事件条件:',eligible_children,'子模块事件:',subevents)
    print('子模块:',str(prefix)+'_submodules.tsv','事件:',str(prefix)+'_submodule_events.tsv','基因:',str(prefix)+'_event_genes.tsv')
    print('提示: chromosome_only 模式不会用字母排序冒充基因物理顺序。结果都是候选，不等于独立演化事件。')

if __name__=='__main__':main()
