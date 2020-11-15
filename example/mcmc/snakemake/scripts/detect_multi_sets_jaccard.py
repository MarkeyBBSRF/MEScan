import os
import sys
import copy
from collections import defaultdict
from math import ceil
from itertools import combinations as combos
import time
def load_geneset_count(f,tg_filter,debug = False):
    v = defaultdict(defaultdict)
    count = 0
    fo_file = f+".filtered"
    start = time.time()
    if os.path.isfile(fo_file):
        with open(fo_file,'r') as fi:
            for line in fi:
                count += 1
                tmp = line.strip().split('\t')

                tg = float(tmp[-1])
                if tg < tg_filter:
                    continue
                if debug:
                    if count > 10000:
                        break
                gs = tmp[0].split(' ')
            
                v[tmp[0]]['genes'] = gs
                v[tmp[0]]['tg'] = tg
    else:
        fo = open(f+".filtered",'w')
        with open(f,'r') as fi:
            for line in fi:
                count += 1
                tmp = line.strip().split('\t')
                tg = float(tmp[-1])
                if tg < tg_filter:
                    continue
                fo.write(line)
                gs = tmp[0].split(' ')
                v[tmp[0]]['genes'] = gs
                v[tmp[0]]['tg'] = tg
        fo.close()
    end = time.time()
    return v


def load_validation_sets(f,tg_filter):
    v = defaultdict(defaultdict)
    count = 0
    fo = open(f+".filtered",'w')
    with open(f,'r') as fi:
        for line in fi:
            count += 1
            tmp = line.strip().split('\t')
            tg = float(tmp[-1])
            if tg < tg_filter:
                continue
            fo.write(line)
            gs = tmp[0].split(' ')
            v[tmp[0]]['genes'] = gs
            v[tmp[0]]['tg'] = tg
    fo.close()
    return v




def get_validated_genesets(large_set,small_set):
    start = time.time()
    ssize = len (next(iter(small_set.values()))['genes'])
    keep = []
    keep_cb = set()
    for ls in large_set:

        cb = combos(large_set[ls]['genes'], ssize)
        cb = list(cb)
        save = True
        for i in cb:
            i = ' '.join(i)
            if i not in small_set:
                save = False
                break

        if save:
            keep_cb |= set([' '.join(x) for x in cb]) 
            keep.append(ls)
    end = time.time()
    print("[valid size = %d] get_validated_genesets time = %d"%(ssize+1,(end-start)/60))
    if len(keep)==0:
        print ("[valid size = %d] nothing left!"%(ssize+1))
        return None
    else:
        out = defaultdict()
        print ("[valid size = %d] %d pass validation !"%(ssize+1,len(keep)))
        for i in keep:
            
            out[i] = {"tg":large_set[i]['tg'],'genes':large_set[i]['genes']}
        return {"comb": keep_cb,"valid":out}

    

gs = defaultdict(defaultdict)

for size in range(2,6):
    p = "%s/gs%d.merge.txt" %(sys.argv[2],size)
    if not os.path.isfile(p):
        print("%s is not exist"%(p))
        sys.exit(1)
    gs[size]['raw'] = p




with open(sys.argv[3],'r') as fi:
    for l in fi:
        cols = l.strip().split()
        if int(cols[0]) in gs:
            gs[int(cols[0])]['tg'] = float(cols[1])


sizes = sorted(gs.keys(),reverse=False)
print(sizes)
size = sizes[0]
gs[size]['filtered'] = load_geneset_count(gs[size]['raw'],gs[size]['tg'],debug=False)

for size in sizes[1:]:
    print(size)
    gs[size]['filtered'] = load_geneset_count(gs[size]['raw'],gs[size]['tg'],debug=False)
    small = size - 1
    tmp =  get_validated_genesets(gs[size]['filtered'],gs[small]['filtered'])
    if tmp is not None:
        gs[size]['validated'] = tmp['valid']
        gs[size]['comb_validated'] = tmp['comb']
    else:
        for x in range(size,(sizes[-1]+1)):
            del gs[x]
        break

del gs[2]

cands = defaultdict(float)

iter = 0
print ("Render candiates: ...")

count = 0
while True:

    if len(gs) < 1:
        break
    contin = True
    for index,i in enumerate(list(gs.keys())[:-1]):
        if i+1 != list(gs.keys())[index+1]:
            contin = False
            print ("STOPPED: sm = %d, missing lg = %d" %(i,i+1))
            break
    if not contin:
        break

    # print ("\n--------------Iter: %d--------------" %iter)
    
    
    for small in list(gs.keys())[:-1]:
        print(small)
        large = small + 1
        if gs[large]['comb_validated'] is None:
            print ('Update combinations')
            keep_cb = set()
            for i in gs[large]['validated']:
                cb = combos(gs[large]['validated'][i]['genes'],small)
                cb =list(cb)
            keep_cb |= set([' '.join(x) for x in cb])
            gs[large]['comb_validated'] = keep_cb
            print ("Updated the combinations for size = %d" %(small+1))

        gs[small]['large_clique'] = copy.copy(gs[small]['validated'])
        for gene_set in gs[small]['validated']:

            if gene_set in gs[large]['comb_validated']:

                del gs[small]['large_clique'][gene_set]

    if(len(gs[list(gs.keys())[-1]]['validated'])==0):
      break
    gs[list(gs.keys())[-1]]['large_clique'] = copy.copy(gs[list(gs.keys())[-1]]['validated'])
    lsize = list(gs.keys())[-1]
    large_clique = defaultdict(float)
    tmp = gs[lsize]['large_clique']
    gsid = sorted(tmp.keys(), key=lambda S: tmp[S]["tg"], reverse=True)[0]

    cands[gsid] = tmp[gsid]['tg']
    lc_id = gsid
    lc_tg = tmp[gsid]['tg']
    print("[render] size: %d -> %s tg = %.3f" %(lsize,lc_id,lc_tg))
    print("[render] filtering similar sets ...")
    lc_genes = set(lc_id.split(" "))
    if lsize ==3:
        count += 1
    if count> 100000:
        break

    if lsize  == 3:
        del gs[3]['validated'][lc_id]
        continue
    for size in list(gs.keys()):

        tmp = copy.copy(gs[size]['validated'])

        for gene_set in gs[size]['validated']:
            cg = set(gs[size]['validated'][gene_set]['genes']).intersection(lc_genes)
            common = len(cg)

            if common > 0:
                min_diff = (size-common)*1.0
                dist = min_diff/size
                if size == 3 and dist ==0:
                    del tmp[gene_set]
                elif size >3 and dist <0.5:
                    del tmp[gene_set]

        if len(tmp) == 0:
            tmp = None
            gs[size]['validated'] = None
            print( "[render] size: %d has nothing left!" %size)
            del gs[size]
        else:
            gs[size]['validated'] = tmp
            gs[size]['comb_validated'] = None
            print ("[render] size: = %d, valid = %d " %(size,len(gs[size]['validated'])))

    iter += 1

with open(sys.argv[1], 'w') as ostream:
    for i in cands:
        ostream.write(" ".join([i,str("count"),str(cands[i])]))
        ostream.write("\n")









