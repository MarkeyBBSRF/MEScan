import os
import sys
import copy
from collections import OrderedDict
from collections import defaultdict 
from math import ceil
from itertools import combinations as combos
import time
class MyDict(OrderedDict):
    def __missing__(self, key):
        val = self[key] = MyDict()
        return val
def load_geneset_count(f,tg_filter,debug = False):
    v = MyDict()
    count = 0
    fo_file = f+".filtered"
    start = time.time()
    if os.path.isfile(fo_file):
    #print("%s file exist"%(fo_file))
        with open(fo_file,'r') as fi:
            for line in fi:
                count += 1
                tmp = line.strip().split('\t')

                tg = float(tmp[-1])
                if tg < tg_filter:
                    continue
                if debug:
                    if count > 10000000:
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
    #print("load_validation_sets time = %d"%((start-end)/60))
    return v


def first(s):
    '''Return the first element from an ordered collection
       or an arbitrary element from an unordered collection.
       Raise StopIteration if the collection is empty.
    '''
    return next(iter(s))

def get_validated_genesets(large_set,small_set):
    """
    Using size n gs to validate n+1 gs.
    1. get all combination in large gs
    2. if all combinations in size_n, keep it, otherwise through this size_n+1
    3. cb contains all combinations, if keep size_n+1, cb will be saved otherwise not.
    """
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
            # save all the possible combination to a set
            keep_cb |= set([' '.join(x) for x in cb]) 
            keep.append(ls)
    end = time.time()
    print("[valid size = %d] get_validated_genesets time = %d"%(ssize+1,(end-start)/60))
    if len(keep)==0:
        print ("[valid size = %d] nothing left!"%(ssize+1))
        return None
    else:
        out = MyDict()
        print ("[valid size = %d] %d pass validation !"%(ssize+1,len(keep)))
        for i in keep:
            
            out[i] = {"tg":large_set[i]['tg'],'genes':large_set[i]['genes']}
            #sorted(out.keys(), key=lambda S: out[S]["tg"], reverse=True)
        return {"comb": keep_cb,"valid":out}

    

gs = defaultdict(defaultdict)


for size in range(2,8):
    p = "%s/gs%d.merge.txt" %(sys.argv[2],size)
    if not os.path.isfile(p):
        print("%s is not exist"%(p))
        sys.exit(1)
    gs[size]['raw'] = p
    print(p)


# add file names to directory
# for size in range(2,11):
#     p = "%s/gs%d.merge.txt" %(sys.argv[2],size)
#     if not os.path.isfile(p):
#         sys.exit(1)
#     gs[size]['raw'] = p

# load tg filter
with open(sys.argv[3],'r') as fi:
    for l in fi:
        cols = l.strip().split()
        if int(cols[0]) in gs:
            gs[int(cols[0])]['tg'] = float(cols[1])


# filter geneset by tgscore
# for size in gs:
#     gs[size]['filtered'] = load_geneset_count(gs[size]['raw'],gs[size]['tg'],debug=True)

# load the gene set while validating.
# if failed at larger size
# stop there
## load the smallest
sizes = sorted(gs.keys(),reverse=False)
print(sizes)
size = sizes[0]
gs[size]['filtered'] = load_geneset_count(gs[size]['raw'],gs[size]['tg'],debug=True)

# validate larger set by smaller ones
for size in sizes[1:]:
    gs[size]['filtered'] = load_geneset_count(gs[size]['raw'],gs[size]['tg'],debug=True)
    small = size - 1
    tmp =  get_validated_genesets(gs[size]['filtered'],gs[small]['filtered'])
    if tmp is not None:
        gs[size]['validated'] = tmp['valid']
        gs[size]['comb_validated'] = tmp['comb']
    else:
        #print ("Validate %d by %d: nothing left" %(size,size - 1))
        for x in range(size,(sizes[-1]+1)):
            del gs[x]
        break



del gs[2]

cands = MyDict()
iter = 0
print ("Render candiates: ...")
print ("Prefer ones with high TG score")

MAXCAND = 100

for i in gs:
    print("size = ",i)
    cn = 0
    skip = list()
    tmp = gs[i]['validated']
    save=False
    for j in list(tmp.keys()):
        if cn > MAXCAND-1:
            break
        print(j)
        if i != 3:
            if cn==0:
                cands[i][j] = tmp[j]['tg']
                cn+=1
                continue
            for key in cands[i]:
                save = True
                ks = tmp[key]['genes']
                cs = tmp[j]['genes']
                cg = set(ks).intersection(set(cs))
                common = len(cg)
                #print(set(gs[size]['validated'][gene_set]['genes']))
                #print(lc_genes)
                if common > 0:
                    # add jaccard distance calculation
                    min_diff = (i-common)*1.0
                    dist = min_diff/i
                    if dist < 0.5:
                        #print (("         %s dist = %.2f"%(" ".join(cs),dist)))
                        save=False
                        break

            if save:
                cands[i][j] = tmp[j]['tg']
                cn+=1
        else:

            cands[i][j] = tmp[j]['tg']
            cn+=1

with open(sys.argv[1], 'w') as ostream:
    for i in cands:
        for j in cands[i]:
            ostream.write(" ".join([j,str("count"),str(cands[i][j])]))
            ostream.write("\n")









