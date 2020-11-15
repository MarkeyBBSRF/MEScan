import os
import re
import sys
import glob
import operator
from collections import defaultdict
from itertools import repeat
try:
    from itertools import izip
except ImportError:
    izip = zip
def getGenesetFreq(MCMC_MY_File, burning, fo):
    #if not os.path.isfile(fo):
    if True:
        results = dict()
        with open(MCMC_MY_File,'r') as lines:
            lc = 0
            for line in lines:
                lc += 1
                if(lc < burning):
                    continue
                cols = line.strip().split("\t")
                geneset = cols[1]
                tg_score = cols[0]
                if geneset in results: 
                    results[geneset]['freq'] += 1
                    results[geneset]['tg'] = float(tg_score)
                else:
                    results[geneset] = dict(freq=1,tg=float(tg_score))
                if lc % 1000000 == 0:
                    print ("Processed %d lines\n"%(lc))

        # results_sorted = sorted(results.items(), key=operator.itemgetter(1),reverse=True)
        key_sorted = sorted(results.keys(), key=lambda S: results[S]["tg"], reverse=True)
        with open(fo,'w') as foutput:
            for i in  key_sorted:
                tg = results[i]['tg']
                if tg < 0:
                    continue
                freq = results[i]['freq']
                foutput.write("\t".join(map(str,[i,freq,tg])))
                foutput.write('\n')


if __name__ == '__main__':
    if len(sys.argv)<2:
        print ("python geneset_freq_single.py input burning output\n")
        sys.exit(0)
    getGenesetFreq(sys.argv[1],int(sys.argv[2]),sys.argv[3])