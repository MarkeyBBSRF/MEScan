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

def mergeFreq(root_dir,minCutoff,size):
	minCutoff = int(minCutoff)
	MCMC_MY_Files = glob.glob(root_dir+"/freq/gs"+str(size)+".s*.freq")
	combdir = root_dir+"/comb/"
	print(MCMC_MY_Files)
	try: 
		os.makedirs(combdir)
	except OSError:
		if not os.path.isdir(combdir):
			raise
	fn = re.sub("\.s.*\.freq$",'.merge.txt', os.path.split(MCMC_MY_Files[0])[-1])
	fn =  os.path.join(combdir,fn)
	results = dict()
	for f in MCMC_MY_Files:
		with open(f,'r') as lines:
			for line in lines:
				cols = line.strip().split("\t")
				geneset = cols[0]
				counts = int(cols[1])
				tg_score = float(cols[2])
				if tg_score < minCutoff:
					continue
				if geneset in results: 
					results[geneset]['freq'] += counts
					
				else:
					results[geneset] = dict(freq=counts,tg=float(tg_score))

	key_sorted = sorted(results.keys(), key=lambda S: results[S]["tg"], reverse=True)
	with open(fn,'w') as foutput:
		for i in  key_sorted:
			tg = results[i]['tg']
			freq = results[i]['freq']
			foutput.write("\t".join(map(str,[i,freq,tg])))
			foutput.write('\n')


if __name__ == '__main__':
	if len(sys.argv)<3:
		print ("python mergeFreq.py root-dir minCutoff gize\n")
		print(len(sys.argv))
		sys.exit(0)
	mergeFreq(sys.argv[1],sys.argv[2],sys.argv[3])
