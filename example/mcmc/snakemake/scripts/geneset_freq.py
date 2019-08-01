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
def getGenesetFreq(MCMC_MY_File, burning, outdir):
	fn = re.sub("raw$",'freq', os.path.split(MCMC_MY_File)[-1])
	fo = os.path.join(outdir,fn)
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
	




def get_geneset_freq_wrapper(args):
	return getGenesetFreq(*args)

def run_getGenesetFreqParallel(root_dir,burning,npool=10,size=2):
	BURNING = burning
	MCMC_MY_Files = glob.glob(root_dir+"/raw/gs"+str(size)+".s*.raw")
	outdir = root_dir+"/freq/"
	combdir = root_dir+"/comb/"
	try: 
		os.makedirs(outdir)
	except OSError:
		if not os.path.isdir(outdir):
			raise
	try: 
		os.makedirs(combdir)
	except OSError:
		if not os.path.isdir(combdir):
			raise

	from multiprocessing import Pool
	p = Pool(len(MCMC_MY_Files))
	p.map(get_geneset_freq_wrapper,izip(MCMC_MY_Files,repeat(int(burning)),repeat(outdir)))




if __name__ == '__main__':
	if len(sys.argv)<3:
		print ("python geneset_freq.py root-dir burning npool\n")
		sys.exit(0)
	run_getGenesetFreqParallel(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
