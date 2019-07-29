# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 14:45:56 2017

@author: S s L
"""

import sys
import os

import random
import math
import itertools
import numpy as np
##############
def singlemeasure(genes_collection1):
	#coverage of genes in genes_collection1
	out1 = 0
	#total number of mutations in genes_collection1
	inside1 = 0

	for sampleID in sample_mutatedGenes:
		genes_in_sample = sample_mutatedGenes[sampleID]
		inside_genes1 = genes_collection1.intersection(genes_in_sample)
		if len(inside_genes1)>0:
			out1 += 1
		num_ig1 = len(inside_genes1)
		inside1 += num_ig1
	c = 0.5
	return c*float(2*out1 - inside1)


######################
########################## test simuIndex = "9"
simuIndex = "65"
###reading gene matrix and mutated gene matrix
allgenes = np.genfromtxt("test/Gnames100Times"+simuIndex+".CSV",dtype='str',delimiter=",")
allmutgenes = np.genfromtxt("test/mut"+simuIndex+"geneMat.CSV",dtype='str',delimiter=",")
#################################

rankDen=[]
for k in range(100):
    #########################
    ####temp genes and mut genes
    mut_gene = allmutgenes[:,k]
    genes = allgenes[:,k]
    sample_mut_f = open("test/Simu"+simuIndex+"Matation_matrix"+str(k+1)+".txt" ,'r')
    sample_mutatedGenes = dict()
    gene_mutatedSamples = dict()
    print "Loading mutations..."
    all_samples = set()
    ##################
    print "Load genes..."
    for line in list(genes):
    	v =line.split()
    	gene_mutatedSamples[v[0]]=set()
    ##################
    for line in sample_mut_f:
    	v = line.split("\t")
    	sampleID = v[0]
    	all_samples.add( sampleID )
    	sample_mutatedGenes[sampleID]=set()
    	for i in range(len(v) - 1):
    		gene = v[i+1].strip("\n")
    		if gene in genes:
    			if gene not in gene_mutatedSamples:
    				gene_mutatedSamples[gene]=set()
    			gene_mutatedSamples[gene].add(sampleID)
    			sample_mutatedGenes[sampleID].add(gene)
    
    sample_mut_f.close()
        
    #############
    ##############
    #genes is a list of all genes
    allvalue = []
    AllComb = []
    for it in itertools.combinations(genes,3):
        to_set = set(it)
        temp = singlemeasure(to_set)
        allvalue.append(temp)
        AllComb.append(to_set)
        
    ################################################
    INDEX = sorted(range(len(allvalue)),reverse=True,key = lambda k:allvalue[k])
    aftersort=sorted(allvalue,reverse=True)
    c = [AllComb[i] for i in INDEX[:10]]
    mutset = set(mut_gene)
    mutscore = singlemeasure(mutset)
    pos = aftersort.index(mutscore)+1
    rankDen.append(pos)
    print(pos)
    np.savetxt('Denrix'+simuIndex+str(k)+'.csv', c, delimiter = ',',fmt='%s') 
    
    
np.savetxt('Denrix'+simuIndex+'.csv', rankDen, delimiter = ',') 