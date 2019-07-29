# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 15:18:59 2017

@author: S s L
"""

from wext import wre_test, re_test, SADDLEPOINT, EXACT
import numpy as np
import itertools
import sys
outdir = sys.argv[1]+'/'
simuIndex = sys.argv[2]
########################## test simuIndex = "9"
###reading Observed mutation matrix
ObsMutMatAll = list()
for j in range(1,101):
    print(j)
    tmpiter = str(j)
    ObsMutMatTmp = np.genfromtxt("/home/tzh232/newrunCoMex/Simu"+simuIndex+"ObsMat"+tmpiter+".CSV",delimiter=",")
    ObsMutMatAll.append(ObsMutMatTmp)
#################################
#mutid = np.genfromtxt("mut"+simuIndex+"IdgeneMat.CSV",delimiter=",")-1
gene100idMat = np.genfromtxt("/home/tzh232/newrunCoMex/Gid100Times"+simuIndex+".CSV",delimiter=",")-1
#################
#mutID = list(np.sort(mutid.astype(int)))
mutID = np.array([1,2,3])-1
#####################
rankingVec = []
N = ObsMutMatTmp.shape[1]
Glength = ObsMutMatTmp.shape[0]
for simu100 in range(100):
    #geneId1 = gene100idMat[:,simu100]
    #geneId1 = list(geneId1.astype(int))
    #RatMat = BigRatMat[geneId1,:]
    MutMat = ObsMutMatAll[simu100]  #each time Obs matrix are different
    ############################
    tempComb = itertools.combinations(range(Glength),3)
    AllComb = np.array(list(tempComb))
    CombNum = AllComb.shape[0]
    PvalueVec = []
    for tim in range(CombNum):
        EachId = list(AllComb[tim,:])
        wj = np.sum(MutMat[EachId,:],axis=1)/N
        subRatMat = np.transpose(np.tile(wj,(N,1)))
        subMutMat = MutMat[EachId,:]
        P = subRatMat
        X = np.sum(subMutMat,axis=1) # num of mutation for each gens over all patients
        X = X.astype(int)
        mutornot = np.sum(subMutMat,axis=0) # num of patient
        T = sum(mutornot==1)                        
        Tbvec= []
        for i in [0,1]:
            for j in [0,1]:
                for k in [0,1]:
                    TempV0 = np.tile([i,j,k],(N,1))
                    TempV = np.transpose(TempV0)
                    indexofcat = np.sum(subMutMat - TempV,axis=0)
                    TempNum = np.sum(indexofcat==0)
                    Tbvec = np.append(Tbvec,TempNum)             
        tbl = Tbvec.astype(int)
        PvalueTp = re_test( T, X, tbl, method=SADDLEPOINT )
        PvalueVec = np.append(PvalueVec,PvalueTp) 		
        #############################################################
    mutIDmat = np.tile(mutID,(CombNum,1))
    diffTemp = abs(AllComb-mutID)
    mutpos = np.where(diffTemp.sum(axis=1)==0)
    mutPvalue = PvalueVec[mutpos]
    aftersortP = np.sort(PvalueVec)
    INDEX = sorted(range(len(PvalueVec)),key = lambda k:PvalueVec[k])
    subsetTop = AllComb[INDEX[:100],:]
    subsetTop = np.column_stack([subsetTop,aftersortP[:100]])
    tprank = np.where( mutPvalue == aftersortP)[0]
    rankingTemp = tprank[0]
    print(rankingTemp)
    rankingVec = np.append(rankingVec,rankingTemp)
    np.savetxt(outdir+'genesetstr_test'+simuIndex+str(simu100)+'.csv', subsetTop, delimiter = ',') 



np.savetxt(outdir+'resultstr_test'+simuIndex+'.csv', rankingVec, delimiter = ',') 





