'''
Created on Dec 19, 2012

@author: akshaykakumanu
'''

import CalculateDistribution
import Clustering
import numpy as np
import os

def main():
    ''' importing the stored numpy read distribution array'''
    distrimat = np.genfromtxt(open("output/3hr_run1_100_5_distribution.mat"),delimiter="\t")
    annotation= np.genfromtxt(open("output/3hr_run1_100_5_matrix_row_annotation.tab"),dtype="str",delimiter="\t")
    activeclusters=Clustering.CreateActiveArray(distrimat)
    hashevents={}
    for j in xrange(annotation.shape[0]):
        hashevents.setdefault(annotation[j,1],annotation[j,0])
    stack=Clustering.InitializeStack(activeclusters)
    fromclus=stack[-1]
    minbic=100000000.0
    minbicclus=[]
    path = os.getcwd()
    fileout = open(path+"/"+"output"+"/"+"3hr_run1"+"_100"+"_"+"5_"+"clusters.tab","w")
    while len(activeclusters)>1:
        (frompos,neighpos)=Clustering.NearestNeighbor(distrimat, activeclusters, fromclus)
        (stack, activeclusters)=Clustering.MergeCluster(frompos, neighpos, activeclusters, stack)
        fromclus=stack[-1]
        bic=Clustering.BIC(activeclusters, distrimat)
        print len(activeclusters),bic,len(minbicclus),minbic
        
        #if len(activeclusters) == 46: 
        #if int(bic) < int(minbic) and len(activeclusters)<600:
        if len(activeclusters) == 5:  
            minbicclus=[]
            for i in xrange(len(activeclusters)):
                minbicclus.extend([activeclusters[i]])
            minbic=bic
    print len(minbicclus),"no of final clusters"
    for c in minbicclus:
        order=Clustering.Postorder(c)
        fileout.write(str(c.name)+"\n")
        for event in order:
            strorder=[str(i) for i in distrimat[event]]
            out2="\t".join(strorder)
            output=str(hashevents[str(event)])+"\t"+out2
            fileout.write(output+"\n")
    print minbic,"minimum bic"
    print len(fromclus.nodes)
    
        
if __name__=="__main__":
    main()
