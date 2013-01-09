'''
Created on Dec 17, 2012

@author: akshaykakumanu
'''

import numpy as np
import os

path= os.getcwd()


mergedclus = np.genfromtxt(open("output/clustermerged_rep-1_top10000.mat"),dtype="str",delimiter="\t")
blacklistbed = np.genfromtxt(open("temp/blackist_events_top10000.bed"),dtype="str",usecols=(3), delimiter="\t")
blacklisthash={}
print blacklistbed[0]
for b in xrange(blacklistbed.shape[0]):
    blacklisthash.setdefault(blacklistbed[b],1)
    #print blacklistbed[b]
for r in xrange(mergedclus.shape[0]):
    #mergedclus[r,1:]=[int(i) for i in mergedclus[r,1:]]
    if mergedclus[r,0] in blacklisthash.keys():
        mergedclus[r,1:]="640"
fout = open(path+"/"+"output"+"/"+"clustermerged_"+"rep-1_"+"top10000_"+"blacklist.mat","w")

for o in xrange(mergedclus.shape[0]):
    strout = "\t".join(mergedclus[o,:])
    fout.write(strout+"\n")

        
    
    