'''
Created on Dec 1, 2012

@author: akshaykakumanu
'''
import os
import subprocess
from subprocess import call
import fileinput
import FormatConversion
import numpy as np
import itertools as itr


def PerformIntersect(bedfilepos,bedfileneg,bambedfile,loc,rep,intervalsize,binsize):
    path = os.getcwd()
    ''' Change the output filename based on the number of events you use to intersect
        eg:- top1000_intersect, top10000_intersect, all_intersect'''
    fileoutpos = open(path+"/"+"temp"+"/"+rep+"_+all_intersect_"+str(intervalsize)+"_"+str(binsize)+".tab","w")
    ppos = subprocess.Popen(["bedtools","intersect","-c","-s","-a",bedfilepos,"-b",bambedfile],stdout=subprocess.PIPE)
    for line in ppos.stdout.readlines():
        fileoutpos.write(line)
    fileoutpos.close()
    fileoutneg = open(path+"/"+"temp"+"/"+rep+"_-all_intersect_"+str(intervalsize)+"_"+str(binsize)+".tab","w")
    pneg = subprocess.Popen(["bedtools","intersect","-c","-s","-a",bedfileneg,"-b",bambedfile],stdout=subprocess.PIPE)
    for linen in pneg.stdout.readlines():
        fileoutneg.write(linen)
    fileoutneg.close()
    return None
        
       
def CalculateMatDis(tabintersectfilepos,tabintersectfileneg,loc,rep,intervalsize,binsize):
    path = os.getcwd()
    filepos = np.genfromtxt(open(tabintersectfilepos),dtype="str",usecols=(3,6),delimiter='\t')
    nelemts = 2*int(intervalsize)/int(binsize)
    nrows = filepos.shape[0]/nelemts
    ncols=(4*intervalsize)/binsize
    matrix = np.arange(nrows*ncols).reshape(nrows,ncols)
    hashevents ={}
    c=itr.count(ncols/2)
    r=itr.count(0)
    for i in xrange(filepos.shape[0]):
        ccurrent=c.next()
        if ccurrent == ncols/2:
            rcurrent=r.next()
            hashevents.update({filepos[i,0]:rcurrent})
            c=itr.count(0)
            ccurrent = c.next()
            matrix[rcurrent,ccurrent]=filepos[i,1]
        else:
            matrix[rcurrent,ccurrent]=filepos[i,1]
        print rcurrent,ccurrent
    print "ppos done"
    fileneg = np.genfromtxt(open(tabintersectfileneg),dtype="str",usecols=(3,6),delimiter='\t')
    ccurrent=-1*ncols/2
    for j in xrange(fileneg.shape[0]):
        ccurrent=ccurrent-1
        if ccurrent == (-1*ncols/2-1):
            ccurrent=-1
        matrix[int(hashevents[fileneg[j,0]]),ccurrent]=fileneg[j,1]
        print hashevents[filepos[j,0]],ccurrent
    np.savetxt(path+"/"+loc+"/"+rep+"_"+str(intervalsize)+"_"+str(binsize)+"_distribution.mat",matrix, fmt="%0.1f", delimiter="\t")
    filehashout=open(path+"/"+loc+"/"+rep+"_"+str(intervalsize)+"_"+str(binsize)+"_matrix_row_annotation.tab","w")
    for j in hashevents.keys():
        filehashout.write(str(j)+"\t"+str(hashevents[j])+"\n")
    return None

  
    
