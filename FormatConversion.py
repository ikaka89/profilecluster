'''
Created on Nov 27, 2012

@author: akshaykakumanu
'''
from subprocess import call
import fileinput
import os
import subprocess
import numpy as np

def MakeDirectories():
    ''' makes 3 directories output, temp, log in the directory from which you call this script'''
    path = os.getcwd() 
    call(["mkdir",path+"/output"])
    call(["mkdir",path+"/temp"])
    call(["mkdir",path+"/log"])
    return None

def StripExtension(filename):
    ''' Removes the path of the file and returns only the filename with file type extension'''
    tabs = filename.split(".")
    return tabs[0]

def GetFilename(filename):
    ''' Given a filename it removes the extension and resturns the base name without the extension'''
    tabs = filename.split("/")
    return tabs[-1]
    
    
def EventsToBed(events,loc):
    ''' events file must contain a header'''
    ''' The output Bed file is not a true bed file as the strand tab contains * '''
    
    filename = GetFilename(events)
    base = StripExtension(filename)
    path = os.getcwd()
    fileout = open(path+"/"+loc+"/"+base+".bed","w")
    for line in fileinput.input([events]):
        if not fileinput.isfirstline():
            line.rstrip()
            tabs = line.split("\t")
            coordinates = tabs[0].split(":")
            chromosome="chr"+coordinates[0]
            #chromosome=coordinates[0]
            start=str(int(coordinates[1])-1)
            end=coordinates[1]
            name=tabs[0]
            strand="*"
            output=chromosome+"\t"+start+"\t"+end+"\t"+name+"\t"+"*\t"+strand+"\n"
            fileout.write(output)
    fileout.close
    return None

def SortBam(bamfile,loc):
    ''' Sorts the bam file and adds _sorted tag to the filename'''
    filename = GetFilename(bamfile)
    base = StripExtension(filename)
    path = os.getcwd()
    call(["samtools","sort",bamfile,path+"/"+loc+"/"+base+"_sort"])
    call(["samtools","index",path+"/"+loc+"/"+base+"_sort"+".bam"])
    return None
def SortBed(bedfile,loc):
    ''' Sorts the Bed file and adds _sorted tag to the filename '''
    filename = GetFilename(bedfile)
    base = StripExtension(filename)
    path = os.getcwd()
    outfile =open(path+"/"+loc+"/"+base+"_sorted.bed","w")
    p = subprocess.Popen(["bedtools","sort","-i",bedfile],stdout=subprocess.PIPE)
    for line in p.stdout.readlines():
        outfile.write(line)
    return None
    
def BamToBed(bamfile,loc):
    '''  Converts a bam file to a bed file'''
    filename = GetFilename(bamfile)
    base = StripExtension(filename)
    path = os.getcwd()
    outfile = open(path+"/"+loc+"/"+base+".bed","w")
    p = subprocess.Popen(["bedtools","bamtobed","-i",bamfile],stdout=subprocess.PIPE)
    for line in p.stdout.readlines():
        outfile.write(line)
    return None

def idxtobed(idxfile,loc):
    ''' converts a idx file to a onebed.bed file'''
    path = os.getcwd()
    filename = GetFilename(idxfile)
    base = StripExtension(filename)
    idxtab = np.genfromtxt(open(idxfile),dtype="str",delimiter="\t",usecols=(0,1,2,3))
    fileout = open(path+"/"+loc+"/"+base+"_onebed.bed","w")
    for i in xrange(idxtab.shape[0]):
        chromosome = idxtab[i,0]
        start = idxtab[i,1]
        end = str(int(idxtab[i,1])+1)
        for p in xrange(int(idxtab[i,2])):
            outp = chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+"*"+"\t"+"*"+"\t"+"+"
            fileout.write(outp+"\n")
        for n in xrange(int(idxtab[i,3])):
            outn = chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+"*"+"\t"+"*"+"\t"+"-"
            fileout.write(outn+"\n")
    return None
    
    
    
    