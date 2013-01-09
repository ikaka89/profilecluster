'''
Created on Nov 28, 2012

@author: akshaykakumanu
'''

import os
import subprocess
from subprocess import call
import fileinput
import FormatConversion

def HashGenomeSeqLength(length):
    ''' Given a tab delimited genome-info file (chr_name\tlength) it creates a dic of chr lengths'''
    chrlength={}
    for line in fileinput.input(length):
        line.strip()
        tabs = line.split("\t")
        chrlength[tabs[0]]=tabs[1]
    return chrlength


def OneAlignmentBed(bedfile,loc):
    ''' bed file has no header'''
    ''' Given a bed read alignment file it colapses all the reads to the 5' alignment position and 
        retuns another onebed.bed file'''
    path = os.getcwd()
    filename = FormatConversion.GetFilename(bedfile)
    base = FormatConversion.StripExtension(filename)
    fileout = open(path+"/"+loc+"/"+base+"_onebed.bed","w")
    for line in fileinput.input([bedfile]):
        line=line.rstrip()
        tabs = line.split("\t")
        chromosome = tabs[0]
        strand = tabs[5]
        start=""
        end=""
        if strand == "+":
            start = tabs[1]
            end = str(int(tabs[1])+1)
        if strand  == "-":
            start = str(int(tabs[2])-1)
            end = tabs[2]
        name = tabs[3]
        score = tabs[4]
        output = chromosome+"\t"+start+"\t"+end+"\t"+name+"\t"+score+"\t"+strand+"\n"
        fileout.write(output)
    fileout.close()
    return None

def IntervalBed(bedfile,loc,intervalsize,strand,length):
    ''' Bed file has no header'''
    ''' Actual interval size is 2*intervalsize'''
    path = os.getcwd()
    filename = FormatConversion.GetFilename(bedfile)
    base = FormatConversion.StripExtension(filename)
    chrlength = HashGenomeSeqLength(length)
    fileout = open(path+"/"+loc+"/"+base+"_"+strand+"_interval.bed","w")
    for line in fileinput.input(bedfile):
        line.rstrip()
        tabs = line.split("\t")
        chromosome = tabs[0]
        start = str(int(tabs[1])-intervalsize)
        end = str(int(tabs[1])+intervalsize)
        if int(start) >= 0 and int(end) < chrlength[chromosome]:
            name = tabs[3]
            score = tabs[4]
            strand=strand
            output = chromosome+"\t"+start+"\t"+end+"\t"+name+"\t"+score+"\t"+strand+"\n"
            fileout.write(output)
    fileout.close()
    return None

def RemoveOverlap(bedfile,loc,eventstxt):
    ''' eventstxt file must contain a header'''
    ''' It is better if the bedfile is sorted however it is not required'''
    path = os.getcwd()
    filename = FormatConversion.GetFilename(bedfile)
    base = FormatConversion.StripExtension(filename)
    fileout = open(path+"/"+loc+"/"+base+"_NoOverlapp.bed","w")
    p = subprocess.Popen(["bedtools","intersect","-wo","-a",bedfile,"-b",bedfile],stdout=subprocess.PIPE)
    overlaps = {}
    ''' The following for loop is to construct a dictionary with each event as the key and the 
    values of the dict is a list of all overlapping events'''
    for line in p.stdout.readlines():
        line=line.strip()
        tabs = line.split("\t")
        if int(tabs[12]) < 200:
            overlaps.setdefault(tabs[3],[]).append(tabs[9])
    ipcount = {}
    ''' The following for loop is to construct a dict with events as keys and the Ip count as values'''
    for row in fileinput.input(eventstxt):
        if not fileinput.isfirstline():
            row=row.strip()
            tabs = row.split("\t")
            ipcount[tabs[0]] = float(tabs[1])
    ''' The following for is to remove all the overlapping intervals'''
    for it in overlaps.keys():
        ''' iter over all the overlapping events '''
        if overlaps[it][-1] != 0 and overlaps[it][-1] != 1:
            ''' Making sure the current event was not visited in a direct or a indirect way '''
            for i in overlaps[it]:
                ''' iter over all the events that 'it' overlaps with'''
                ''' Check is 1 if 'it' has the max ip count and 0 if 'it' does not have (Initialized with 1)'''
                check=1
                ''' making sure not itering over 1 and 0 but only over events'''
                if i != 0 and i != 1:
                    if ipcount[it] < ipcount[i]:
                        check=0
            
            if check == 0:
                ''' If 'it' does no have the max ip count'''
                overlaps[it].append(int(0))
            else:
                ''' If it has the maximum ipcount then make it 1 and make all the other events in its valu
                0 irrespective of what it had earlier. Note: these events will never be visited again'''
                overlaps[it].append(int(1))
                for j in overlaps[it]:
                    if j !=0 and j != 1:
                        overlaps[j].append(int(0))
    ''' For loop that prints all the events that have the max ip count among overlapping events ''' 
    for line in fileinput.input(bedfile):
        line.strip()
        tabs = line.split("\t")
        if tabs[3] in overlaps.keys():
            if overlaps[tabs[3]][-1] == 1:
                fileout.write(line)
        else:
            fileout.write(line)
    return None
    
def BinInterval(bedfile,loc,binsize):
    path =os.getcwd()
    filename = FormatConversion.GetFilename(bedfile)
    base = FormatConversion.StripExtension(filename)
    fileout = open(path+"/"+loc+"/"+base+"_binned.bed","w")
    for line in fileinput.input(bedfile):
        line=line.strip()
        tabs = line.split("\t")
        bins = list(range(0,200,binsize))
        for i in bins:
            chromosome = tabs[0]
            start = str(int(tabs[1])+i)
            end = str(int(tabs[1])+i+binsize)
            name = tabs[3]
            score = tabs[4]
            strand = tabs[5]
            output = chromosome+"\t"+start+"\t"+end+"\t"+name+"\t"+score+"\t"+strand+"\n"
            fileout.write(output)
    return None
       
    
