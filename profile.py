'''
Created on Nov 27, 2012

@author: akshaykakumanu
'''

import FormatConversion
import IntervalSettings
import os
import CalculateDistribution
import time

def main():
    #FormatConversion.MakeDirectories()
    start = time.time()
    print ''' Converting events.txt into a bed file and sorting it'''
    FormatConversion.EventsToBed("input/Gata1_inputCtrl-rep1_Gata1_3h.events","temp")
    FormatConversion.SortBed("temp/Gata1_inputCtrl-rep1_Gata1_3h.bed", "temp")
    end1 = time.time()
    print "Time elapsed: "+str(end1-start)
    #print ''' Sorting and converting BAM alignment files into BED files '''
    #FormatConversion.SortBam("input/Hela_1_bowtie_unique+colorspace_hg18.bam", "temp")
    #FormatConversion.SortBam("input/Hela_2_bowtie_unique+colorspace_hg18.bam", "temp")
    ##FormatConversion.SortBam("input/Hela_3_bowtie_unique+colorspace_hg18.bam", "temp")
    #FormatConversion.BamToBed("temp/Hela_1_bowtie_unique+colorspace_hg18_sorted.bam", "temp")
    #FormatConversion.BamToBed("temp/Hela_2_bowtie_unique+colorspace_hg18_sorted.bam", "temp")
    #FormatConversion.BamToBed("temp/Hela_3_bowtie_unique+colorspace_hg18_sorted.bam", "temp")
    print ''' Sorting and converting idx files into BED files '''
    FormatConversion.idxtobed("input/GATA1_3_ER4_3h_run1.idx", "temp")
    FormatConversion.SortBed("temp/GATA1_3_ER4_3h_run1_onebed.bed", "temp")
    end2=time.time()
    print "Time elapsed: "+str(end2-end1)
    print ''' Converting Bam files reads location to a single point (5' position) as a bed file and then sorting it'''
    #IntervalSettings.OneAlignmentBed("temp/p53_6h_merged-U2OS_sorted.bed", "temp")
    #IntervalSettings.OneAlignmentBed("temp/Hela_2_bowtie_unique+colorspace_hg18_sorted.bed", "temp")
    #IntervalSettings.OneAlignmentBed("temp/Hela_3_bowtie_unique+colorspace_hg18_sorted.bed", "temp")
    #FormatConversion.SortBed("temp/p53_6h_merged-U2OS_sorted_onebed.bed","temp")
    #FormatConversion.SortBed("temp/Hela_2_bowtie_unique+colorspace_hg18_sorted_onbed.bed","temp")
    #FormatConversion.SortBed("temp/Hela_3_bowtie_unique+colorspace_hg18_sorted_onbed.bed","temp")
    end3=time.time()
    print "Time elapsed: "+str(end3-end2)
    print ''' Extend the event location by 100 bp (left and right) in both positive and negative strands'''
    IntervalSettings.IntervalBed("temp/Gata1_inputCtrl-rep1_Gata1_3h_sorted.bed", "temp", 100, "+", "input/mm9.info")
    IntervalSettings.IntervalBed("temp/Gata1_inputCtrl-rep1_Gata1_3h_sorted.bed", "temp", 100, "-", "input/mm9.info")
    end4=time.time()
    print "Time elapsed: "+str(end4-end3)
    print ''' Remove overlapping intervals '''
    IntervalSettings.RemoveOverlap("temp/Gata1_inputCtrl-rep1_Gata1_3h_sorted_+_interval.bed", "temp", "input/Gata1_inputCtrl-rep1_Gata1_3h.events")
    IntervalSettings.RemoveOverlap("temp/Gata1_inputCtrl-rep1_Gata1_3h_sorted_-_interval.bed", "temp", "input/Gata1_inputCtrl-rep1_Gata1_3h.events")
    end5=time.time()
    print "Time elapsed: "+str(end5-end4)
    print ''' Binning intervals '''
    IntervalSettings.BinInterval("temp/Gata1_inputCtrl-rep1_Gata1_3h_sorted_+_interval_NoOverlapp.bed", "temp", 5) 
    IntervalSettings.BinInterval("temp/Gata1_inputCtrl-rep1_Gata1_3h_sorted_-_interval_NoOverlapp.bed", "temp", 5)
    end6=time.time()
    print "Time elapsed: "+str(end6-end5)
    print ''' Calculating the distribution in a vector form '''
    CalculateDistribution.PerformIntersect("temp/Gata1_inputCtrl-rep1_Gata1_3h_sorted_+_interval_NoOverlapp_binned.bed", 
                                           "temp/Gata1_inputCtrl-rep1_Gata1_3h_sorted_-_interval_NoOverlapp_binned.bed", 
                                           "temp/GATA1_3_ER4_3h_run1_onebed_sorted.bed", "temp", "3hr_run1", 100, 5)
    #CalculateDistribution.PerformIntersect("temp/CTCF_HeLa_Pugh_exo_GEM_events_sorted_+_interval_NoOverlapp_binned.bed", 
    #                                       "temp/CTCF_HeLa_Pugh_exo_GEM_events_sorted_-_interval_NoOverlapp_binned.bed", 
    #                                       "temp/Hela_2_bowtie_unique+colorspace_hg18_sorted_onebed.bed", "temp", "rep-2", 100, 5)
    #CalculateDistribution.PerformIntersect("temp/CTCF_HeLa_Pugh_exo_GEM_events_sorted_+_interval_NoOverlapp_binned.bed", 
    #                                       "temp/CTCF_HeLa_Pugh_exo_GEM_events_sorted_-_interval_NoOverlapp_binned.bed", 
    #                                       "temp/Hela_3_bowtie_unique+colorspace_hg18_sorted_onebed.bed", "temp", "rep-3", 100, 5)
    #CalculateDistribution.CalculateVector("temp/rep-1_+_intersect_100_5.tab", "temp/rep-1_-_intersect_100_5.tab", 
    #                                     "output", "rep-1", 100, 5)
    #CalculateDistribution.CalculateVector("temp/rep-2_+_intersect_100_5.tab", "temp/rep-2_-_intersect_100_5.tab", 
    #                                      "output", "rep-2", 100, 5)
    #CalculateDistribution.CalculateVector("temp/rep-3_+_intersect_100_5.tab", "temp/rep-3_-_intersect_100_5.tab", 
    #                                      "output", "rep-3", 100, 5)
    CalculateDistribution.CalculateMatDis("temp/3hr_run1_+all_intersect_100_5.tab", "temp/3hr_run1_-all_intersect_100_5.tab", 
                                          "output", "3hr_run1", 100, 5)
    end7=time.time()
    print "Time elapsed: "+str(end7-end6)
    
if __name__ == "__main__":
    main()
