'''
Created on Dec 16, 2012

@author: akshaykakumanu
'''
import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as m
cdict = {
  'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
  'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
  'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
}

cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

#data=np.genfromtxt(open("output/clustermerged_top10000.mat"),usecols=(2,3,4,5,6,7,8,9,10,11,12,13,14,
#                                                             15,16,17,18,19,20,21,22,23,24,
#                                                             25,26,27,28,29,30,31,32,33,34,35
#                                                             ,36,37,38,39,40,41,42,43,44,45,46,47
#                                                             ,48,49,50,51,52,53,54,55,56,57,58,59,70,71,72,73,74,75,76,77,78,79,80),delimiter="\t")

data=np.genfromtxt(open("output/Cluster4_random.mat"),delimiter="\t")


fig = plt.figure()
#plt.gray()
plt.pcolormesh(data,cmap=cm)
#plt.pcolormesh(data[:,1:],cmap=cm)
#plt.pcolormesh(data[:,1:],cmap=cm)
plt.colorbar() 
plt.savefig('Cluster_4_random_heatmap.png')
