'''
Created on Jan 4, 2013

@author: akshaykakumanu
'''
import matplotlib.pyplot as plt
import numpy as np

data=np.genfromtxt(open("output/Cluster4.mat"),delimiter="\t")

composite = np.arange(data.shape[1]-1).reshape(1,data.shape[1]-1)

for i in xrange(composite.shape[1]):
    composite[0,i]=np.sum(data[:,i+1])
print composite
bins = composite.shape[1]/2.0
poscomp = np.arange(bins).reshape(1,bins)
negcomp = np.arange(bins).reshape(1,bins)

poscomp = composite[0,0:bins]


for i in xrange(int(bins)):
    j=2*bins-1-i
    print composite[0,j],j
    negcomp[0,i]=-1.0*composite[0,j]

xaxis = range(-20,20)
#xaxis = range(-100,100)
plt.plot(xaxis,poscomp[0:40])
plt.plot(xaxis,negcomp[0,0:40])
#plt.plot(xaxis,poscomp[0:200])
#plt.plot(xaxis,negcomp[0,0:200])
plt.savefig('Cluster_4_composite_plot.eps',fmt='eps')