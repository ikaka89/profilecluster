'''
Created on Dec 21, 2012

@author: akshaykakumanu
'''

import numpy as np
import random
import math

class cluster:
    Totalcount=0
    def __init__(self,name,leaf,childl,childr,active,nodes,pos):
        self.leaf=leaf
        self.childl=childl
        self.childr=childr
        self.active=active
        self.nodes=nodes
        self.name=name
        self.pos=pos
        cluster.Totalcount += 1

def CreateActiveArray(distrimat,scansize,intervalsize,binsize):
    activeclusters=[]
    scanstart=(distrimat.shape[1]/2-(2*scansize/binsize))/2
    for i in xrange(distrimat.shape[0]):
        activeclusters.extend([cluster("clus"+str(i),1,"","",1,[i],[scanstart])])
    return activeclusters

def NearestNeighbor(distrimat,activeclusters,fromclus):
    maxP=-2
    #neigh = cluster("dummy",0,"","",0,[])
    averagex=np.zeros((1,distrimat.shape[1]))
    frompos=0
    neighpos=0
    for xl in list(fromclus.nodes):
        averagex=distrimat[xl,:]+averagex
    averagex=averagex/(len(fromclus.nodes)*1.0)
    xadjusted = (averagex-np.mean(averagex[0,:]))/(np.std(averagex[0,:])*1.0)    
    for cpos in xrange(len(activeclusters)):
        if activeclusters[cpos] != fromclus:
            averagey = np.zeros((1,distrimat.shape[1]))
            for l in list(activeclusters[cpos].nodes):
                averagey=np.transpose(distrimat[l,:])+averagey
            averagey=averagey/(len(activeclusters[cpos].nodes)*1.0)
            averageyt=np.transpose(averagey)
            yadjusted=(averageyt-np.mean(averageyt[:,0]))/(np.std(averageyt[:,0])*1.0)                
            product = np.dot(xadjusted,yadjusted)
            currentP=product[0,0]/(averagey.shape[1]-1.0)
            if currentP > maxP:
                maxP=currentP
                neighpos=cpos
        else:
            frompos=cpos
    #print frompos, neighpos,maxP
    return frompos,neighpos

def NearestNeighborScan(distrimat,activeclusters,fromclus,scansize,intervalsize,binsize):
    maxP=-2
    maxX=0
    maxY=0
    frompos=0
    neighpos=0
    averagex=np.zeros((1,2*2*scansize/binsize))
    if fromclus.leaf == 1:
        averagex[0,0:2*scansize/binsize]=[distrimat[fromclus.nodes[0],fromclus.pos[0]:fromclus.pos[0]+(scansize*2/binsize)+1]]
        averagex[0,2*scansize/binsize:]=[distrimat[fromclus.nodes[0],-1*fromclus.pos[0]-scansize*2/binsize:-1*fromclus.pos[0]]]
    else:
        leftchildx = fromclus.childl
        rightchildx=fromclus.childr
        leftx=fromclus.leftvec[0,fromclus.leftpos:fromclus.leftpos+intervalsize/binsize]*len(leftchildx.nodes)
        rightx=fromclus.righvec[0,fromclus.rightpos:fromclus.rightpos+intervalsize/binsize]*len(rightchildx.nodes)
        averagex=(leftx+rightx)/(1.0*(len(leftchildx.nodes)+len(rightchildx.nodes)))
    for cpos in xrange(len(activeclusters)):
        if activeclusters[cpos] != fromclus:
            averagey = np.zeros((1,distrimat.shape[1]))
            if activeclusters[cpos].leaf == 1:
                averagey[0,0:2*scansize/binsize]=[distrimat[activeclusters[cpos].nodes[0],activeclusters[cpos].pos[0]:activeclusters[cpos].pos[0]+scansize*2/binsize+1]]
                averagey[0,2*scansize/binsize:]=[distrimat[activeclusters[cpos].nodes[0],-1*activeclusters[cpos].pos[0]-scansize*2/binsize:-1*activeclusters[cpos].pos[0]]]
            else:
                leftchildy=activeclusters[cpos].childl
                rightchildy=activeclusters[cpos].cjildr
                lefty=activeclusters[cpos].leftvec[0,activeclusters[cpos].leftpos:activeclusters[cpos].
                                                   lefpos+intervalsize/binsize]*len(leftchildy.nodes)
                righty=activeclusters[cpos].rightvec[0,activeclusters[cpos].rightpos:activeclusters[cpos].
                                                     rightpos+intervalsize/binsize]*len(rightchildy.nodes)
                averagey=(lefty+righty)/(1.0*(len(leftchildy.nodes)+len(rightchildy.nodes)))
                averageyt=np.transpose(averagey)
                endscanx = intervalsize/(binsize*1.0)
                startscanx = 0
                maxinterval=-2
                maxxstart=0
                maxystart=0
                while endscanx < len(scansize):
                    pcx = averagex[0,startscanx:endscanx]
                    pcxadjusted=pcx-np.mean(pcx[0,:])/(1.0*np.std(pcx[0,:]))
                    endscany=intervalsize/(binsize*1.0)
                    startscany=0
                    while endscany < len(scansize):
                        pcy = averageyt[startscany:endscany,0]
                        pcyadjusted = pcy-np.mean(pcy[:,0])/(1.0*np.std(pcy[:,0]))
                        product = np.dot(pcxadjusted,pcyadjusted)
                        currentinterP = product[0,0]/(pcyadjusted.shape[1]*1.0)
                        if currentinterP > maxinterval:
                            maxinterval = currentinterP
                            maxxstart = startscanx
                            maxystart = startscany
                            maxxvec = averagex[0,startscanx:]
                            maxyvec = averageyt[startscany:,0]
                            if scansize/binsize - maxxvec.shape[0] > scan
            if maxinterval > maxP:
                    maxP =maxinterval
                    maxX = maxxstart
                    maxY = maxystart
                    neighpos=cpos
        else:
            frompos=cpos
        
    return frompos,neighpos,maxX,maxY

def InitializeStack(activeclusters):
    stack = []
    index = random.randint(0,len(activeclusters)-1)
    stack.extend([activeclusters[index]])
    return stack

def MergeCluster(frompos,neighpos,activeclusters,stack):
    if len(stack) >=2:
        if activeclusters[neighpos] != stack[-2]:
            stack.extend([activeclusters[neighpos]])
        else:
            a=stack.pop()
            b=stack.pop()
            if a!= activeclusters[frompos] and b != activeclusters[neighpos]:
                print "TROUBLE"
            merged = cluster("clusmerged",0,activeclusters[frompos],activeclusters[neighpos],1,[])
            merged.nodes.extend(activeclusters[frompos].nodes)
            merged.nodes.extend(activeclusters[neighpos].nodes)
            if frompos > neighpos:
                del(activeclusters[frompos])
                del(activeclusters[neighpos])
            else:
                del(activeclusters[neighpos])
                del(activeclusters[frompos])
            stack.extend([merged])
            activeclusters.extend([merged])
    else:
        stack.extend([activeclusters[neighpos]])
    return stack,activeclusters

def MergeClusterScan(frompos,neighpos,maxX,maxY,activeclusters,stack):
    if len(stack) >= 2:
        if activeclusters[neighpos] != stack[-2]:
            stack.extend([activeclusters[neighpos]])
        else:
            a=stack.pop()
            b=stack.pop()
            if a != activeclusters[frompos] and b!=activeclusters[neighpos]:
                print "Trouble"
            merged = cluster("clusmerged",0,activeclusters[frompos],activeclusters[neighpos],1,[],[])
            merged.nodes.extend(activeclusters[frompos].nodes)
            merged.nodes.extend(activeclusters[neighpos].nodes)
            merged.leftpos = []
            merged.rightpos = []
            
            
    

def CheckActive(activeclusters,fromclus):
    value =0
    for c in activeclusters:
        if c.active == 1 and c!= fromclus:
            value = 1
            break
    return value
 
def BIC(activeclusters,distrimat):
    k=len(activeclusters)
    n=distrimat.shape[0]
    totaldistance=0.000000001
    for cpos in xrange(len(activeclusters)):
        averagevec=np.zeros((1,distrimat.shape[1]))
        for i in activeclusters[cpos].nodes:
            averagevec = distrimat[i]+averagevec
        averagevec = averagevec/(len(activeclusters[cpos].nodes)*1.0)
        for i in activeclusters[cpos].nodes:
            distancei = distrimat[i,:]-averagevec[0,:]
            distancei=distancei**2
            totaldistance=totaldistance+np.sum(distancei)
            #print np.sum(distancei)
    var=totaldistance/n
    #print var
    bic=n*math.log(var)+k*math.log(n)
    return bic

def Postorder(root):
    order=[]
    parentstack=[]
    if root.leaf!=1:
        current=root
        while len(order) != len(root.nodes):
            if current.leaf!=1:
                parentstack.extend([current])
                current=current.childl
            else:
                order.extend([current.nodes[0]])
                if len(parentstack) != 0:
                    temp=parentstack.pop()
                    current=temp.childr
    else:
        order.extend([root.nodes[0]])
    return order