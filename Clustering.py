'''
Created on Dec 12, 2012

@author: akshaykakumanu
'''
import numpy as np
import random
import math

class cluster:
    Totalcount=0
    def __init__(self,name,leaf,childl,childr,active,nodes):
        self.leaf=leaf
        self.childl=childl
        self.childr=childr
        self.active=active
        self.nodes=nodes
        self.name=name
        cluster.Totalcount += 1

def CreateActiveArray(distrimat):
    activeclusters=[]
    for i in xrange(distrimat.shape[0]):
        activeclusters.extend([cluster("clus"+str(i),1,"","",1,[i])])
    return activeclusters

def NearestNeighbor(distrimat,activeclusters,fromclus):
    maxP=-2
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
    return frompos,neighpos

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
    var=totaldistance/n
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

    
    
        

            
                    
            
        
        
    
        
        