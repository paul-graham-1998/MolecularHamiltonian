# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 22:11:53 2019

@author: Paul
"""
import numpy as np
import cmath
#N is number of unit cells
def Labels(N):
    LabelsList = []
    for i in range(1,N+1):
        LabelsList.append((i, 'a'))
    return LabelsList 

def LabelToR(myLabel):
    if myLabel[1] == 'a':
        length = (myLabel[0] - 1)
    return length

    
    
def LabelToIndex(myLabel):
    if myLabel[1] == 'a':
        index = (myLabel[0] - 1) 
    return(index)
    
def momentumLabels(n):
    labels = []
    for i in range(1,n):
        labels.append(i)
    
def momentumLabelToK(myLabel , n):
#    if myLabel[1] == 'a':
#        length = 2*np.pi*myLabel[0] / (2*n)
#    if myLabel[1] == 'b':
#        length = 2*np.pi*myLabel[0] /(2*n + 0.8)
    length = 2*np.pi*myLabel / (n)
    return length

#labels = Labels(5)
#for label in labels:
#    print(LabelToR(label))
#    print(LabelToIndex(label))