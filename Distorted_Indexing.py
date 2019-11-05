# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:30:38 2019

@author: Paul
"""

import numpy as np
import cmath
#n is number of unit cells
def Distorted_Labels(n):
    LabelsList = []
    for i in range(1,n+1):
        LabelsList.append((i, 'a'))
        LabelsList.append((i, 'b'))
    return LabelsList


def Distorted_LabelToR(myLabel):
    if myLabel[1] == 'a':
        length = 2*(myLabel[0] - 1)
    if myLabel[1] == 'b':
        length = 2*(myLabel[0] - 1) + 0.8     
    return length
    
    
def Distorted_LabelToIndex(myLabel, UnitCell):
    if myLabel[1] == 'a':
        index = 2*(myLabel[0] - 1)
    if myLabel[1] == 'b':
        index = 2*(myLabel[0] - 1) + 1    
    return(index)
    
#Label = (3,'b')
#print(LabelToIndex(Label, 'DH'))
#n ATOMS
def Distorted_momentumLabels(n):
    labels = []
    for i in range(n):
        labels.append((i,'a'))
        labels.append((i, 'b'))
    return labels

def Distorted_momentumLabelToK(myLabel , n):
    length = 2*np.pi*myLabel[0] / (2*n)
    return length