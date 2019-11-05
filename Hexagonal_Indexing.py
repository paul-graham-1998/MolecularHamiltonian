# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:41:06 2019

@author: Paul
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
#GRAPHENE LATTICE STUFF



def Hexagonal_Labels(n1 , n2):
    labels = []
    for x in range(n1):
        for y in range(n2):
            labels.append(((x,y),'a'))
            labels.append(((x,y),'b'))
    return labels
#def Labels(n1 , n2):
#    labels = []
#    for y in range(n2):
#        for x in range(n1):
#            labels.append(((x,y),'a'))
#            labels.append(((x,y),'b'))
#    return labels

def Hexagonal_LabelToR(myLabel):
#    x = (myLabel[0][0]-1) + 0.5*(myLabel[0][1]-1)
#    y = (np.sqrt(3)/2) * (myLabel[0][1]-1)
    x = myLabel[0][0] + 0.5*myLabel[0][1]
    y = (np.sqrt(3)/2) * myLabel[0][1]
    if myLabel[1] == 'b':
       y += 1/np.sqrt(3)
    return((x,y))
    
def Hexagonal_LabelToIndex(myLabel , n1 , n2):
    index = 0
    index += 2 * myLabel[0][1]
    index += 2 * n2 * myLabel[0][0]
    if myLabel[1] == 'b':
        index += 1
    return index

#labels =  Hexagonal_Labels(3,3)
##print(labels)
#for i in range(len(labels)):
##    print([i , Hexagonal_LabelToIndex(labels[i] , 3,3)])
#    print(labels[i])
#

#Momentum space graphene lattice    
def Hexagonal_momentumLabels(n1 , n2):
    momLabels = []
    for m in range(n1):
        for n in range(n2):
            momLabels.append((m, n) , 'a')
            momLabels.append((m, n) , 'b')      
    return momLabels

def Hexagonal_momentumLabelToK(myLabel, n1, n2):
    x = (2*np.pi/n1) * (myLabel[0][0])
    y = (-2*np.pi/(n1*np.sqrt(3))) * (myLabel[0][0]) + (4*np.pi/(n2*np.sqrt(3))) * (myLabel[0][1])

    return((x,y))
#def BrillouinZone():



def Plotter(n1, n2):
    xCoords = []
    yCoords = []
    markers = []
    indexes = []
    labels = Hexagonal_Labels(n1, n2)
    
    for label in labels:
        xCoords.append(Hexagonal_LabelToR(label)[0])
        yCoords.append(Hexagonal_LabelToR(label)[1])
        markers.append(str(label))
        indexes.append(str(Hexagonal_LabelToIndex(label,n1,n2)))
    plt.figure()
    plt.scatter(xCoords, yCoords, s =500 , c = 'black')
    for i in range(len(labels)):
        plt.annotate(markers[i] , (xCoords[i], yCoords[i]), fontsize = 50)
#NEXT JUNK------------------------------
    for label in labels:
        if label[1] == 'a':
            label_New1 = ( ( (label[0][0]+1)%n1 , (label[0][1]-1)%n2 ) , 'a' )
            label_New2 = ( ( label[0][0] , (label[0][1]+1)%n2 ) , 'a' )
            label_New3 = ( ( (label[0][0]+1)%n1 , label[0][1] ) , 'a' )
        else:
            label_New1 = ( ( (label[0][0]+1)%n1 , (label[0][1]-1)%n2 ) , 'b' )
            label_New2 = ( ( label[0][0] , (label[0][1]+1)%n2 ) , 'b' )
            label_New3 = ( ( (label[0][0]+1)%n1 , label[0][1]) , 'b' )
        print([label , label_New1 , label_New2 , label_New3])
        initial = Hexagonal_LabelToR(label)
        first = Hexagonal_LabelToR(label_New1)
        second = Hexagonal_LabelToR(label_New2)
        third = Hexagonal_LabelToR(label_New3)
#        if label[1] == 'a':
#            plt.arrow(initial[0] , initial[1] , first[0] - initial[0] , first[1] - initial[1] , length_includes_head = True, width = 0.01 , head_width = 0.1, facecolor = 'r')
#            plt.arrow(initial[0] , initial[1] , second[0] - initial[0] , second[1] - initial[1] , length_includes_head = True, width = 0.01 , head_width = 0.1, facecolor = 'b')
#            #plt.arrow(initial[0] , initial[1] , third[0] - initial[0] , third[1] - initial[1], length_includes_head = True , width = 0.01 , head_width = 0.1, facecolor = 'g')
#            plt.arrow(third[0] , third[1] , initial[0] - third[0], initial[1] - third[1], length_includes_head = True , width = 0.01 , head_width = 0.1, facecolor = 'g')        
#        if label[1] == 'b':
#            plt.arrow(first[0], first[1] , initial[0] - first[0], initial[1] - first[1], length_includes_head = True, width = 0.01 , head_width = 0.1, facecolor = 'r')
#            plt.arrow(second[0] , second[1] , initial[0] - second[0], initial[1] - second[1], length_includes_head = True, width = 0.01 , head_width = 0.1, facecolor = 'b')
#            #plt.arrow(third[0] , third[1] , initial[0] - third[0], initial[1] - third[1], length_includes_head = True , width = 0.01 , head_width = 0.1, facecolor = 'g')
#            plt.arrow(initial[0] , initial[1] , third[0] - initial[0] , third[1] - initial[1], length_includes_head = True , width = 0.01 , head_width = 0.1, facecolor = 'g')
#
    
    
def MomentumPlotter(n1, n2):
    xCoords = []
    yCoords = []
    markers = []
    labels = Hexagonal_momentumLabels(n1, n2)
    print(labels)
    for label in labels:
        xCoords.append(Hexagonal_momentumLabelToK(label,n1,n2)[0])
        yCoords.append(Hexagonal_momentumLabelToK(label,n1,n2)[1])
        markers.append(str(label))
    plt.figure()
    plt.scatter(xCoords, yCoords, s = 350)
    for i in range(len(labels)):
        plt.annotate(labels[i] , (xCoords[i], yCoords[i]), fontsize = 10)
    plt.title('Momentum Space')
    plt.xlabel('kx', fontsize = 20)
    plt.ylabel('ky', fontsize = 20)
#Plotter(4,4)
#MomentumPlotter(9,9)

