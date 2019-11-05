# -*- coding: utf-8 -
"""
Created on Sun Apr 21 01:43:38 2019

@author: Paul
"""

#Hydrogen chain hamiltonian

import numpy as np
import cmath
from Hydrogen_Indexing import momentumLabelToK
from Hydrogen_Indexing import LabelToR
import matplotlib.pyplot as plt
#n is number of unit cells
def Hamiltonian(n):
    hamiltonian = np.zeros((n,n))
    for i in range(n-1):
        hamiltonian[i][i+1] = 1#-1
        hamiltonian[i+1][i] = 1#-1
    hamiltonian[0][n-1] = 1#-1
    hamiltonian[n-1][0] = 1#-1
    return hamiltonian

def Fourier(n):
    fourier = np.zeros((n,n), dtype = complex)
    for x in range(n):
        for y in range(n):
            myLabel_X = (x , 'a')
            myLabel_Y = y
            k_y = momentumLabelToK(myLabel_Y , n)
            r_x = LabelToR(myLabel_X)
            fourier[x][y] = (1/np.sqrt(n)) * cmath.exp(1j * k_y * r_x)
            #print(fourier[x][y])
    return fourier
def Junk(n):
    fourier = Fourier(n)
    hamiltonian = Hamiltonian(n)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.set_xticklabels(labels, fontsize = 24)   
    ax1.xaxis.set_tick_params(labelsize=40)
    ax1.yaxis.set_tick_params(labelsize=40)
    #ax1.set_yticklabels(fontsize = 24)
    ax1.set_title('Hamiltonian', fontsize = 128)
    matshow1 = ax1.matshow(hamiltonian, cmap = plt.cm.rainbow)
    
    #np.conjugate(np.transpose(M))
    #print(np.dot(np.conjugate(np.transpose(fourier)) , fourier))
    #plt.matshow(np.absolute(np.dot(np.conjugate(np.transpose(fourier)) , fourier)))
    transform = np.matmul(np.matmul(np.conjugate(np.transpose(fourier)) , np.asmatrix(hamiltonian)) , np.asmatrix(fourier))
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.set_xticklabels(labels, fontsize = 24)   
    ax1.xaxis.set_tick_params(labelsize=40)
    ax1.yaxis.set_tick_params(labelsize=40)
    #ax1.set_yticklabels(fontsize = 24)
    ax1.set_title('Hydrogen Fourier Transform Hamiltonian', fontsize = 128)
    matshow1 = ax1.matshow(np.real(transform), cmap = plt.cm.rainbow)
    
    #print(transform)
    #print(transform.shape)
    Energies = []
    Distances = []
    for i in range(n):
        Energies.append(transform[i,i])
        Distances.append( momentumLabelToK(i , n))
    plt.figure()
    plt.plot(Distances , Energies , 'y--' , linewidth = 10 , label = 'Hamiltonian Basis Energies')
    plt.hlines((min(Energies)+max(Energies))/2 , 0 , 2*np.pi, linewidth = 5 , label = 'Fermi Level')
    x = np.linspace(0, 2 * np.pi,100)
    y = 2 * np.cos(x)
    plt.plot(x,y,'b:', linewidth = 5, markersize = 10 , label = '2*cosine()')
    plt.show()
    plt.title('Hydrogen Chain Fermi Diagram', fontsize = 100 )
    plt.ylabel('Energy' , fontsize  = 60)
    plt.xlabel('Reciprocal Lattice Location' , fontsize = 60)
    plt.tick_params(labelsize = 40)
    plt.legend(prop={'size': 50})


Junk(50)
    
   
            