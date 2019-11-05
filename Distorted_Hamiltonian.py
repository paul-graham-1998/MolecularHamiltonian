# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:03:51 2019

@author: Paul
"""

import numpy as np
import cmath
from Distorted_Indexing import Distorted_momentumLabelToK
from Distorted_Indexing import Distorted_LabelToR
from Distorted_Indexing import Distorted_Labels
from Distorted_Indexing import Distorted_momentumLabels
import matplotlib.pyplot as plt
#n is number of unit cells
def Distorted_Hamiltonian(n):
    Atom_Number = 2 * n
    hamiltonian = np.zeros((Atom_Number,Atom_Number))
    for i in range(Atom_Number-1):
        if (i%2 == 0):
            hamiltonian[i][i+1] = -1
            hamiltonian[i+1][i] = -1
        else:
            hamiltonian[i][i+1] = -0.1
            hamiltonian[i+1][i] = -0.1
    hamiltonian[0][Atom_Number-1] = -0.1
    hamiltonian[Atom_Number-1][0] = -0.1
    plt.matshow(hamiltonian, cmap = 'rainbow')
    return hamiltonian
Distorted_Hamiltonian(6)

def Distorted_Fourier(n):
    Atom_Number = 2 * n
    fourier = np.zeros((Atom_Number,Atom_Number), dtype = complex)
    labels = Distorted_Labels(n)
    Mom_labels = Distorted_momentumLabels(n)
    for x in range(Atom_Number):
        for y in range(Atom_Number):
            myLabel_X = labels[x]
            myLabel_Y = Mom_labels[y]
            k_y = Distorted_momentumLabelToK(myLabel_Y , n)
            r_x = Distorted_LabelToR(myLabel_X)
            fourier[x][y] = (1/np.sqrt(n)) * int(myLabel_X[1] == myLabel_Y[1]) * cmath.exp(1j * k_y * r_x)
    return fourier

def Distorted_Junk(n):
    
    Atom_Number = 2*n
    fourier = Distorted_Fourier(n)
    hamiltonian = Distorted_Hamiltonian(n)
    #plt.matshow(np.abs(fourier))
    #plt.matshow(hamiltonian)  #, cmap = 'rainbow'
    #print(np.dot(np.conjugate(np.transpose(fourier)) , fourier))
    #plt.matshow(np.real(np.dot(np.conjugate(np.transpose(fourier)) , fourier)))
    #print(hamiltonian[20,10])
    #print(np.real(np.dot(np.conjugate(np.transpose(fourier)) , fourier))[20,10])
    transform = np.matmul(np.matmul(np.conjugate(np.transpose(fourier)) , hamiltonian) , fourier)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.set_xticklabels(labels, fontsize = 24)   
    ax1.xaxis.set_tick_params(labelsize=40)
    ax1.yaxis.set_tick_params(labelsize=40)
    #ax1.set_yticklabels(fontsize = 24)
    ax1.set_title('Distorted Fourier Transform Hamiltonian', fontsize = 128)
    matshow1 = ax1.matshow(np.real(transform), cmap = plt.cm.rainbow)
    #transform = np.dot(np.dot(fourier , hamiltonian) , np.conjugate(np.transpose(fourier)))
    #transform = np.abs(transform)
    #plt.matshow(np.abs(transform))
    Energies1 = []
    Energies2 = []
    Distances = []
    for i in range(0,Atom_Number,2):
        #print(np.linalg.eigh(transform[i:i+2 , i:i+2]))
        E,v = np.linalg.eigh(transform[i:i+2 , i:i+2])
        
        Energies1.append(E[0])
        Energies2.append(E[1])
        Distances.append( Distorted_momentumLabelToK((i/2,'a') , n))


    plt.figure()
    plt.plot(Distances, Energies1, '-o' , linewidth = 5, markersize = 20 , label = 'Upper Band')
    plt.plot(Distances, Energies2, '-o' , linewidth = 5, markersize = 20 , label = 'Lower Band')
    plt.hlines((min(Energies1)+max(Energies2))/2 , 0 , np.pi, linewidth = 5 , label = 'Fermi Level')
    plt.title('Distorted Hydrogen Chain Fermi Diagram', fontsize = 100 )
    plt.ylabel('Energy' , fontsize  = 60)
    plt.xlabel('Reciprocal Lattice Location' , fontsize = 60)
    plt.tick_params(labelsize = 40)
    plt.legend(prop={'size': 50} , loc = 'best')

    
#Distorted_Junk(27)
    
   
            