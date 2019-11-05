# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 18:29:23 2019

@author: Paul
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:30:31 2019

@author: Paul
"""

import numpy as np
import cmath
from Hexagonal_Indexing import Hexagonal_momentumLabelToK
from Hexagonal_Indexing import Hexagonal_LabelToR
from Hexagonal_Indexing import Hexagonal_Labels
from Hexagonal_Indexing import Hexagonal_LabelToIndex

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#n is number of unit cells
def Hamiltonian(n1 , n2 , M):
    labels = Hexagonal_Labels(n1 , n2)
    hamiltonian = np.zeros((2 * n1 * n2 , 2 * n1 * n2))
    for label in labels:
        Index = Hexagonal_LabelToIndex(label , n1 , n2)
        if label[1] == 'a':
            label_New1 = ( label[0]    , 'b' )
            label_New2 = ( ( label[0][0] , (label[0][1]-1)%n2 ) , 'b' )
            label_New3 = ( ( (label[0][0]+1)%n1 , (label[0][1]-1)%n2 ) , 'b' )
            hamiltonian[Index , Index] = M

        else:
           label_New1 = ( label[0]    , 'a' )
           label_New2 = ( ( label[0][0] , (label[0][1]+1)%n2 ) , 'a' )
           label_New3 = ( ( (label[0][0]-1)%n1 , (label[0][1]+1)%n2 ) , 'a' )
           hamiltonian[Index , Index] = -1*M

        Index_New1 = Hexagonal_LabelToIndex(label_New1 , n1 , n2)
        Index_New2 = Hexagonal_LabelToIndex(label_New2 , n1 , n2)
        Index_New3 = Hexagonal_LabelToIndex(label_New3 , n1 , n2)
        hamiltonian[Index , Index_New1] = -1
        hamiltonian[Index , Index_New2] = -1
        hamiltonian[Index , Index_New3] = -1
    return hamiltonian
#print(Hamiltonian(5,5))

def Fourier(n1 , n2):
    Atom_Number = 2*n1*n2
    fourier = np.zeros((Atom_Number,Atom_Number), dtype = complex)
    labels = Hexagonal_Labels(n1, n2)
    for i in range(Atom_Number):
        for j in range(Atom_Number):
            i_label = labels[i]
            j_label = labels[j]
            k_y = Hexagonal_momentumLabelToK(j_label , n1 , n2)
            r_x = Hexagonal_LabelToR(i_label)
            fourier[i][j] = (1/np.sqrt(n1*n2)) * int(i_label[1] == j_label[1]) * cmath.exp(1j * (k_y[0] * r_x[0] + k_y[1] * r_x[1]))
            #print(fourier[x][y])
    return fourier
def Junk(n1 ,n2 , M):
    fourier = Fourier(n1, n2)
    hamiltonian = Hamiltonian(n1, n2 , M)
    #plt.matshow(hamiltonian)
    #plt.matshow(np.absolute(np.matmul(np.conjugate(np.transpose(fourier)) , fourier)))
    transform = np.dot(np.dot(np.conjugate(np.transpose(fourier)) , np.asmatrix(hamiltonian)) , np.asmatrix(fourier))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.set_xticklabels(labels, fontsize = 24)   
    ax1.xaxis.set_tick_params(labelsize=40)
    ax1.yaxis.set_tick_params(labelsize=40)
    #ax1.set_yticklabels(fontsize = 24)
    ax1.set_title('Nitride Fourier Transform Hamiltonian', fontsize = 128)
    matshow1 = ax1.matshow(np.real(transform), cmap = plt.cm.rainbow)


    #transform = np.real(transform)
    #plt.matshow(abs(transform))
    Energies1 = []
    Energies2 = []
    k_x = []
    k_y = []
#    for i in range(0,Atom_Number,2):
#        #print(np.linalg.eigh(transform[i:i+2 , i:i+2]))
#        E,v = np.linalg.eigh(transform[i:i+2 , i:i+2])
#        
#        Energies1.append(E[0])
#        Energies2.append(E[1])
#        Distances.append( Distorted_momentumLabelToK((i/2,'a')[0] , n))
#        Distances.append( Distorted_momentumLabelToK((i/2,'a')[1] , n))
    colorList = ['black' , 'red' , 'gold' , 'darkgreen' , 'aqua' , 'navy' , 'greenyellow' , 'hotpink' , 'saddlebrown']
    fig2 = plt.figure()
    for i in range(n1):
        tempLoc = 0
        tempLocs = []
        tempEnergies1 = []
        tempEnergies2 = []
        for j in range(n2):
            label = ((i , j) , 'a')
            C = Hexagonal_LabelToIndex(label , n1,n2)
            E,v = np.linalg.eigh(transform[C:C+2 , C:C+2])
            Energies1.append(E[0])
            Energies2.append(E[1])       
            kx = Hexagonal_momentumLabelToK(label , n1 , n2)[0]
            ky = Hexagonal_momentumLabelToK(label , n1 , n2)[1]
            k_x.append(kx)
            k_y.append(ky)
            #print([kx , ky])
            tempLocs.append(ky)
            tempLoc = kx
            tempEnergies1.append(E[0])
            tempEnergies2.append(E[1])
        plt.plot(tempLocs , tempEnergies1 , c = colorList[i%9], linewidth = 5, label = 'kx of ' + str(tempLoc))
        plt.plot(tempLocs , tempEnergies2 , c = colorList[i%9], linewidth = 5)
    plt.legend(prop={'size': 20} , loc = 'best')
    plt.title('Slices of  Graphene Lattice Fermi Diagram', fontsize = 100 )
    plt.ylabel('Energy' , fontsize  = 60)
    plt.xlabel('y Momentum Space' , fontsize = 60)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.plot_trisurf(k_x , k_y , Energies2)
    ax.plot_trisurf(k_x , k_y , Energies1)
    ax.set_title('Boron Nitride Lattice Fermi Diagram, M=1', fontsize = 100)
    ax.set_xlabel('x Momentum Space', fontsize = 50, labelpad = 60)
    ax.set_ylabel('y Momentum Space', fontsize = 50, labelpad = 60)
    ax.set_zlabel('Energy', fontsize = 50, labelpad = 30)
    ax.tick_params(labelsize = 40)
    
    
#Junk(3 ,3  , 1)
    