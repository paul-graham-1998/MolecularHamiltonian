# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:57:48 2019

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
#fig = plt.figure()
#ax = fig.add_subplot(111, projection = '3d')
#ax.plot_trisurf(all_kx , all_ky , E)
#n is number of unit cells
def Hamiltonian(n1 , n2 , t, theta , M):
    labels = Hexagonal_Labels(n1 , n2)
    hamiltonian = np.zeros((2 * n1 * n2 , 2 * n1 * n2), dtype = complex)
    dictionary = {}
    for label in labels:
        dictionary[label] = 0
    counter = 0
    for label in labels:
        Index = Hexagonal_LabelToIndex(label , n1 , n2)
        if label[1] == 'a':
            label_New1 = ( ( (label[0][0]+1)%n1 , (label[0][1]-1)) , 'a' )
            label_New2 = ( ( label[0][0] , (label[0][1]+1)) , 'a' )
            label_New3 = ( ( (label[0][0]+1)%n1 , label[0][1] ) , 'a' )
            
            label_New4 = ( label[0]    , 'b' )
            label_New5 = ( ( label[0][0] , (label[0][1]-1)) , 'b' )
            label_New6 = ( ( (label[0][0]+1)%n1 , (label[0][1]-1)) , 'b' )
            multiplier = 1#SHOULD BE 1
            hamiltonian[Index , Index] = M
        else:
            label_New1 = ( ( (label[0][0]+1)%n1 , (label[0][1]-1)) , 'b' )
            label_New2 = ( ( label[0][0] , (label[0][1]+1)) , 'b' )
            label_New3 = ( ( (label[0][0]+1)%n1 , label[0][1]) , 'b' )
            
            label_New4 = ( label[0]    , 'a' )
            label_New5 = ( ( label[0][0] , (label[0][1]+1)) , 'a' )
            label_New6 = ( ( (label[0][0]-1)%n1 , (label[0][1]+1)) , 'a' )
            multiplier = -1#SHOULD BE -1
            hamiltonian[Index , Index] = -1*M   

#        dictionary[label] += 3
#        dictionary[label_New1] += 1
#        dictionary[label_New2] += 1
#        dictionary[label_New3] += 1

        Index_New1 = Hexagonal_LabelToIndex(label_New1 , n1 , n2)
        Index_New2 = Hexagonal_LabelToIndex(label_New2 , n1 , n2)
        Index_New3 = Hexagonal_LabelToIndex(label_New3 , n1 , n2)
        
        Index_New4 = Hexagonal_LabelToIndex(label_New4 , n1 , n2)
        Index_New5 = Hexagonal_LabelToIndex(label_New5 , n1 , n2)
        Index_New6 = Hexagonal_LabelToIndex(label_New6 , n1 , n2)
#        print(label)
#        print([Index_New1 , Index_New2 , Index_New3, Index_New4 , Index_New5 , Index_New6])
        if (label_New1[0][1]  < n1 and label_New1[0][1] > -1):
            hamiltonian[Index , Index_New1] -= 1 * t * cmath.exp(multiplier * 1j * theta)
            hamiltonian[Index_New1 , Index] -= 1 * t * cmath.exp(multiplier * -1 * 1j * theta)

        if (label_New2[0][1]  < n1 and label_New2[0][1] > -1):
            hamiltonian[Index , Index_New2] -= 1 * t * cmath.exp(multiplier * 1j * theta)
            hamiltonian[Index_New2 , Index] -= 1 * t * cmath.exp(multiplier * -1 * 1j * theta)

        if (label_New3[0][1]  < n1 and label_New3[0][1] > -1):
            hamiltonian[Index , Index_New3] -= 1 * t * cmath.exp(multiplier * -1 * 1j * theta)
            hamiltonian[Index_New3 , Index] -= 1 * t * cmath.exp(multiplier * 1j * theta)

        if (label_New4[0][1]  < n1 and label_New4[0][1] > -1):
            hamiltonian[Index , Index_New4] = -1
            hamiltonian[Index_New4 , Index] = -1
        if (label_New5[0][1]  < n1 and label_New5[0][1] > -1):
            hamiltonian[Index , Index_New5] = -1
            hamiltonian[Index_New5 , Index] = -1
        if (label_New6[0][1]  < n1 and label_New6[0][1] > -1):
            hamiltonian[Index , Index_New6] = -1
            hamiltonian[Index_New6 , Index] = -1
#        hamiltonian[Index_New4 , Index] = -1
#        hamiltonian[Index_New5 , Index] = -1
#        hamiltonian[Index_New6 , Index] = -1
        
#        hamiltonian2[Index , Index_New1] = 1
#        hamiltonian2[Index , Index_New2] = 1
#        hamiltonian2[Index , Index_New3] = 1
#        hamiltonian3[Index , Index_New1] = -1
#        hamiltonian3[Index , Index_New2] = -1
#        hamiltonian3[Index , Index_New3] = -1
        counter += 3
        #plt.matshow(np.real(hamiltonian))
        
#    #plt.matshow(np.real(hamiltonian2), cmap = 'rainbow')
#    #plt.matshow(np.real(hamiltonian3), cmap = 'rainbow')
#    fig = plt.figure()
#    ax1 = fig.add_subplot(111)
#    #ax1.set_xticklabels(labels, fontsize = 24)   
#    ax1.xaxis.set_tick_params(labelsize=40)
#    ax1.yaxis.set_tick_params(labelsize=40)
#    #ax1.set_yticklabels(fontsize = 24)
#    ax1.set_title('Hamiltonian', fontsize = 128)
#    matshow1 = ax1.matshow(np.real(hamiltonian))
#    fig.colorbar(matshow1)

    
    return hamiltonian

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
            #print([k_y , r_x])
            fourier[i][j] = (1/np.sqrt(n1*n2)) * int(i_label[1] == j_label[1]) * cmath.exp(1j * (k_y[0] * r_x[0] + k_y[1] * r_x[1]))
            #print(fourier[x][y])
    return fourier

def Junk(n1 , n2, t , theta , M,  G1 , G2):

    fourier = Fourier(n1, n2)
    hamiltonian = Hamiltonian(n1, n2 , t, theta, M)
#    for i in range(2*n1*n2):
#        for j in range(2*n1*n2):
#            print([hamiltonian[i , j] , '        ',np.conjugate(np.transpose(hamiltonian))[i,j]])
    #transform = np.matmul(np.matmul(fourier,hamiltonian) , np.conjugate(np.transpose(fourier)) )
    transform = np.matmul(np.matmul(np.conjugate(np.transpose(fourier)),hamiltonian) , fourier)

    #plt.np.real(transform)
    
    if G1:
        plt.matshow(np.absolute(fourier))
        plt.matshow(np.real(hamiltonian), cmap = 'rainbow')
        plt.matshow(np.imag(hamiltonian), cmap = 'rainbow')
        plt.matshow(np.absolute(np.matmul(np.conjugate(np.transpose(fourier)) , fourier)))
        plt.matshow(np.absolute(np.subtract(np.conjugate(np.transpose(hamiltonian)) , hamiltonian)))        
        plt.matshow(abs(transform))
        plt.matshow(np.real(transform))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.set_xticklabels(labels, fontsize = 24)   
    ax1.xaxis.set_tick_params(labelsize=40)
    ax1.yaxis.set_tick_params(labelsize=40)
    #ax1.set_yticklabels(fontsize = 24)
    ax1.set_title('Edge Fourier Transform Hamiltonian', fontsize = 128)
    matshow1 = ax1.matshow(np.real(transform), cmap = plt.cm.rainbow)


    
    Energies1 = []
    Energies2 = []
    EnergiesDiff = []
    Vectors1 = []
    Vectors2 = []
    k_x = []
    k_y = []
    
#    for i in range(n2):
#        label = ((i , 0) , 'a')
#        C = Hexagonal_LabelToIndex(label , n1,n2)
#        E,v = np.linalg.eigh(transform[C:C+2*n2 , C:C+2*n2]) 
#        for j in E:
#            Energies1.append(j)
#            k_x.append(Hexagonal_momentumLabelToK(label , n1 , n2)[0])

        #k_x.append(Hexagonal_momentumLabelToK(label , n1 , n2)[0])
    for j in range(n2):
        label = ((0 , j) , 'a')
        C = Hexagonal_LabelToIndex(label , n1,n2)
        E,v = np.linalg.eigh(transform[C:C+2*n2 , C:C+2*n2]) 
        for i in E:
            Energies1.append(i)
            k_y.append(Hexagonal_momentumLabelToK(label , n1 , n2)[1])

        #k_x.append(Hexagonal_momentumLabelToK(label , n1 , n2)[0])
    if G2:
        plt.figure()
        
        plt.plot(k_y , Energies1 , "o")
        plt.title('Haldane open x boundary, M=2.0', fontsize = 100 )
        plt.ylabel('Energy' , fontsize  = 60)
        plt.xlabel('ky - Reciprocal Lattice Location' , fontsize = 60)
        plt.tick_params(labelsize = 40)
    #print(Vectors1)
    #return (Vectors1)
    #return vector
#Junk(6 ,6 , 0.3, 0.7 , 2, False, True)

def Junk2(n1, n2 , t , theta):
    Energies = []
    Ms = []
    for M in np.linspace(0, 2, 50):
        print(M)
        Energies.append(Junk(n1 ,n2, t, theta, M, False, False))
        Ms.append(M)
    plt.figure()
    plt.scatter(Ms, Energies)
    plt.show()
#Junk2(9, 9, 0.3, 0.7)
def phaseGetter(vector1 , vector2):  
    Phase = np.angle(np.dot(np.transpose(vector1) , vector2))
#    if Phase > 2*np.pi:
#        Phase -= 2*np.pi
#    if Phase < 2*np.pi:
#        Phase += 2*np.pi
    return Phase

def Berry(n1 , n2 , t , theta, M , G1):
    Ms = []
    kxs = []
    kys = []
    Phases = []
    vectorList = Junk(n1 , n2 , t , theta , M , False , False)
    print(len(vectorList))
    for kx in range(n1):
        for ky in range(n2):
            #Total_Phase = 0
            i1 = kx * n1 + ky
            i2 = kx * n1 + (ky + 1)
            i4 = (kx + 1)%n1 * n1 + ky
            i3 = (kx + 1)%n1 * n1 + (ky + 1)
#            print([i1 , i2 , i3 , i4])
#            print( vectorList[i1] )
#            print( vectorList[i2] )
#            print( vectorList[i3] )
#            print( vectorList[i4] )
#            print( "            " )
            TotalPhase = 1
            if i2 < n2 & i4 < n2:
                TotalPhase *= np.dot(np.transpose(np.conjugate(vectorList[i1])) , vectorList[i2])
                TotalPhase *= np.dot(np.transpose(np.conjugate(vectorList[i2])) , vectorList[i3])
                TotalPhase *= np.dot(np.transpose(np.conjugate(vectorList[i3])) , vectorList[i4])
                TotalPhase *= np.dot(np.transpose(np.conjugate(vectorList[i4])) , vectorList[i1])
                TotalPhase = np.angle(TotalPhase)
                while TotalPhase > np.pi:
                    TotalPhase -= 2*np.pi
                while TotalPhase < -1 * np.pi:
                    TotalPhase += 2*np.pi
                #print([kx , ky , "       " , TotalPhase])
                Phases.append(TotalPhase)
                kxs.append(kx)
                kys.append(ky)
    if G1:
        fig1 = plt.figure()
        ax = fig1.add_subplot(111, projection = '3d')
        ax.plot_trisurf(kxs , kys, Phases)
    #print(Vectors1)
    return(sum(Phases) / (2*np.pi))

#Berry(9 , 9 , 0.3 , 0.7 , 1.1, True)

def Berry2(n1 , n2, t , theta):
    Ms = []
    Chern = []
    for M in np.linspace(0 , 2 , 20):
        print(M)
        Ms.append(M)
        Chern.append(Berry(n1 , n2 , t , theta , M , False))
    plt.figure()
    plt.plot(Ms , Chern)
        
#Berry2(9 , 9 , 0.3 , 0.7)
