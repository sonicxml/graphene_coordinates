#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import pylab as mpl

def plotgraphene(gc):
    # MATLAB: plot(graphenecoordinates(:,1),graphenecoordinates(:,2),'o')'
    mpl.scatter(gc[:,0],gc[:,1])
    mpl.grid(True)
    mpl.xlabel('Length (angstroms)')
    mpl.ylabel('Width (angstroms)')
    mpl.title('Graphene Lattice')
    mpl.show()
    
def Hamiltonian(gc):
    Norb = gc.shape[0]
    numH = Norb
    x = np.tile(gc[:,0], (1, Norb))
    y = np.tile(gc[:,1], (1, Norb))
    z = np.tile(gc[:,2], (1, Norb))
    R = np.sqrt(np.square(x.T-x) + np.square(y.T-y) + np.square (z.T-z))
    
    H = np.zeros(shape=(numH, numH))
    print numH
    for k in xrange(0, (numH)):
        for l in xrange(0, (numH)):
            if R[k,l] <= 1.45 and R[k,l] >= 1.38:
                H[k,l]=-2.7
#            else:
#                H[k,l]=0
    print H
    np.savetxt('pythonHam.txt', H, delimiter='\t', fmt='%f')
    return H
    
def main():
    # Load coordinate file as array
    gc = np.loadtxt("graphenecoordinates.txt")
    gc.view('i8,i8,i8').sort(order=['f0','f1'], axis=0) # For 64-bit systems - for 32-bit, change 'i8' to 'i4'
    
    # Run functions
    plotgraphene(gc)
    Hamiltonian(gc)
    return 0

if __name__ == '__main__':
	main()

