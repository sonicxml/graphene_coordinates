#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from pylab import *

def plotgraphene(gc):
    # MATLAB: plot(graphenecoordinates(:,1),graphenecoordinates(:,2),'o')'
    scatter(gc[:,0],gc[:,1])
    grid(True)
    xlabel('Length (angstroms)')
    ylabel('Width (angstroms)')
    title('Graphene Lattice')
    show()
    
def Hamiltonian(gc):
    Norb = gc.shape[0]
    numH = Norb
    x = tile(gc[:,0], (1, Norb))
    y = tile(gc[:,1], (1, Norb))
    z = tile(gc[:,2], (1, Norb))
    R = sqrt(square(x.T-x) + square(y.T-y) + square (z.T-z))
    
    H = zeros(shape=(numH, numH))
    for k in xrange(0, (numH)):
        for l in xrange(0, (numH)):
            if R[k,l] <= 1.45 and R[k,l] >= 1.38:
                H[k,l]=-2.7
#            else:
#                H[k,l]=0
    savetxt('pythonHam.txt', H, delimiter='\t', fmt='%f')
    return H
    
def main():
    # Load coordinate file as array
    gc = loadtxt("graphenecoordinates.txt")
    gc.view('i8,i8,i8').sort(order=['f0','f1'], axis=0) # For 64-bit systems - for 32-bit, change 'i8' to 'i4'
    
    # Run functions
    plotgraphene(gc)
    Hamiltonian(gc)
    return 0

if __name__ == '__main__':
	main()

