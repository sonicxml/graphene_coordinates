#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from pylab import *
import scipy as sp
import scipy.linalg
#from IPython.core.debugger import Tracer; debug_here = Tracer()

def grmagnus(alpha, beta, betad):
    # From J. Phys. F Vol 14, 1984, 1205, M P Lopez Sancho, J. Rubio
    # 20-50 % faster than gravik
    tmp = linalg.inv(alpha)   # Inverse part of Eq. 8
    t=dot(-tmp, betad)              # Eq. 8 (t0)
    tt=dot(-tmp, beta)              # Eq. 8 (t0 tilde)
    T=t.copy()                       # First term in Eq. 16
    Id=eye(alpha.shape[0])       # Save the identity matrix
    Toldt=Id.copy()                  # Product of tilde t in subsequent terms in Eq. 16
    change=1                  # Convergence measure
    counter=0                 # Just to make sure no infinite loop

    eta = 0.000001
    etan = 0.0000000000001
    while linalg.norm(change) > eta and counter < 100:
        counter += 1
        Toldt = dot(Toldt, tt)
        print (1/(linalg.cond(Id - dot(t,tt) - dot(tt,t))))
        if (1/(linalg.cond(Id - dot(t,tt) - dot(tt,t)))) < etan:
            g = 0
            nan_inf_flag = 1
            return g, nan_inf_flag
        tmp = linalg.inv(Id - dot(t,tt) - dot(tt,t))
        
        t = dot(dot(tmp, t), t)
        tt = dot(dot(tmp, tt), tt)
        change = dot(Toldt, t)
        T = T + change
        
        if isnan(change).sum() or isinf(change).sum():
            g = 0
            nan_inf_flag = 1
            return g, nan_inf_flag
    
    print (1/(linalg.cond(alpha + dot(beta, T))))
    if (1/(linalg.cond(alpha + dot(beta, T)))) < etan:
        g = 0
        nan_inf_flag = 1
        return g, nan_inf_flag
    
    g = linalg.inv(alpha + dot(beta, T))
    
    gn = abs(g - linalg.inv(alpha - dot(dot(beta, g), betad)))
    if gn.max() > 0.001 or counter > 99:
        g = 0
        nan_inf_flag = 1
        return g, nan_inf_flag
    
    nan_inf_flag = 0
    return g, nan_inf_flag

def gravik(alpha, beta, betad):
    ginit = linalg.inv(alpha)
    g = ginit
    eps = 1
    it = 1
    while eps > 0.0001:
        it = it + 1
        S = g.copy()
        # print g
        # debug_here()
        g = alpha - dot(dot(beta, S), betad)
        try:
            g = linalg.inv(g)
        except LinAlgError:
            # print g
            pass

        g = dot(g, 0.5) + dot(S, 0.5)
        # eps = sum(sum(abs(g-S))) / sum(sum(abs(g+S)))
        eps = abs(g-S).sum() / abs(g+S).sum()
        # print "it  = " + str(it)
        # print "eps = " + str(eps)
        if it>200:
            #if eps > 0.01:
            #    debug_here()
            eps = -eps
    return g


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
            else:
                H[k,l]=-0
    # savetxt('pythonHam.txt', H, delimiter='\t', fmt='%f')
    return H

def Transmission(H):
    Hon=H[1:10,1:10]
    Hoff=H[1:10,11:20]
    Hoffd = Hoff.conj().T
    # Hoffd = Hoff.T
	
    eta = 0.001

    I = eye(Hon.shape[0])
    Ne = 101
    E = linspace(-2, 2, Ne)
    
    T = zeros((Ne))
    for kp in xrange(0, Ne):
        EE = E[kp]
        
        a = dot(dot((EE+1j), eta), I)-Hon
        b = dot(dot((EE+1j), eta), I)-Hoff
        bd = dot(dot((EE+1j), eta), I)-Hoffd
        
        (g1, nan_inf_flag) = grmagnus(a,bd,b)
        (g2, nan_inf_flag_2) = grmagnus(a,b,bd)
        
        sig1 = dot(dot(bd, g1), b)
        sig2 = dot(dot(b, g1), bd)
        
        gam1 = dot(1j, (sig1 - sig1.conj().T))
        gam2 = dot(1j, (sig2 - sig2.conj().T))
        
        G = linalg.inv(dot((EE + 1j) * eta, I)-Hon-sig1-sig2)
        Tr = dot(dot(dot(gam1, G), gam2), G.conj().T)
        
        T[kp] = Tr.trace(offset=0).real
    # print Hon
    # print Hoff
    # print Hoffd
    plot(E,T)
    grid(True)
    xlabel('Energy (eV)')
    ylabel('Density of States (DoS)')
    title('Density of States')
    show()
    
def main():
	
    # Load coordinate file as array
    gc = loadtxt("graphenecoordinates.txt")
    gc.view('i8,i8,i8').sort(order=['f0','f1'], axis=0) # For 64-bit systems - for 32-bit, change 'i8' to 'i4'

    # Run functions
    # plotgraphene(gc)
    # print Hamiltonian(gc)
    Transmission(Hamiltonian(gc)) #Just running Hamiltonian() twice for now
    return 0

if __name__ == '__main__':
	main()

