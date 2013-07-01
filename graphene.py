#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
graphene.py - Atomic Coordinate Generator, Hamiltonian Generator, and Transmission Calculator for graphene lattices

Trevin Gandhi, based off work by Frank Tseng
"""

from __future__ import division     # Use Python3 division
import string                       # May or may not be here for a reason
import sys                          # General python import
from pylab import *                 # Matplotlib
import timeit                       # Used for timing functions - see randomcode()
import scipy.linalg                 # Used for linear algebra operations
import scipy.sparse as sparse       # Used for sparse matrices
import scipy.sparse.linalg as spla  # Used for sparse matrix calculations
import scipy.io as sio              # Used for saving arrays as MATLAB files
import gc as garcol                 # Used for garbage collection to help save memory

#
# Getting Started
#

# Before using this program, you must install:
#   Python 2.7.x
#   Numpy python module
#   Scipy python module
#   Matplotlib python module
#   iPython (optional)
# To run, in a python shell type %run filename.py (for Windows machines)
# To change parameters of the graphene lattice, open the parameters() function and change away


#
# To-Do
#

# 1. Finalize hamiltonian and transmission calculators - MOSTLY DONE, NEED TO CHECK ANTIDOT CALCULATIONS
# 2. Build in an antidot generator - DONE
# 3. Possibly use sparse matrix calculations for transmission - MOSTLY DONE, STILL A WIP
# 4. Possibly use Cython to optimize calculations and improve calculation speed - NOT STARTED
# 5. LOW PRIORITY: Define xdist, ydist, and zdist for zigzag orientation - and make sure generator works for zigzag - NOT STARTED
# 6. LOW PRIORITY: Finish translator() - NOT STARTED

# Speed boost options: PyPy, numexpr, scipy sparse matrix, pytables, cython
# http://technicaldiscovery.blogspot.com/2011/06/speeding-up-python-numpy-cython-and.html
# http://www.physics.udel.edu/~bnikolic/teaching/phys824/MATLAB/
# http://deeplearning.net/software/theano/
# "Another way of computing inverses involves gram-schmidt orthogonalization and then transposing the matrix, the transpose of an orthogonalized matrix is its inverse!"

def parameters():
    #
    # Graphene Lattice Constants
    #

    global a, armchair, wleg, dwleg, hleg, dhleg, nanometers, x, y, numxtrans, numytrans, cuttype, rectx, recty, recth, rectw, xdist, ydist, zdist

    armchair = True      # True if armchair orientation, false if zigzag
    a = 1.42             # Carbon-carbon bond length in Angstroms (1.42 is an average of carbon single bonds (C-C) and carbon double bonds (C=C)

    if armchair:
        wleg = 0.71      # Horizontal leg: Hypotenuse (1.42 Angstroms) / 2 because it's a 30-60-90 triangle
        dwleg = 2.84     # wleg * 4 (width of hexagon)
        hleg = 1.2306    # Vertical leg: wleg * sqrt(3) (for the other leg of the triangle)
        dhleg = 2.4612   # hleg * 2  (height of hexagon)
        xdist = 2.13     # 3a/2
        ydist = hleg     # (sqrt(3)*a)/2
        zdist = a        # a
    else:
        hleg = 0.71      # Vertical leg: Hypotenuse (1.42 Angstroms) / 2 because it's a 30-60-90 triangle
        dhleg = 2.84     # hleg * 4 (width of hexagon)
        wleg = 1.2306    # Horizontal leg: hleg * sqrt(3) (for the other leg of the triangle)
        dwleg = 2.4612   # wleg * 2  (height of hexagon)

    #
    # Units
    #

    nanometers = True   # True if parameter units are in nanometers, false if in Angstroms

    #
    # General Lattice Paramters
    #

    x = 15           # Width of the unit cell
    y = 15           # Height of the unit cell
    numxtrans = 1    # Number of times to translate unit cell along the x-axis
    numytrans = 2    # Number of times to translate unit cell along the y-axis
    cuttype = 1      # 0 if no antidots, 1 if rectangular

    #
    # Rectangular Antidot Parameters
    #

    rectx = 5        # x-coordinate of the bottom left corner of the antidot
    recty = 5        # y-coordinate of the bottom left corner of the antidot
    recth = 5        # Height of the antidot
    rectw = 5        # Width of the antidot

    return 0

def grmagnus(alpha, beta, betad, kp):
    # From J. Phys. F Vol 14, 1984, 1205, M P Lopez Sancho, J. Rubio
    # 20-50 % faster than gravik
    # From Huckel IV Simulator
    tmp = linalg.inv(alpha)                     # Inverse part of Eq. 8
    t = (-1*tmp)*(betad)                        # Eq. 8 (t0)
    tt = (-1*tmp)*(beta)                        # Eq. 8 (t0 tilde)
    T = t.copy()                                # First term in Eq. 16
    Id = eye(alpha.shape[0])                    # Save the identity matrix
    Toldt = Id.copy()                           # Product of tilde t in subsequent terms in Eq. 16
    change = 1                                  # Convergence measure
    counter = 0                                 # Just to make sure no infinite loop

    etag = 0.000001
    etan = 0.0000000000001
    while linalg.norm(change) > etag and counter < 100:
        counter += 1
        Toldt = Toldt.dot(tt) # Product of tilde t in subsequent terms in Eq. 16
        if (1/(linalg.cond(Id - dot(t,tt) - dot(tt,t)))) < etan:
            g = 0
            print "1: tmp NaN or Inf occured, return forced. Kp: " + str(kp)
            return g

        tmp = linalg.inv(Id - t.dot(tt) - tt.dot(t)) # Inverse part of Eq. 12

        t = tmp.dot(t).dot(t)       # Eq. 12 (t_i)
        tt = tmp.dot(tt).dot(tt)    # Eq. 12 (t_i tilde)
        change = Toldt.dot(t)       # Next term of Eq. 16
        T = T + change              # Add it to T, Eq. 16

        if isnan(change).sum() or isinf(change).sum():
            g = 0
            print "2: tmp NaN or Inf occured, return forced. Kp: " + str(kp)
            return g

    if (1/(linalg.cond(alpha + beta.dot(T)))) < etan:
        g = 0
        print "3: tmp NaN or Inf occured, return forced. Kp: " + str(kp)
        return g

    g = linalg.inv(alpha + beta.dot(T))

    gn = abs(g - linalg.inv(alpha - beta.dot(g).dot(betad)))

    if gn.max() > 0.001 or counter > 99:
        g = 0
        print "4: Attention! not correct sgf. Kp: " + str(kp)
        return g

    nan_inf_flag = 0

    # Help save memory
    del tmp, t, tt, T, Id, Toldt, change, counter, etag, etan, gn

    return g

def grmagnus2(alpha, beta, betad, kp):
    # grmagnus with sparse matrices
    # From J. Phys. F Vol 14, 1984, 1205, M P Lopez Sancho, J. Rubio
    # 20-50 % faster than gravik
    # From Huckel IV Simulator
    global I
    tmp = linalg.inv(alpha.todense())           # Inverse part of Eq. 8
    t = (-1*tmp)*(betad).todense()              # Eq. 8 (t0)
    tt = (-1*tmp)*(beta).todense()              # Eq. 8 (t0 tilde)
    T = t.copy()                                # First term in Eq. 16
    Toldt = I.copy()                            # Product of tilde t in subsequent terms in Eq. 16
    change = 1                                  # Convergence measure
    counter = 0                                 # Just to make sure no infinite loop

    etag = 0.000001
    etan = 0.0000000000001
    while linalg.norm(change) > etag and counter < 100:
        counter += 1
        Toldt = Toldt*tt # Product of tilde t in subsequent terms in Eq. 16
        tmp = I - t*tt - tt*t
        if (1/(linalg.cond(tmp))) < etan:
            g = 0
            print "1: tmp NaN or Inf occured, return forced. Kp: " + str(kp)
            return g

        tmp = faster_inverse(tmp)                      # Inverse part of Eq. 12
        t = (tmp*t*t)                                  # Eq. 12 (t_i)
        tt = (tmp*tt*tt)                               # Eq. 12 (t_i tilde)
        change = Toldt*t                               # Next term of Eq. 16
        T = T + change                                 # Add it to T, Eq. 16

        if isnan(change).sum() or isinf(change).sum():
            g = 0
            print "2: tmp NaN or Inf occured, return forced. Kp: " + str(kp)
            return g

    g = (alpha + beta*T)

    if (1/(linalg.cond(g))) < etan:
        g = 0
        print "3: tmp NaN or Inf occured, return forced. Kp: " + str(kp)
        return g

    g = faster_inverse(g)
    gn = (abs(g - linalg.inv(alpha - beta*g*betad)))

    if gn.max() > 0.001 or counter > 99:
        g = 0
        print "4: Attention! not correct sgf. Kp: " + str(kp)
        return g

    # Help save memory
    del tmp, t, tt, T, Toldt, change, counter, etag, etan, gn

    return g

def gravik(alpha, beta, betad):
    # From J. Phys. F Vol 14, 1984, 1205, M P Lopez Sancho, J. Rubio
    # From Huckel IV Simulator
    ginit = alpha.I
    g = ginit.copy()
    eps = 1
    it = 1
    while eps > 0.000001:
        it = it + 1
        S = g.copy()
        g = alpha - beta*S*betad
        try:
            g = g.I
        except LinAlgError:
            pass

        g = g*0.5 + S*0.5
        eps = (abs(g-S).sum())/(abs(g+S).sum())
        if it>200:
            #if eps > 0.01:
            #    debug_here()
            eps = -eps
    return g

def plotgraphene(coord2,plotb):
    # Plot the graphene sheet
    if plotb:
        # Plot xy coordinates
        plot(coord2[:,0],coord2[:,1],marker='o') # plot() for line graph, scatter() for point graph
        grid(True)
        xlabel('Length (Angstroms)')
        ylabel('Width (Angstroms)')
        title('Graphene Lattice')
        show()

def hamiltonian(coord):
    # DEPRECATED

    Norb = coord.shape[0]
    numH = Norb

    x = copy(coord[:,0])
    x = vertconvert(x)

    y = copy(coord[:,1])
    y = vertconvert(y)

    z = copy(coord[:,2])
    z = vertconvert(z)

    x = tile(x, Norb)
    y = tile(y, Norb)
    z = tile(z, Norb)

    x = square(x.T-x)
    y = square(y.T-y)
    z = square(z.T-z)
    Q = x + y + z

    # Help save memory
    del x, y, z

    R = sqrt(Q)

    # Help save memory
    del Q
    garcol.collect()

    H = zeros(shape=(numH, numH))
    for k in xrange(0, (numH)):
        for l in xrange(0, (numH)):
            if R[k,l] <= 1.45 and R[k,l] >= 1.38:
                H[k,l]=-2.7
            else:
                H[k,l]=-0

    savetxt('pythonHam.txt', H, delimiter='\t', fmt='%f')
    atoms = H.shape[0]
    garcol.collect()
    return H, atoms

def transmission(H, atoms):
    global I, Ic
    # Calculate the transmission across the graphene sheet using a recursive Non-Equilibrium Green's Function
    atomsh = atoms // 2
    if atoms%2 != 0:
        atoms -= 1
    print str(atoms)
    print str(atomsh)

    # Make Hon and Hoff each half of H
    H = asmatrix(H)                        # Convert H to a matrix
    Hon = sparse.dia_matrix(H[0:atomsh,0:atomsh])
    Hoff = sparse.dia_matrix(H[0:atomsh,(atomsh):atoms])
    Hoffd = sparse.dia_matrix(Hoff.H)      # Conjugate Transpose of Hoff
    del H

    eta = -0.003
    etao = 0.001

    I = eye(Hon.shape[0])     # I is an identity matrix that is the size of Hon
    Ic = eye(Hon.shape[0], dtype=cfloat)


    Ne = 5        # Number of data points

    E = linspace(-2, 2, Ne)     # Energy Levels to calculate transmission at

    T = [None] * Ne     # Initialize T

    for kp in xrange(0, Ne):
        EE = E[kp]
        print str(EE)

        alpha = sparse.coo_matrix((EE + 1j*eta) * I - Hon)
        beta = sparse.coo_matrix((EE + 1j*eta) * I - Hoff)
        betad = sparse.coo_matrix((EE + 1j*eta) * I - Hoffd)

        # Use grmagnus
        g1 = grmagnus2(alpha,betad,beta,E[kp])
        g2 = grmagnus2(alpha,beta,betad,E[kp])

        # Use gravik
        # g1 = gravik(alpha,betad,beta)
        # g2 = gravik(alpha,beta,betad)

        #
        # Equations Used
        #

        # Non-Equilibrium Green's Function: G = [EI - H - Sig1 - Sig2]^-1
        #   EI = 0.003i
        # Transmission: T = Trace[Gam1*G*Gam2*G.H]
        # Gam1 (Broadening Function of lead 1) = i(Sig1 - Sig1.H)
        # Gam2 (Broadening Function of lead 2) = i(Sig2 - Sig2.H)

        sig1 = betad*g1*beta
        sig2 = beta*g2*betad

        # Help save memory
        del alpha, beta, betad, g1, g2

        gam1 = (1j * (sig1 - sig1.H))
        gam2 = (1j * (sig2 - sig2.H))

        G = ((EE - 1j*0.003)*I - Hon - sig1 - sig2)
        n_eq = G.shape[1]
        results = linalg.lapack_lite.zgesv(n_eq, G.shape[0], G, n_eq, zeros(n_eq, intc), Ic, n_eq, 0)
        G = asmatrix(Ic)
        Ic = eye(Hon.shape[0], dtype=cfloat)

        # Help save memory
        del sig1, sig2, n_eq

        T[kp] = trace(gam1*G*gam2*(G.H)).real

        # Help save memory
        del gam1, gam2, G
    data = column_stack((E,T))
    savetxt('pythonData.txt', data, delimiter='\t', fmt='%f')
    plot(E,T)
    grid(True)
    xlabel('Energy (eV)')
    ylabel('Transmission')
    title('Transmission vs Energy')
    show()

def maingenerator():
    # Initialize parameters
    parameters()
    global x, y, xlimit, ylimit, xtimes, ytimes, xdist, ydist, zdist

    # If x and y are in nanometers, convert to angstroms
    if nanometers:
        x *= 10
        y *= 10

    # Limits for the loops
    # xlimit
    if x > dwleg:
        xdiff = (x%dwleg)
    xlimit = x # - xdiff

    # ylimit
    if y > dhleg:
        ydiff = (y%dhleg)
    ylimit = y # - ydiff

    # coord = generatorloop()  - Deprecated

    xtimes = round(xlimit // xdist)
    ytimes = round(ylimit / ydist)

    print "xlimit:" + str(xlimit)
    print "ylimit:" + str(ylimit)
    if xdiff:
        print "xdiff:" + str(xdiff)
    if ydiff:
        print "ydiff:" + str(ydiff)
    print "xdist:" + str(xdist)
    print "ydist:" + str(ydist)
    print "Xtimes:" + str(xtimes)
    print "Ytimes:" + str(ytimes)
    (coord, coord2) = vgenerator()

    return coord, coord2

def generatorloop():
    # DEPRECATED

    # Initializations
    global rectx, recty, recth, rectw
    pointy = 1
    holey = -1
    coord = array([])
    checker = 0
    hunits = 0
    counter = 1

    if cuttype == 1:
        # Convert to Angstroms
        if nanometers:
            rectx *= 10
            recty *= 10
            recth *= 10
            rectw *= 10

        # Get upper left Y value and bottom right x value of rectangle
        oppx = rectx + rectw
        oppy = recty + recth

    while hunits <= ylimit:
        if pointy > 0:
            wunits = wleg
        else:
            wunits = 0

        while wunits <= xlimit:
            if cuttype == 1:
                if (((hunits >= recty) and (hunits <= oppy)) and ((wunits >= rectx) and (wunits <= oppx))):
                    cut = True
                else:
                    cut = False
            else:
                cut = False

            if not cut:
                try:
                    coord = vstack((coord,[wunits,hunits,0]))
                except NameError:
                    coord = [wunits,hunits,0]

            if (checker == 0) and (pointy == 1):
                holey = -1
                checker = 1
            elif (checker == 0) and (pointy == -1):
                holey = 1
                checker = 1
            else:
                holey *= -1

            if (holey > 0):
                wincrement = wleg*4
            else:
                wincrement = wleg*2

            # In C code, had if (armchair) and else { wincrement = dwleg} - haven't added that in yet
            wunits += wincrement

        hincrement = hleg
        # In C code, had if (armchair) and else ifs - haven't added that yet

        if (hunits == 0):
            pointy = -1
        else:
            pointy += sign(pointy)
            if abs(pointy) > 1:
                pointy = -sign(pointy)

        checker = 0

        hunits += hincrement

    savetxt('graphenecoordinates.txt', coord, delimiter='\t', fmt='%f')

    return coord

def vgenerator():
    # Creates a 3-dimensional array (n x 0 x 2 - since the index starts at 0) 'coord'
    # Where n is the number of atoms in the sheet
    # Defined by the unit cell consisting of two atoms, one a horizontal translation of the other
    # Each atomic coordinate is defined by 3 numbers:
    #   A x-value (coordx): As the unit cell is translated horizontally, the x-value increments by 1
    #   A y-value (coordy): As the unit cell is translated vertically, the y-value increments by 1
    #   A unit cell value (coordu): Defines which point in the unit cell the atom is
    #       0 = the atom on the left
    #       1 = the atom on the right
    # To convert these numbers into xyz coordinates,
    #   If a = 1.42
    #   ((coordx, coordy, coordu).((3a/2), 0, a), (coordx, coordy, coordu).(0, (sqrt(3)*a)/2, 0))
    #   Where . represents the dot product

    global xtimes, ytimes, rectx, recty, rectw, recth, cuttype

    #
    # Coordinate Generator
    #

    # Build upwards, then horizontally
    for j in xrange(0,int(xtimes)):
        for i in xrange(0,int(ytimes),2):
            try:
                if j%2 != 0:
                    coord = vstack((coord,[[[j,i+1,0]]]))
                    coord = vstack((coord,[[[j,i+1,1]]]))
                else:
                    coord = vstack((coord,[[[j,i,0]]]))
                    coord = vstack((coord,[[[j,i,1]]]))
            except NameError:
                coord = [[[j,i,0]]]
                coord = vstack((coord,[[[j,i,1]]]))

    # Build horizontally, then upwards
    # for i in xrange(0,int(ytimes),2):
    #     for j in xrange(0,int(xtimes)):
    #         try:
    #             if j%2 != 0:
    #                 coord = vstack((coord,[[[j,i+1,0]]]))
    #                 coord = vstack((coord,[[[j,i+1,1]]]))
    #             else:
    #                 coord = vstack((coord,[[[j,i,0]]]))
    #                 coord = vstack((coord,[[[j,i,1]]]))
    #         except NameError:
    #             coord = [[[j,i,0]]]
    #             coord = vstack((coord,[[[j,i,1]]]))
    #


    #
    # Antidot Generator
    #

    # Can possibly have this calculate coord2 while it is looking for antidots
    if cuttype == 1:
        # Convert to Angstroms
        if nanometers:
            rectx *= 10
            recty *= 10
            recth *= 10
            rectw *= 10

        # Get upper left Y value and bottom right x value of rectangle
        oppx = rectx + rectw
        oppy = recty + recth

    (x,y,z) = coord.shape
    a = 0
    while a < x:
        # Translate vector form of coord into xyz points of coord2
        # For xy points, remove the ", 0" from the end of the lines
        cx = coord[a,0,0]*xdist + coord[a,0,2]*zdist
        cy = coord[a,0,1]*ydist

        cut = False
        # Check to see if antidot at that location
        if (cuttype == 1) and ((cx >= rectx and cx <= oppx) and (cy >= recty and cy <= oppy)):
            coord = delete(coord, a, 0)
            (x,y,z) = coord.shape # Redefine x since coord just got shortened
            cut = True
            a -= 1 # Prevent while loop from skipping a line

        # Build coord2 - array of xyz atomic coordinates
        if not cut:
            try:
                coord2 = vstack((coord2, [cx, cy, 0]))
            except NameError:
                coord2 = [cx, cy, 0]

        # Increment a
        a += 1


    # Save as a MATLAB file for easy viewing and to compare MATLAB results with Python results
    sio.savemat('coord.mat', {'coord':coord})

    # Save xyz coordinates to graphenecoordinates.txt
    savetxt('graphenecoordinates.txt', coord2, delimiter='\t', fmt='%f')

    return coord, coord2

def vhamiltonian(coord):
    global xdist, ydist, zdist

    # Generates the Hamiltonian of the sheet
    # Uses t = -2.7 eV as the interaction (hopping) parameter
    # Only does nearest-neighbor calculations

    (x,y,z) = coord.shape
    num = x*y # Number of atoms in the lattice
    Ham = zeros(shape=(num, num))

    for i in xrange(num):
        for j in xrange(num):

            if coord[i,0,0] == coord[j,0,0] and coord[i,0,1] == coord[j,0,1] and coord[i,0,2] != coord[j,0,2]: # Check if they are part of the same translated unit cell
                Ham[i,j] = -2.7

            elif abs(coord[j,0,1]-coord[i,0,1]) == 1: # Check to see if y-coordinates are off by one
                if coord[i,0,0]-coord[j,0,0] == 1 and coord[i,0,2] != coord[j,0,2]: # Check that ix is greater than jx and that unit cell numbers aren't equal
                    if absolute((coord[i,0,0]*xdist + coord[i,0,2]*zdist)-(coord[j,0,0]*xdist + coord[j,0,2]*zdist)) == wleg: # Check that their x-distance is a/2
                        Ham[i,j] = -2.7

                elif coord[j,0,0]-coord[i,0,0] == 1 and coord[j,0,2] != coord[i,0,2]: # Check that jx is greater than ix and that unit cell numbers aren't equal
                    if absolute((coord[i,0,0]*xdist + coord[i,0,2]*zdist)-(coord[j,0,0]*xdist + coord[j,0,2]*zdist)) == wleg: # Check that their x-distance is a/2
                        Ham[i,j] = -2.7

    # Save as a MATLAB file for easy viewing and to compare MATLAB results with Python results
    sio.savemat('Ham.mat', {'Ham':Ham})

    # Help save memory
    garcol.collect()

    return Ham, num

def sign(x):
    # Returns 1 if x > 0, 0 if x == 0, -1 if x < 0
    return (( x > 0 ) - ( x < 0 ))

def vertconvert(x):
    # Convert a 1xn array to nx1
    x = atleast_2d(x)
    x = column_stack((x))
    return x

def faster_inverse(A):
    # http://stackoverflow.com/questions/11972102/is-there-a-way-to-efficiently-invert-an-array-of-matrices-with-numpy
    # http://nullege.com/codes/show/src%40n%40u%40numpy-refactor-HEAD%40numpy%40linalg%40linalg.py/18/numpy.core.zeros/python
    # http://www.netlib.org/lapack/double/dgesv.f
    # http://www.netlib.org/lapack/complex16/zgesv.f

    # Even faster inverse
    # numpy/scipy's linalg.inv(A) essentially does linalg.solve(A, identity(A.shape[0])
    # Looking into linalg.solve(), one can see that there are many safeguards to ensure the correct input
    # Removing those safeguards greatly speeds up the code

    global Ic

    b = Ic.copy()
    n_eq = A.shape[1]
    results = linalg.lapack_lite.zgesv(n_eq, A.shape[0], A, n_eq, zeros(n_eq, intc), b, n_eq, 0)
    if results['info'] > 0:
        raise LinAlgError('Singular matrix')
    return asmatrix(b)

def translator(gc, numxtrans, numytrans):
    # WIP
    # Cannot yet translate in both x and y directions

    # Translate an array gc numxtrans times horizontally and numytrans vertically

    global xlimit, ylimit

    gc1 = copy(gc[:,0])
    gc2 = copy(gc[:,1])
    gc3 = copy(gc[:,2])

    atoms = gc.shape[0]
    atomsh = atoms / 2
    if atoms%2!=0:
        atoms -= 1

    for counter in xrange(4):
        for i in xrange(numxtrans if (numxtrans > numytrans) else numytrans):
            gc1trans = copy(gc1)
            if i < numxtrans and (counter == 2 or counter == 3):
                gc1trans = [x+((xlimit+1.24)*(i+1)) for x in gc1trans]
            gc1trans = vertconvert(gc1trans)

            gc2trans = copy(gc2)
            if i < numytrans and (counter == 1 or counter == 3):
                gc2trans = [y+((ylimit+1.24)*(i+1)) for y in gc2trans]
            gc2trans = vertconvert(gc2trans)

            gc3trans = copy(gc3)
            gc3trans = vertconvert(gc3trans)

            gc1trans = column_stack((gc1trans,gc2trans))
            gc1trans = column_stack((gc1trans,gc3trans))

            gc = concatenate((gc,gc1trans),axis=0)

    return gc

def main():
    # Check to make sure garbage collection is enabled
    garcolboolean = garcol.isenabled()
    print garcolboolean

    # Plot the graphene lattice?
    plot = False

    # Generate Coordinates
    (coord, coord2) = maingenerator()

    # Plot graphene and help save memory
    plotgraphene(coord2, plot)
    del coord2

    # Generate Hamiltonian and help save memory
    (H, atoms) = vhamiltonian(coord)
    del coord

    # Calculate Transmission
    transmission(H, atoms)

    return 0

def randomcode():
    # Just random code that I want to keep that doesn't belong anywhere else

    # Coordinate plotting in MATLAB:
    #   plot(graphenecoordinates(:,1),graphenecoordinates(:,2),'o')'

    # Don't truncate printed arrays
    #   set_printoptions(threshold=nan)

    # Load coordinate file as array
    #   gc = loadtxt("graphenecoordinates.txt")
    #   gc.view('i8,i8,i8').sort(order=['f1'], axis=0) # For 64-bit systems - for 32-bit, change 'i8' to 'i4'

    # Hamiltonian Speed Testing
    #   gc = loadtxt("graphenecoordinates.txt")
    #   t1 = timeit.Timer(lambda: hamiltonian(gc))
    #   t2 = timeit.Timer(lambda: vhamiltonian(coord))
    #   print t1.timeit(number=1)
    #   print t2.timeit(number=1)

    # Slow faster version of faster_inverse

    # (Slow) Fast version
    #   lapack_routine = lapack_lite.zgesv
    #   print A.shape
    #   b = eye(A.shape[0], dtype=A.dtype)
    #   n_eq = A.shape[1]
    #   n_rhs = A.shape[0]
    #   pivots = zeros(n_eq, intc)
    #   identity = eye(n_eq)
    #   b = copy(identity)
    #   results = linalg.lapack_lite.zgesv(n_eq, n_rhs, A, n_eq, pivots, b, n_eq, 0)
    #   if results['info'] > 0:
    #       raise LinAlgError('Singular matrix')

    # For speed testing, in IPython type (after commenting out show plot of Transmission and graphene)
    #   %timeit -n 20 %run graphenex-xx.py
    #   Where 20 is the number of times you want it to loop / 3

    # Translator Remove Zeroes
    #   gc2 = filter(lambda a: a != 0, gc2) - Was used to remove zeroes, but adding (xlimit+1.24) solves that problem
    return 0

if __name__ == '__main__':
    main()
