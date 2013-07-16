def random_code():
    # Just random code that I want to keep that doesn't belong anywhere else

    # Coordinate plotting in MATLAB:
    #   plot(graphenecoordinates(:,1),graphenecoordinates(:,2),'o')'

    # Imports
    # import string # May or may not be here for a reason
    # import sys    # General python import
    # import timeit # Used for timing functions - see random_code()
    # import scipy.sparse.linalg as spla  # Used for sparse matrix calculations

    # Don't truncate printed arrays
    #   set_printoptions(threshold=nan)

    # DeMorgan's Law on a line in generator_loop()
    # if not (not (cut_type == 1) or not (rect_y <= h_units <= opp_y) or not (rect_x <= w_units <= opp_x)):

    # Load coordinate file as array
    #   gc = loadtxt("graphenecoordinates.txt")
    #   gc.view('i8,i8,i8').sort(order=['f1'], axis=0) # For 64-bit systems - for 32-bit, change 'i8' to 'i4'

    # Hamiltonian Speed Testing
    #   gc = loadtxt("graphenecoordinates.txt")
    #   t1 = timeit.Timer(lambda: hamiltonian(gc))
    #   t2 = timeit.Timer(lambda: v_hamiltonian(coord))
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
    #   gc2 = filter(lambda a: a != 0, gc2) - Was used to remove zeroes, but adding (x_limit+1.24) solves that problem

    # # noinspection PyArgumentList
    # def hamiltonian(coord):
    #     """
    #     DEPRECATED
    #     """
    #     Norb = coord.shape[0]
    #     numH = Norb
    #
    #     x = np.copy(coord[:, 0])
    #     x = vert_convert(x)
    #
    #     y = np.copy(coord[:, 1])
    #     y = vert_convert(y)
    #
    #     z = np.copy(coord[:, 2])
    #     z = vert_convert(z)
    #
    #     x = np.tile(x, Norb)
    #     y = np.tile(y, Norb)
    #     z = np.tile(z, Norb)
    #
    #     x = np.square(x.T - x)
    #     y = np.square(y.T - y)
    #     z = np.square(z.T - z)
    #     Q = x + y + z
    #
    #     # Help save memory
    #     del x, y, z
    #
    #     R = np.sqrt(Q)
    #
    #     # Help save memory
    #     del Q
    #     garcol.collect()
    #
    #     H = np.zeros(shape=(numH, numH))
    #     for k in xrange(0, numH):
    #         for l in xrange(0, numH):
    #             if 1.38 <= R[k, l] <= 1.45:
    #                 H[k, l] = -2.7
    #             else:
    #                 H[k, l] = 0
    #
    #     np.savetxt('pythonHam.txt', H, delimiter='\t', fmt='%f')
    #     atoms = H.shape[0]
    #     garcol.collect()
    #     return H, atoms

    # def generator_loop():
    #     """
    #     DEPRECATED
    #     """
    #
    #     global rect_x, rect_y, rect_h, rect_w
    #     pointy = 1
    #     holey = -1
    #     coord = linalg.array([])
    #     checker = 0
    #     h_units = 0
    #
    #     if cut_type == 1:
    #         # Convert to Angstroms
    #         if nanometers:
    #             rect_x *= 10
    #             rect_y *= 10
    #             rect_h *= 10
    #             rect_w *= 10
    #
    #         # Get upper left Y value and bottom right x value of rectangle
    #         opp_x = rect_x + rect_w
    #         opp_y = rect_y + rect_h
    #
    #     while h_units <= y_limit:
    #         if pointy > 0:
    #             w_units = w_leg
    #         else:
    #             w_units = 0
    #
    #         while w_units <= x_limit:
    #             if cut_type == 1 and rect_y <= h_units <= opp_y and rect_x <= w_units <= opp_x:
    #                 cut = True
    #             else:
    #                 cut = False
    #
    #             if not cut:
    #                 try:
    #                     coord = np.vstack((coord, [w_units, h_units, 0]))
    #                 except NameError:
    #                     coord = [w_units, h_units, 0]
    #
    #             if (checker == 0) and (pointy == 1):
    #                 holey = -1
    #                 checker = 1
    #             elif (checker == 0) and (pointy == -1):
    #                 holey = 1
    #                 checker = 1
    #             else:
    #                 holey *= -1
    #
    #             if holey:
    #                 w_increment = w_leg * 4
    #             else:
    #                 w_increment = w_leg * 2
    #
    #             # In C code, had if (armchair) and else {w_increment = dw_leg} - haven't added that in yet
    #             w_units += w_increment
    #
    #         h_increment = h_leg
    #         # In C code, had if (armchair) and else ifs - haven't added that yet
    #
    #         if not h_units:
    #             pointy = -1
    #         else:
    #             pointy += sign(pointy)
    #             if abs(pointy) > 1:
    #                 pointy = -sign(pointy)
    #
    #         checker = 0
    #
    #         h_units += h_increment
    #
    #     np.savetxt('graphenecoordinates.txt', coord, delimiter='\t', fmt='%f')
    #
    #     return coord

    # def grmagnus(alpha, beta, betad, kp):
    #     """
    #     From J. Phys. F Vol 14, 1984, 1205, M P Lopez Sancho, J. Rubio
    #     20-50 % faster than gravik
    #     From Huckel IV Simulator
    #     """
    #     tmp = linalg.inv(alpha)                    # Inverse part of Eq. 8
    #     t = (-1 * tmp) * betad              # Eq. 8 (t0)
    #     tt = (-1 * tmp) * beta              # Eq. 8 (t0 tilde)
    #     T = t.copy()                        # First term in Eq. 16
    #     Id = np.eye(alpha.shape[0])            # Save the identity matrix
    #     Toldt = Id.copy()                   # Product of tilde t in subsequent terms in Eq. 16
    #     change = 1                          # Convergence measure
    #     counter = 0                         # Just to make sure no infinite loop
    #
    #     etag = 0.000001
    #     etan = 0.0000000000001
    #     while linalg.norm(change) > etag and counter < 100:
    #         counter += 1
    #         Toldt = linalg.dot(Toldt, tt)    # Product of tilde t in subsequent terms in Eq. 16
    #         if (1 / (np.linalg.cond(Id - linalg.dot(t, tt) - linalg.dot(tt, t)))) < etan:
    #             g = 0
    #             print "1: tmp NaN or Inf occurred, return forced. Kp: " + str(kp)
    #             return g
    #
    #         tmp = linalg.inv(Id - t.dot(tt) - tt.dot(t))    # Inverse part of Eq. 12
    #
    #         t = tmp.dot(t).dot(t)       # Eq. 12 (t_i)
    #         tt = tmp.dot(tt).dot(tt)    # Eq. 12 (t_i tilde)
    #         change = Toldt.dot(t)       # Next term of Eq. 16
    #         T = T + change              # Add it to T, Eq. 16
    #
    #         if np.isnan(change).sum() or np.isinf(change).sum():
    #             g = 0
    #             print "2: tmp NaN or Inf occurred, return forced. Kp: " + str(kp)
    #             return g
    #
    #     if (1 / (np.linalg.cond(alpha + beta.dot(T)))) < etan:
    #         g = 0
    #         print "3: tmp NaN or Inf occurred, return forced. Kp: " + str(kp)
    #         return g
    #
    #     g = linalg.inv(alpha + beta.dot(T))
    #
    #     gn = abs(g - linalg.inv(alpha - beta.dot(g).dot(betad)))
    #
    #     if gn.max() > 0.001 or counter > 99:
    #         g = 0
    #         print "4: Attention! not correct sgf. Kp: " + str(kp)
    #         return g
    #
    #     # Help save memory
    #     del tmp, t, tt, T, Id, Toldt, change, counter, etag, etan, gn
    #
    #     return g
    pass
