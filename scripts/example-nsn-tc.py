#!/usr/bin/env python
"""
Compute the critical temperature of a NSN structure, comparing to analytical
results.

Perform (embarassingly) parallel computation using mpi4py
(http://mpi4py.scipy.org/).

"""
import sys, os

from scipy import *
from scipy import linalg, optimize, special

import usadel1 as u
from mpi4py import MPI

MAX_NE = 100
E_MAX = 300
OMEGA_D = 1000

def main():
    Ls = linspace(3, 16, 16)

    # Parallel computation
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nslaves = comm.Get_size() - 1

    assert nslaves > 0

    def fmt_line(L, Tc):
        return "%g %g %g" % (L, Tc, nsn_tc_an(L, OMEGA_D))

    if rank == 0:
        # Master
        fn = 'tc.dat'
        f = open(fn, 'w')
        print >> f, "% L_s, T_c"

        # Seed slaves
        msgs = []
        for j, L in enumerate(Ls):
            msg = comm.send(L, dest=1 + j % nslaves)
            msgs.append(msg)

        # Collect results
        nresults = 0
        while nresults < len(Ls):
            L, Tc = comm.recv(source=MPI.ANY_SOURCE)
            nresults += 1
            print >> f, fmt_line(L, Tc)
            f.flush()
        f.close()

        # Tell slaves to exit
        for j in xrange(nslaves):
            msg = comm.send(None, dest=1 + j)
            msgs.append(msg)
    else:
        # Slave
        while True:
            L = comm.recv(source=0)
            if L is None:
                break
            Tc = compute(L)
            print "#%d %s" % (rank, fmt_line(L, Tc))
            comm.send([L, Tc], dest=0)

def compute(L):
    # Make the overlaps a part of the superconducting wire
    geometry = nsn_geometry(L=L, T=1.0, omega_D=OMEGA_D)

    #print "\n### L_s = %g" % L
    Tc = tc_num(geometry, tol=1e-3,
                max_ne=MAX_NE, E_max=E_MAX)

    return Tc

def nsn_geometry(L, T, omega_D):
    """
    Generate a geometry describing a NSN structure

       N==========N

    where the wire is superconducting.

    """

    g = u.Geometry(nnode=2, nwire=1)

    g.t_type = [u.NODE_CLEAN_N_TERMINAL,
                u.NODE_CLEAN_N_TERMINAL]
    g.w_type = [u.WIRE_TYPE_S]

    g.w_conductance = 1.0
    g.w_length = L
    g.w_ends[0,:] = 0, 1

    g.omega_D = omega_D
    g.w_delta[0] = 1.0
    g.coupling_lambda = 1/log(2*omega_D)

    g.t_t = T

    return g


#------------------------------------------------------------------------------
# Numerical Tc computation
#------------------------------------------------------------------------------

def tc_num(geometry, tol=1e-2, matsubara=True, **kw):
    last_gtr_delta = [1.0]

    def func(T):
        if T > 0.6:
            return -1
        #print "    @ T =", T, ":"
        geometry.t_t = T
        geometry.w_delta = last_gtr_delta[0]

        it = u.self_consistent_matsubara_iteration(geometry,
                                                   output_func=None,
                                                   **kw)

        for k, d, v in it:
            rel_err = d.residual_norm() / abs(d.get()).max()
            #sys.stdout.write("%g(%.1g)  " % (rel_err, abs(d.get()).max()))
            #sys.stdout.flush()
            if rel_err < 1e-3 and v < 1e-4:
                break
            if abs(d.get()).max() < tol:
                geometry.w_delta = 0.0
                geometry.w_phase = 0.0
                break

        dmax = abs(geometry.w_delta).max()
        if dmax == 0:
            #sys.stdout.write("=> N\n")
            return -1
        last_gtr_delta[0] = geometry.w_delta.copy()
        #sys.stdout.write("=> S (%g)\n" % dmax)
        return dmax

    Tc = brenth_right(func, 0., 0.9, xtol=tol, max_pow=3)
    if isnan(Tc):
        Tc = 0.0
    return Tc


#------------------------------------------------------------------------------
# Analytical Tc of a NSN structure
#------------------------------------------------------------------------------

@vectorize
def nsn_tc_an(L_s, omega_D, tol=1e-4):
    r"""
    Compute the Tc of a NSN structure (where N are N-terminals).

    """
    if L_s == 0:
        return 0.0

    lambda_0 = 1/log(2*omega_D)
    def func(T):
        return (
            special.polygamma(0, omega_D/(2*pi*T) + .25*(pi/(T*L_s**2) + 6))
            - special.polygamma(0, .25*(pi/(T*L_s**2) + 2))
            - 1/lambda_0)
    
    Tc = brenth_right(func, 0., 1.01, xtol=tol, max_pow=8)
    if isnan(Tc):
        Tc = 0.0
    return Tc

#------------------------------------------------------------------------------
# Utility routines
#------------------------------------------------------------------------------

def brenth_right(func, a, b, xtol=1e-3, max_pow=10):
    """Find root of func(x) in [a, b], avoiding evaluation of func(a)"""
    _memo = {}

    def memo_func(x):
        x = float(x)
        if x not in _memo:
            _memo[x] = func(x)
        return _memo[x]

    v0 = memo_func(b)
    x = b
    for j in xrange(max_pow):
        last_x = x
        x = (a + x)/2.
        v = memo_func(x)
        if v*v0 < 0:
            break
    else:
        return nan

    return optimize.brenth(memo_func, x, last_x, xtol=xtol)

if __name__ == "__main__":
    main()
