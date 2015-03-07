#!/usr/bin/env python
"""
singlewire.py
-------------

Usage: Run singlewire.py and input values in the format

    NAME = VALUE

after all values have been input, type

    RUN

If you type values in by hand, use Control+D to exit.
Of course, it may be more useful to write them into a file, and do

    ./singlewire.py < input.inp > output.dat

Values that you need to input before running are

    E_T         Thouless energy [Kelvin]
    Delta       Energy gap [Kelvin]
    R           Resistance [Kelvin]
    phi         Phase difference

The following are optional:

    min_T       minimum temperature  [default = 1e-6]
    max_T       maximum temperature  [default = T_c]
    n_T         Number of temperature points to plot [default = 50]
    Delta_Tdep  Use a temperature-dependent Delta (more work)
                If 0, Delta is assumed temperature-independent.
                If not 0, Delta is assumed to have the same temperature
                dependence as in bulk.
    omega_D     Debye temperature [Kelvin]
                Needed only if Delta is T-dependent.
    n_E         Number of energy discretization points [Default=500]

The progam then prints the corresponding temperature-dependence of
supercurrent.
"""
from __future__ import division
import usadel1 as u
import sys
import re
import scipy as S
import traceback
S.pkgload('integrate', 'optimize')

__revision__ = "$Id: singlewire.py 3152 2006-09-29 10:26:22Z pauli $"

ec = 1.60217733e-19
"""Electron charge       ec = 1.60217733e-19             [C]"""
Rgas = 8.31451
"""Ideal gas constant    R = 8.31451                     [J / K * mol]"""
N_A = 6.0221367e23
"""Avogadro's number     N_A = 6.0221367e23              [1]"""
k_B = Rgas / N_A
"""Boltzmann's constant  k_B = Rgas/N_A = 1.380658e-23   [J / K]"""
EulerGamma = 0.577215664901532860606512090083
"""Euler's Gamma"""


def main():
    if len(sys.argv) > 1:
        print __import__('__main__').__doc__
        return

    allvars = ['e_t', 'delta', 'r', 'phi', 'min_t', 'max_t', 'n_t',
               'delta_tdep', 'omega_d', 'n_e']

    var = {}
    while True:
        try:
            line = raw_input()
        except EOFError:
            break
        line = re.sub(r'\s+', '', line)
        line = re.sub(r'[%#].*$', '', line)
        if not line.strip():
            continue
        vals = line.split('=')
        try:
            if len(vals) == 1 and vals[0].lower() == "run":
                try:
                    run(var)
                except Exception, e:
                    traceback.print_exc()
                    #print >> sys.stderr, "Error while running: %s" % e
            elif len(vals) == 2:
                vals[0] = vals[0].lower()
                if vals[0] not in allvars:
                    raise RuntimeError("unknown variable %s" % vals[0])
                var[vals[0]] = float(vals[1])
            else:
                raise RuntimeError("not understood")
        except Exception, e:
            print >> sys.stderr, "Invalid input '%s': %s" % (line, e)

def run(var):

    ## Collect parameters
    
    if not 'e_t' in var:
        raise RuntimeError("Thouless energy not given")
    elif not 'delta' in var:
        raise RuntimeError("Delta not given")
    elif not 'phi' in var:
        raise RuntimeError("phi not given")
    elif not 'r' in var:
        raise RuntimeError("R not given")

    phi     = var['phi']
    E_T     = var['e_t']
    R       = var['r']
    Delta   = var['delta']

    min_T = var.get('min_t', 1e-5)
    max_T = var.get('max_t', 1.02 * S.exp(EulerGamma)/S.pi*Delta)
    n_T   = var.get('n_t', 50)

    n_E   = var.get('n_e', 500)

    Delta_Tdep = bool(var.get('delta_tdep', False))

    if Delta_Tdep and not 'omega_d' in var:
        raise RuntimeError("Omega_D not given")

    omega_D = var.get('omega_d', 394)

    ## Run

    g = geometry_SNS(-phi, Delta/E_T)

    print "% Quasiclassical calculation for I_S in a single wire."
    print "% :author: Pauli Virtanen <pauli@ltl.tkk.fi>, 2006"
    print "% :organization: Low Temperature Laboratory, TKK"
    print "%"
    print '%% Calculating temperature range T = %s ... %s K' % (min_T, max_T)
    print "%% R = %g Ohm" % R
    print "%% E_T = %g K" % E_T
    print "%% Phase difference = %g" % phi
    print "%% Delta = %g K" % Delta
    if Delta_Tdep:
        print '%% Omega_D = %g' % omega_D
        print '% Delta is temperature-dependent'
    else:
        print '% Delta is not temperature-dependent'
    print "%"
    print "%", ("%15s  " * 7) % ('T (K)', 'Delta (K)', 'I_S (A)',
                                 'phi', 'R (Ohm)', 'E_T (K)', 'omega_D (K)')
    
    solver = u.CurrentSolver(g, ne=n_E, chunksize=300)

    T_in_Ks = S.linspace(min_T, max_T, n_T)

    for i, T_in_K in enumerate(T_in_Ks):
        T = T_in_K / E_T
        
        if Delta_Tdep:
            print >> sys.stderr, "%"
            print >> sys.stderr, "%% T = %g K  (%d / %d)" % (T_in_K, i,
                                                             len(T_in_Ks))
            g.t_delta = get_bulk_Delta(omega_D, Delta, T_in_K) / E_T
            print >> sys.stderr, "%% Delta = %g K" % (g.t_delta[0] * E_T)
            if g.t_delta[0] > 1e-4:
                solver.solve_spectral()
            sys.stdout.flush()
        else:
            solver.solve_spectral_if_needed(calculate_G=False)

        if g.t_delta[0] > 1e-4:
            E = solver.E
            Ic = S.integrate.trapz(
                E, solver.coefficient.ijE[:,0,0] * S.tanh(.5*E/T))
        else:
            Ic = 0

        print " ", ("%15.7g  " * 7) % (T_in_K, g.t_delta[0]*E_T,
                                       Ic*E_T*k_B/ec/R,
                                       phi, R, E_T, omega_D)

def geometry_SNS(phi, Delta):
    """N-wire between S terminals."""
    g = u.Geometry(1, 2)
    
    g.t_type = [ u.NODE_CLEAN_S_TERMINAL,
                 u.NODE_CLEAN_S_TERMINAL ]
    
    g.t_delta = [ Delta, Delta ]
    g.t_phase = [ -.5*phi, .5*phi ]
    
    g.t_inelastic = 1e-9
    g.t_spinflip = 0
    g.t_t = 1e-6
    g.t_mu = 0
    
    g.w_type = u.WIRE_TYPE_N
    g.w_length = 1
    g.w_conductance = 1
    g.w_inelastic = 1e-9
    g.w_spinflip = 0
    
    g.w_ends[0,:] = [ 0, 1 ]
    
    return g

def get_bulk_Delta(omega_D, Delta_0, T0, initial_guess=None):
    oomega_D = omega_D
    Delta_0 *= 1000/omega_D
    T0 *= 1000/omega_D
    omega_D = 1000

    T0 += 1e-15

    if not initial_guess:
        initial_guess = Delta_0
    
    c = 1/S.arccosh(omega_D / Delta_0)

    def func(z):
        Delta = z[0]
        def integrand(E):
            return Delta*S.tanh(.5*E/T0) / S.sqrt(E**2 - Delta**2 + 1e-9)
        r = S.integrate.quad(integrand, Delta, omega_D,
                             epsrel=1e-5, epsabs=1e-5*Delta/c + 1e-8)
        x = r[0] - Delta/c
        if abs(x) < 1e-40: return 0
        return x

    z = S.optimize.fsolve(func, initial_guess, xtol=1e-4) * oomega_D/1000
    return abs(z)

if __name__ == "__main__": main()
