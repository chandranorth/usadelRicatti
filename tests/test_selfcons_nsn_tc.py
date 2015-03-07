from testutils import *

from scipy import *
from scipy import linalg, optimize, special

import usadel1 as u

def test_single_one():
    params = [(4.0, 10000, False)]
    return test_single(params)

@slow
def test_single(params=None):
    if params is None:
        params = [(8.0,   300, False),
                  (8.0, 10000, False),
                  (8.0, 10000, True),
                  (4.0, 10000, False),
                  (3.3, 10000, False),
                  (3.3, 10000, True),
                  (3.0, 10000, False),
                  ]

    for L_s, omega_D, force_integral in params:
        geometry = nsn_geometry(1.0, L_s, omega_D)
        an = nsn_tc(L_s, omega_D, tol=1e-5)
        num = nsn_tc_num(geometry, tol=1e-3, force_integral=force_integral)
        print "==", an, "vs", num
        assert allclose(an, num, atol=1e-2), (L_s, omega_D, an, num)

@slow
def test_triple():
    params = [(8.0, 1.0, 0.3,  300,    False),
              (4.0, 1.0, 1.0,  1000,   False),
              (8.0, 0.5, 2.0,  10000,  False),
              (8.0, 0.1, 3.0,  10000,  False),
              (4.0, 10.0, 0.1, 10000,  False),
              (3.3, 3.0, 0.1,  10000,  False),
              (3.3, 1.0, 0.5,  10000,  False),
              (3.1, 0.1, 0.05, 10000,  False),
              (0.1, 3.5, 0.05, 10000,  False),
              ]

    for L_s1, L_s2, L_s3, omega_D, force_integral in params:
        geometry = nsn_geometry_triple(1.0, L_s1, L_s2, L_s3, omega_D)
        an = nsn_tc(L_s1 + L_s2 + L_s3, omega_D, tol=1e-5)
        num = nsn_tc_num(geometry, tol=1e-3,
                         max_ne=150,
                         E_max=500/max([L_s1,L_s2,L_s3])**2,
                         force_integral=force_integral)
        print "==", an, "vs", num
        assert allclose(an, num, atol=1e-2), (L_s1, L_s2, L_s3,
                                              omega_D, an, num)

#------------------------------------------------------------------------------
# Numerical Tc of a NSN structure
#------------------------------------------------------------------------------

def nsn_tc_num(geometry, tol=1e-2, matsubara=True, **kw):
    last_gtr_delta = [1.0]

    def func(T):
        print "    @ T =", T, ":"
        geometry.t_t = T
        geometry.w_delta = last_gtr_delta[0]

        if matsubara:
            it = u.self_consistent_matsubara_iteration(geometry,
                                                       output_func=None,
                                                       **kw)
        else:
            solver = u.CurrentSolver(geometry)
            it = u.self_consistent_realtime_iteration(solver, output_func=None,
                                                      **kw)

        zero_count = 0
        for k, d, v in it:
            rel_err = d.residual_norm() / abs(d.get()).max()
            sys.stdout.write("%g(%.2g)  " % (rel_err, abs(d.get()).max()))
            sys.stdout.flush()
            if abs(d.get()).max() < 1e-3:
                zero_count += 1
            else:
                zero_count = 0
            if zero_count > 4:
                geometry.w_delta = 0.0
                geometry.w_phase = 0.0
                break
            if rel_err < 1e-4 and v < 1e-4:
                break

        dmax = abs(geometry.w_delta).max()
        if dmax == 0:
            sys.stdout.write("=> N\n")
            return -1
        last_gtr_delta[0] = geometry.w_delta.copy()
        sys.stdout.write("=> S (%g)\n" % dmax)
        return dmax

    Tc = brenth_right(func, 0., 0.8, xtol=tol, max_pow=5)
    if isnan(Tc):
        Tc = 0.0
    return Tc

def nsn_geometry(T, L_s, omega_D):
    g = u.Geometry(nnode=2, nwire=1)

    g.t_type = [u.NODE_CLEAN_N_TERMINAL,
                u.NODE_CLEAN_N_TERMINAL]
    g.w_type = [u.WIRE_TYPE_S]

    g.w_conductance = 1.0

    g.w_ends[0,:] = 0, 1
    g.w_length = L_s

    Delta_0 = 1.0
    omega_D = omega_D
    lambda_0 = 1/log(2*omega_D/Delta_0)

    g.omega_D = omega_D
    g.w_delta = Delta_0
    g.coupling_lambda = lambda_0

    g.t_t = T

    return g

def nsn_geometry_triple(T, L_s1, L_s2, L_s3, omega_D):
    """
    Split the S wire into three parts
    """
    g = u.Geometry(nnode=4, nwire=3)

    g.t_type = [u.NODE_CLEAN_N_TERMINAL,
                u.NODE_CLEAN_NODE,
                u.NODE_CLEAN_NODE,
                u.NODE_CLEAN_N_TERMINAL]
    g.w_type = [u.WIRE_TYPE_S, u.WIRE_TYPE_S, u.WIRE_TYPE_S]

    g.w_conductance = 1.0

    g.w_ends[0,:] = 0, 1
    g.w_ends[1,:] = 2, 1
    g.w_ends[2,:] = 2, 3
    g.w_length = L_s1, L_s2, L_s3

    Delta_0 = 1.0
    omega_D = omega_D
    lambda_0 = 1/log(2*omega_D/Delta_0)

    g.omega_D = omega_D
    g.w_delta = Delta_0
    g.coupling_lambda = lambda_0

    g.t_t = T

    return g

#------------------------------------------------------------------------------
# Analytical Tc of a NSN structure
#------------------------------------------------------------------------------

def nsn_tc(L_s, omega_D, tol=1e-4):
    r"""
    Compute the Tc of a NSN structure (where N are N-terminals).

    """
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
        x = (a + x)/3.
        v = memo_func(x)
        if v*v0 < 0:
            break
    else:
        return nan

    return optimize.brenth(memo_func, x, last_x, xtol=xtol)

if __name__ == "__main__":
    run_tests()
