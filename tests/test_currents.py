"""
Some (outdated) tests that use `usadel1.currents`

They function, but the style does not conform well to the current usadel1 code.
"""
from __future__ import division
from testutils import *
from scipy import *
from geometries import *
from usadel1 import *

__revision__ = "$Id: test_currents.py 3279 2007-02-15 14:55:41Z pauli $"

def setup():
    global T0, V0, phi, x, E

    T0 = 1e-6
    V0 = 15
    phi = 1

    x = linspace_weighed(0, 1, 71, ( (0, 1, 0.1), (1, 1, 0.1) ))
    E = linspace_weighed(0, 250, 500, ((0, 1, 20),))

def test_try_SnnNnS_IV():
    """Evaluate the currents for the three-probe setup."""
    global T0, V0, phi, x, E

    geometry = geometry_SnnNnS(V0, T0, phi)

    sc = CurrentSolver(geometry, E, chunksize=300)
    sc.solve_spectral_if_needed(calculate_G=False)

    Vs = linspace(0, 20, 3)
    for V in Vs:
        sc.geometry.t_mu[2] = V
        print "Solving for V = ", V
        Ic, Ie = sc.get_currents(ix=1)
    # FIXME: no asserts

@slow
def test_currents_lazy():
    """Evaluate the currents for the three-probe setup lazily."""
    global T0, V0, phi, x, E

    geometry = geometry_SnnNnS(V0, T0, phi)
    solver = CurrentSolver(geometry, E=E)
    solver.set_solvers(sp_solver=SP_SOLVER_TWPBVP,
                       kin_solver=KIN_SOLVER_TWPBVP)
    solver.solve()
    solver.set_solvers(kin_solver=KIN_SOLVER_BLOCK)
    solver.calculate_G()
    solver.set_solvers(kin_solver=KIN_SOLVER_TWPBVP)

    Vs = [1e-5, 1e-3, 1e-1, 1, 1e1]
    Ts = [1e-5, 1e-3, 1e-1, 1, 1e1]

    def equilibrium_terminal_f(E, mu, T):
        fL = .5*(_n.tanh(.5*(E + mu)/T) + _n.tanh(.5*(E - mu)/T))
        fT = .5*(_n.tanh(.5*(E + mu)/T) - _n.tanh(.5*(E - mu)/T))
        return fL, fT

    for V in Vs:
        for T in Ts:
            geometry.t_T = T
            geometry.t_mu[2] = V

            Ic, IE = solver.get_currents_lazy(epsabs=1e-5, epsrel=1e-3, ix=0)
            Ic2, IE2 = solver.get_currents(ix=0)
            Ic3, IE3 = solver.get_currents_from_G(ix=0)

            assert allclose(Ic, Ic2, atol=1e-2, rtol=1e-2), locals()
            assert allclose(IE, IE2, atol=1e-2, rtol=1e-2), locals()
            assert allclose(Ic, Ic3, atol=1e-2, rtol=1e-2), locals()
            assert allclose(IE, IE3, atol=1e-2, rtol=1e-2), locals()

def test_SnnNnNnnS():
    """Non-trivial test for a four-probe setup: find the thermopower."""
    global T0, V0, phi, x, E

    T = [2.3360300000, 2.3832200000, 2.3596300000]

    geometry = geometry_SnnNnNnnS(T[0], T[1], T[2], phi)
    sc = CurrentSolver(geometry, E, chunksize=100)
    #sc.solver.set_solvers(kin_solver=3, sp_solver=2)
    sc.solve_spectral_if_needed(calculate_G=False)

    sc.geometry.t_t[0] = T[0]
    sc.geometry.t_t[1] = T[1]

    print sc.geometry.get_idstr()

    sc.solver.set_geometry(sc.geometry)
    sc.solver.set_kinetic(sc.coefficient)

    def zero_current():
        Ic, Ie = sc.get_currents_lazy(w_jT=[0,1], w_jL=[])
        return [Ic[0], Ic[1]]

    def adjust_potentials(z):
        print "<----", z
        sc.geometry.t_mu[0] = z[0]
        sc.geometry.t_mu[1] = z[1]

    z0 = array([0, 0], float64)
    z = optimize_parameters_for(z0, zero_current, adjust_potentials,
                                xtol=1e-3, epsfcn=1e-5)

    mu = array(z)
    print "====> MU = ", mu

    ## Compare to an old result
    target_mu = array([-0.22640172484e-2, -0.22631056449e-2])
    assert abs((mu - target_mu)/mu).max() < 1e-1, \
                    "The error in resulting voltages is too large!"

@slow
def test_try_SnnNnNnnS_thompson():
    """Non-trivial test for a four-probe setup: find Thompson effect."""
    global T0, V0, phi, x, E

    #T = [2, 2, 1]
    T = [4.136743, 2.000000, 3.000000]

    geometry = geometry_SnnNnNnnS(T[0], T[1], T[2], -phi)
    sc = CurrentSolver(geometry, E, chunksize=100)
    sc.solve_spectral_if_needed(calculate_G=False)

    def zero_thermal_current():
        Ic, Ie = sc.get_currents_lazy(w_jT=[0], w_jL=[0])
        thcur = Ic[0] * sc.geometry.t_mu[0] + Ie[0]
        print "---->", Ic[0], Ie[0], thcur
        return [ thcur ]

    def adjust_temperature(z):
        print "<----", z
        sc.geometry.t_t[0] = abs(z[0])

    sc.geometry.t_t[0] = T[0]
    sc.geometry.t_t[1] = T[1]

    sc.geometry.t_mu[0] = -4
    sc.geometry.t_mu[1] = -1

    z0 = array([1], float64)
    z = optimize_parameters_for(z0, zero_thermal_current,
                                adjust_temperature,
                                xtol=1e-3, epsfcn=1e-5)

    tt = array(z)
    print "====> T = ", tt
    assert (tt[0] - T[0])/T[0] < 1e-4

@slow
def test_SnnNnNnnS_thompson_compare_old():
    """Non-trivial test for a four-probe setup: find Thompson effect,
        compare to old code."""
    global T0, V0, phi, x, E

    #T = [2, 2, 1]
    T = [4.547578, 4.503250, 3.000000]

    geometry = geometry_SnnNnNnnS(T[0], T[1], T[2], -phi)
    geometry.w_length[0] = 0.5 / (1 + 1 + .75)
    geometry.w_length[1] = 0.5  / (1 + 1 + .75)
    geometry.w_length[2] = 1.0  / (1 + 1 + .75)
    geometry.w_length[3] = 1.0  / (1 + 1 + .75)
    geometry.w_length[4] = 0.75 / (1 + 1 + .75)
    sc = CurrentSolver(geometry, E, chunksize=100)
    sc.solve_spectral_if_needed(calculate_G=False)

    def zero_thermal_current():
        Ic, Ie = sc.get_currents_lazy(w_jT=[0], w_jL=[0])
        thcur = Ic[0] * sc.geometry.t_mu[0] + Ie[0]
        print "---->", Ic[0], Ie[0], thcur
        return [ thcur ]

    def adjust_temperature(z):
        print "<----", z
        sc.geometry.t_t[0] = abs(z[0])

    sc.geometry.t_t[0] = T[0]
    sc.geometry.t_t[1] = T[1]

    sc.geometry.t_mu[0] = 0
    sc.geometry.t_mu[1] = 1

    z0 = array([1], float64)
    z = optimize_parameters_for(z0, zero_thermal_current,
                                adjust_temperature,
                                xtol=1e-3, epsfcn=1e-5)

    tt = array(z)
    print "====> T = ", tt
    assert (tt[0] - T[0])/T[0] < 1e-4

def test_try_SnnNnNnnS_IV():
    """Evaluate the currents for the four-probe setup."""
    global T0, V0, phi, x, E

    T = [2.9780900000, 2.9810700000, 1]
    geometry = geometry_SnnNnNnnS(T[0], T[1], T[2], phi)

    sc = CurrentSolver(geometry, E, chunksize=300)
    sc.solve_spectral_if_needed(calculate_G=False)

    dTs = linspace(-1e-2, 1e-2, 5)
    E = linspace(0, 30, 200)
    print "Interpolating coefficients..."
    for dT in dTs:
        print "Solving for dT = ", dT

        geometry.t_mu[0] = dT
        geometry.t_mu[1] = dT

        sc.solve_kinetic()

        Ic, Ie = sc.get_currents(ix=1)
        print "Ic = ", Ic[:]

    # FIXME: no asserts


def test_SnnNnNnnS_IV_compare_to_old_code():
    """Evaluate the currents for the four-probe setup, comparing \
    to the old code.
    """
    global T0, V0, phi, x, E

    # Setup geometry
    T = [0.4, 0.2, 1]
    geometry = geometry_SnnNnNnnS(T[0], T[1], T[2], -phi)
    geometry.t_mu[0] = 1e-3
    geometry.t_mu[1] = 4e-3

    # Evaluate spectral coefficients
    sc = CurrentSolver(geometry, E, chunksize=300)
    sc.solver.set_solvers(kin_solver=3, sp_solver=2)
    sc.solve_spectral_if_needed(calculate_G=False)
    sc.solve_kinetic()

    # Evaluate currents
    E = sc.kinetic.E
    jL, jT = sc.kinetic.jL[:,:,1], sc.kinetic.jT[:,:,1]

    # Compare solutions
    c = load_ascii_array('test_currents.1.dat')
    cE = c[:,0]
    cjT = c[:,1:6] * 3
    cjL = c[:,6:11] * 3

    for w in range(0,5):
        difference = curve_max_difference(cE, cjT[:,w], E, jT[:,w])
        rel_difference = curve_max_relative_difference(cE, cjT[:,w],
                                                       E, jT[:,w])
        print "max-norm (abs., rel.) difference in jT_%d:" % w, \
              difference, rel_difference
        assert difference < 2e-2
        assert rel_difference < 2e-2

        difference = curve_max_difference(cE, cjL[:,w], E, jL[:,w])
        rel_difference = curve_max_relative_difference(cE, cjL[:,w],
                                                       E, jL[:,w])
        print "max-norm (abs., rel.) difference in jL_%d:" % w, \
              difference, rel_difference
        assert difference < 2e-2
        if w == 2 or w == 3:
            # No point in comparing these, they are zero
            pass
        else:
            assert rel_difference < 2e-2

def test_SnnNnNnnS_IV_compare_to_old_code_2():
    """Evaluate the currents for the four-probe setup, comparing \
    to the old code.
    """
    global T0, V0, phi, x, E

    # Setup geometry
    T = [0.4, 0.2, 1]
    geometry = geometry_SnnNnNnnS(T[0], T[1], T[2], -phi)
    geometry.t_mu[0] = 1e-3
    geometry.t_mu[1] = 4e-3

    # Evaluate spectral coefficients
    sc = CurrentSolver(geometry, E)
    sc.solver.set_solvers(kin_solver=3, sp_solver=2)
    sc.solve_spectral_if_needed(calculate_G=False)
    sc.solve_kinetic()

    # Evaluate currents
    E = sc.kinetic.E
    jL, jT = sc.kinetic.jL[:,:,1], sc.kinetic.jT[:,:,1]

    # Compare solutions
    c = load_ascii_array('test_currents.1.dat')
    cE = c[:,0]
    cjT = c[:,1:6] * 3
    cjL = c[:,6:11] * 3

    for w in range(0,5):
        difference = curve_max_difference(cE, cjT[:,w], E, jT[:,w])
        rel_difference = curve_max_relative_difference(cE, cjT[:,w],
                                                       E, jT[:,w])
        print "max-norm (abs., rel.) difference in jT_%d:" % w, \
              difference, rel_difference
        assert difference < 2e-2
        assert rel_difference < 2e-2

        difference = curve_max_difference(cE, cjL[:,w], E, jL[:,w])
        rel_difference = curve_max_relative_difference(cE, cjL[:,w],
                                                       E, jL[:,w])
        print "max-norm (abs., rel.) difference in jL_%d:" % w, \
              difference, rel_difference
        assert difference < 2e-2
        if w == 2 or w == 3:
            # No point in comparing these, they are zero
            pass
        else:
            assert rel_difference < 2e-2

if __name__ == "__main__":
    run_tests()
