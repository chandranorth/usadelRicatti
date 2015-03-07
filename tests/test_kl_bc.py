"""
Test the Kupriyanov-Lukichev boundary condition

"""
from __future__ import division
from testutils import *
from scipy import *
import usadel1 as u

def test_kl_spectral_almost_clean():
    """
    Compare spectral solutions for a low-resistance tunnel junction
    and a clean node between two S terminals
    """

    g1 = u.Geometry(2, 3)
    g1.t_type = [u.NODE_CLEAN_S_TERMINAL,
                 u.NODE_TUNNEL_NODE,
                 u.NODE_CLEAN_S_TERMINAL]

    g2 = u.Geometry(2, 3)
    g2.t_type = [u.NODE_CLEAN_S_TERMINAL,
                 u.NODE_CLEAN_NODE,
                 u.NODE_CLEAN_S_TERMINAL]

    Delta_0 = 200
    r = 1e-3
    
    for g in (g1, g2):
        g.t_delta = Delta_0, 0, Delta_0
        g.t_phase = pi/4, 0, -pi/4
        g.t_inelastic = 1e-9
        g.w_type = u.WIRE_TYPE_N
        g.w_length = 0.5, 0.5
        g.w_conductance = 1
        g.w_ends[0,:] = 0, 1
        g.w_ends[1,:] = 2, 1
        
    g1.t_resistance = 0, r, 0

    results = []
    for g in (g1, g2):
        solver = u.CurrentSolver(g, ne=500)
        solver.solve_spectral()

        results.append(solver.spectral)

    # Check that the results are close
    assert_allclose(results[0].a, results[1].a,
                    atol=2*r, rtol=2*r)
    assert_allclose(results[0].b, results[1].b,
                    atol=2*r, rtol=2*r)

def test_kl_sc_snins():
    """
    Compute supercurrent in a setup where high-temperature tunnel-junction
    asymptotic expression is known:  S----T----S, i.e., SNINS

    Tests only the spectral equations.

    """

    g = u.Geometry(2, 3)

    g.t_type = [u.NODE_CLEAN_S_TERMINAL,
                u.NODE_TUNNEL_NODE,
                u.NODE_CLEAN_S_TERMINAL]

    Delta_0 = 2000

    g.t_delta = Delta_0, 0, Delta_0
    g.t_phase = pi/4, 0, -pi/4
    g.t_inelastic = 1e-9
    g.w_type = u.WIRE_TYPE_N
    g.w_length = 0.5, 0.5
    g.w_conductance = 1
    g.w_ends[0,:] = 0, 1
    g.w_ends[1,:] = 2, 1

    for r in [1.0, 5.0, 15.0]:
        g.t_resistance = 0, r, 0

        Ts = linspace(5, 50, 150)
        I = []

        # Analytical prediction
        k = sqrt(2*pi*Ts)
        I_asy = 4*k**3 * (
            (4*tan(pi/8)*exp(-0.5*k))**2 * 1/(r*k)
            / (1 + 1/(r*k))**2
            )

        # Slightly wrong analytical expression
        I_asy_wrong = 4*k**3 * (
            (4*tan(pi/8)*exp(-0.5*k))**2 * 1/(r*k)
            )

        # Numerics
        solver = u.CurrentSolver(g, ne=500)
        solver.solve_spectral()
        for T in Ts:
            g.t_t = T
            g.t_mu = 0

            Ic, IE = solver.get_currents(ix=0)
            I.append(-Ic[0])


        # Compare
        assert_allclose(I_asy, I, rtol=0.04, atol=1e-3*I_asy[0],
                        err_msg="R = %r" % r)

        if allclose(I_asy_wrong, I, rtol=0.04, atol=1e-3*I_asy[0]):
            raise AssertionError("Current matches to wrong value")


def test_kl_sc_sins():
    """
    Compute supercurrent in a setup where high-temperature tunnel-junction
    asymptotic expression is known: ST-----S, i.e., SINS

    Tests only the spectral equations.

    """
    Delta_0 = 2000

    # Construct SIN boundary condition from S-T-----S
    g = u.Geometry(2, 3)
    g.t_type = [u.NODE_CLEAN_S_TERMINAL,
                u.NODE_TUNNEL_NODE,
                u.NODE_CLEAN_S_TERMINAL]
    g.t_delta = Delta_0, 0, Delta_0
    g.t_phase = pi/4, 0, -pi/4
    g.t_inelastic = 1e-9
    g.w_type = u.WIRE_TYPE_N
    g.w_length = 0.01, 1.0
    g.w_conductance = 1
    g.w_ends[0,:] = 0, 1
    g.w_ends[1,:] = 2, 1

    # The SIN boundary condition is also implemented separately,
    # so we can also compare against that
    g2 = u.Geometry(1, 2)
    g2.t_type = [u.NODE_TUNNEL_S_TERMINAL,
                 u.NODE_CLEAN_S_TERMINAL]
    g2.t_delta = Delta_0, Delta_0
    g2.t_phase = pi/4, -pi/4
    g2.t_inelastic = 1e-9
    g2.w_type = u.WIRE_TYPE_N
    g2.w_length = 1
    g2.w_conductance = 1
    g2.w_ends[0,:] = 0, 1

    # Check the results for a few different resistance values
    for r in [1.0, 5.0, 15.0]:
        g.t_resistance = 0, r, 0
        g2.t_resistance = r, 0

        Ts = linspace(5, 50, 150)

        # Analytical prediction
        k = sqrt(2*pi*Ts)
        I_asy = 4*k**3 * (
            4*tan(pi/8)*exp(-1.0*k) * 1/(2*r*k)
            )

        # Numerics
        I = []
        I2 = []

        for gx, Ix in [(g, I), (g2, I2)]:
            solver = u.CurrentSolver(gx, ne=500)
            solver.solve_spectral()

            for T in Ts:
                gx.t_t = T
                gx.t_mu = 0

                Ic, IE = solver.get_currents(ix=0)
                Ix.append(-Ic[0])

        # Compare
        assert_allclose(I_asy, I, rtol=0.04, atol=1e-3*I_asy[0],
                        err_msg="R = %r" % r)
        assert_allclose(I_asy, I2, rtol=0.02, atol=1e-3*I_asy[0],
                        err_msg="R = %r" % r)

if __name__ == "__main__":
    run_tests()
