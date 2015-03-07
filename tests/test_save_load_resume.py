"""
Test saving, loading and resuming CurrentSolvers.
"""
from testutils import *
from scipy import *
import usadel1 as u
import geometries
import tables

@in_tempdir
def test_currentsolver():
    data = tables.openFile('test.h5', 'w')
    data2 = tables.openFile('test2.h5', 'w')

    # Save

    g = geometries.geometry_NnN(0, 1e-9, 0)
    solver = u.CurrentSolver(g, E=array([1,2,3,4]))
    assert solver.geometry is g
    g.w_length = 2
    solver.solve_spectral()
    solver.save(data.root)

    g2 = geometries.geometry_SnnNnS(0, 1e-9, 0.5)
    solver2 = u.CurrentSolver(g2, E=array([1,2,3,4]))
    assert solver2.geometry is g2
    solver2.solve_spectral()
    solver2.save(data2.root, 'foo/bar/quux')

    # Load

    solver3 = u.CurrentSolver.load(data)
    assert solver3.geometry is not g
    assert alltrue(solver3.geometry.equal(solver.geometry))
    assert alltrue(solver3.coefficient.TT == solver.coefficient.TT)
    solver3.solve_kinetic()

    solver3 = u.CurrentSolver.load(data2.root.foo.bar.quux)
    assert solver3.geometry is not g2
    assert alltrue(solver3.geometry.equal(solver2.geometry))
    assert alltrue(solver3.coefficient.TT == solver2.coefficient.TT)
    solver3.solve_kinetic()

    solver3 = u.CurrentSolver.load(data2.root, 'foo/bar/quux')
    assert solver3.geometry is not g2
    assert alltrue(solver3.geometry.equal(solver2.geometry))
    assert alltrue(solver3.coefficient.TT == solver2.coefficient.TT)
    solver3.solve_kinetic()

    solver3 = u.CurrentSolver.resume(data2.root, g2, 'foo/bar/quux',
                                     E=array([1,2,3,4]))
    assert solver3.geometry is g2
    assert alltrue(solver3.geometry.equal(solver2.geometry))
    assert alltrue(solver3.coefficient.TT == solver2.coefficient.TT)
    solver3.solve_kinetic()

    data.close()
    data2.close()

    # Try loading invalid files

    data = tables.openFile('test.h5', 'w')
    data.close()

    assert raises(lambda: u.CurrentSolver.load(data), IOError)

if __name__ == "__main__":
    run_tests()
