
"""
Test calculated derivatives in a three-probe structure

This test explicitly checks that the various 'd*' variables stored in
the data files match numerically calculated derivatives of variables.
"""
from __future__ import division
from testutils import *
import scipy.interpolate
import scipy.integrate
import usadel1

from geometries import *
import tables
import usadel1.hdf5pickle as hdf5pickle

from usadel1 import *
from usadel1.data import *
from numpy import *

__revision__ = "$Id$"

#
# Initialize:
#
def setup_module():
     global solver, geometry, E, x
     E = linspace_weighed(0, 10, 200, ((.10, 1, 30),))
     x = linspace(0, 1, 101)

     geometry = geometry_SnnNnNnnS(1e-6, 1e-6, 1e-6, 1)
     geometry.t_delta[2:4] = 50
     geometry.w_length = [ 1, 2.1, 2.2, 0.7, 2.3 ]
     geometry.w_length /= sum(geometry.w_length[2:5])
     geometry.w_conductance = [ 4.4, 2.3, 1.2, 4.2, 0.3 ]
     geometry.w_delta = 0

     solver = CurrentSolver(geometry, ne=300, maxE=100)
     solver.solve_spectral_if_needed(calculate_G=False)

def teardown_module():
     global solver, geometry
     solver = None
     geometry = None

def maxnorm(a, axis=None):
     return abs(a).max(axis=axis)

def _derivative_ok(f, path, varname, dername, rtol=5e-2, atol=1e-2):
     """
     Check that the derivative integrates to something close to the
     actual variable.

     This is less sensitive to numerical noise than differentiating
     the variable and comparing to the derivative.

     """
     g = hdf5pickle.load(f, '/geometry')
     var = f.getNode(path + '/' + varname).read()
     der = f.getNode(path + '/' + dername).read()
     x = f.getNode(path + '/x').read()

     dx = (.5*(r_[0, diff(x)] + r_[diff(x), 0])[None,None,:]
           * g.w_length[None,:,None])
     y = var[...,0][...,None] + cumsum(der * dx, axis=-1)

     d = maxnorm(var - y, axis=2) / (atol/rtol + maxnorm(var, axis=2))
     assert is_zero(d, tolerance=rtol)

@in_tempdir
def test_all():
     solver.save('test.h5')
     f = tables.openFile('test.h5')
     try:
         _derivative_ok(f, '/coefficient',    'DL',    'dDL')
         _derivative_ok(f, '/coefficient',    'DT',    'dDT')
         _derivative_ok(f, '/coefficient',    'TT',    'dTT')
         _derivative_ok(f, '/coefficient',   'rjE',   'drjE', atol=1e-2)
         _derivative_ok(f, '/coefficient',   'ijE',   'dijE', atol=1e-2)
         _derivative_ok(f,    '/spectral',     'a',     'da')
         _derivative_ok(f,    '/spectral',     'b',     'db')
     finally:
          f.close()

if __name__ == "__main__":
    run_tests()
