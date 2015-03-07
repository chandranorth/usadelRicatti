"""
Test utility functions from `usadel1.util`

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *
from numpy import *
import usadel1.util as _u
import scipy

__revision__ = "$Id: test_util.py 3184 2006-10-02 07:05:45Z pauli $"

# 
# Simplest things first:
# 
def test_linspace():
    assert arrays_equal(_u.linspace(0, 1, 29), scipy.linspace(0, 1, 29))

# 
# Weighed spacing:
# 
def test_linspace_weighed():
    x_known = array([0.,  0.23390597,  0.50844397,  0.77746496,  1.])
    x = _u.linspace_weighed(0, 1, 5, ((0, 1, 1), (1, 1, 1)))
    assert arrays_equal(x, x_known, tolerance=1e-8)

# 
# Norm2
# 
def test_norm2():
    x = array([[1, 2], [5, 9.45], [3, 1.32], [6.44+2j, -3+2j]])
    assert is_zero(_u.norm2(x) - 13.730203931478949, tolerance=1e-9)

# 
# "Stretching" arrays:
# 
def test_stretch_shape():
    x = array([1, 2, 3])
    y = array([[1,2,3], [1,2,3], [1,2,3], [1,2,3]])
    z = _u.stretch_shape(x, (4, -1))
    assert arrays_equal(x, z)
 
# 
# Polar coordinates
# -----------------
# 
# Complex polar coordinate transformations. First scalars,
# 
def test_polar():
    z = 1 + 1j
    mag, arg = _u.complex_to_polar(z)
    assert is_zero(mag - sqrt(2), tolerance=1e-10)
    assert is_zero(arg - pi/4, tolerance=1e-10)
# 
# Then arrays
# 
    z = array([[1, 2], [3, 4]]) + array([[4, 3], [2, 1]])*1j
    mag, arg = _u.complex_to_polar(z)
    assert arrays_equal(mag - abs(z), 0, tolerance=1e-10)
    assert arrays_equal(arg - log(z).imag, 0, tolerance=1e-10)
# 
# Then do the opposite
# 
    assert is_zero(1+1j - _u.polar_to_complex(sqrt(2),pi/4),
                    tolerance=1e-10)
# 
# Make phase "continuous":
# 
    x = _u.linspace(0, 10, 20) * pi
    y = array(x, copy=True)
    y[3] -= 1*2*pi
    y[9] -= 6*2*pi
    y[12] += 1*2*pi
    z = _u.continuous_phase(y)
    assert arrays_equal(x, z, tolerance=1e-10)
# 
# Also, "center" the phase:
# 
    z = _u.continuous_phase(y, center=True)
    assert max(abs(z)) <= 6.000001*pi
    assert arrays_equal(exp(1j*z), exp(1j*x), tolerance=1e-10)
    assert arrays_equal(exp(1j*y), exp(1j*x), tolerance=1e-10)
# 
# Then the same, but along a different axis
# 
    x = _u.stretch_shape(_u.linspace(0, 10, 20) * pi, (5, 6, -1))
    y = array(x, copy=True)
    y[:,:,3] -= 2*2*pi
    y[:,:,2] -= 1*2*pi
    y[:,:,5] += 3*2*pi
    z = _u.continuous_phase(y, axis=2)
    assert arrays_equal(x, z, tolerance=1e-10)
 
# 
# Splitting vectors
# -----------------
# 
def test_divide_vector_to_chunks():
# 
    x = array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    y = _u.divide_vector_to_chunks(x, [0.5, 4])
    z = [array([0]), array([1, 2, 3]), array([4, 5, 6, 7, 8, 9])]
    for a, b in zip(y, z):
         assert arrays_equal(a, b)
 
# 
# Assign-by-value array members
# -----------------------------
#
def test_ArrayProperty():
    class Cls(object):
        x = _u.ArrayProperty('x')
        y = _u.ArrayProperty('y', settable=True)

        def __init__(self):
            self.x = array([1,2,3,4,5])
            self.__dict__['y'] = array([1,2,3,4,5])
    a = Cls()
    a.x = 0
    a.y = 0
    assert arrays_equal(a.x, array([0, 0, 0, 0, 0]))
    assert arrays_equal(a.y, array([0, 0, 0, 0, 0]))

# 
# Sliceable instances
# -------------------
# 
def test_Sliceable():
    class Cls(_u.Sliceable):
         def __init__(self):
             self.x = array([ [1,2,3],
                              [4,5,6] ])
             self.y = array([7,8,9])
             _u.Sliceable.__init__(self, {'x': (0, 1), 'y': (0, None)})
    a = Cls()
    b = a[1:2, :]

    assert a.x.shape == (2, 3)
    assert b.x.shape == (1, 3)
    assert a.y.shape == (3,)
    assert b.y.shape == (1,)
    assert arrays_equal(a.x[1:2,:], b.x)
    assert arrays_equal(a.y[1:2], b.y)

    assert raises(lambda: a[1,2,3], KeyError)
    assert raises(lambda: a[:], KeyError)

if __name__ == "__main__":
    run_tests()
