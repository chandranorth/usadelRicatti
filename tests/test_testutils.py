"""
Testing testutils

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *
from numpy import *

__revision__ = "$Id: test_testutils.py 3184 2006-10-02 07:05:45Z pauli $"

# 
# Basics
# ------
# 
# Test some basic routines:
# 
def basic_tests():
    assert is_zero(0, tolerance=1e-99)
    assert is_zero(1, tolerance=1.1)
    assert raises(lambda: is_zero(1, tolerance=0.9), ValueError)
# 
# More of the same:
# 
    x = array([1.5, 2.5, 3.9, 4.5])
    y = array([2.0 + 1j, 2.0 + 3j, 4.0 - 1.1j, 5.2 + 2j])
    assert is_zero(array_norm2(x - y) - 4.026164427839479, tolerance=1e-10)
    assert is_zero(array_norm2(x) - 6.6302337817003103, tolerance=1e-10)
    assert is_zero(array_normmax(x) - 4.5, tolerance=1e-10)
    assert is_zero(array_normmax(y) - 5.5713553108736482, tolerance=1e-10)

    assert not arrays_equal(x, y)
    assert arrays_equal(x, x)

# 
# Curve comparison
# ----------------
# 
# First, define some test data, curves that are close together:
# 
def test_curve_comparison():

    def get_good_test_data():
         data = {}
         data['x1'] = arange(0, 1, 0.01)
         data['x2'] = arange(0, 3, 0.5)
         data['y1'] = sin(cosh(data['x1']))
         data['y2'] = sin(cosh(data['x2']) + 1e-3*data['x2'])
         maxdiff = 0.014546
         norm2diff = 0.00942
         return (data, maxdiff, norm2diff)
# 
# ``maxdiff`` and ``norm2diff`` here are precalculated. Then, specify
# curves that are not close together:
# 
    def get_bad_test_data():
         data = {}
         data['x1'] = arange(0, 1, 0.01)
         data['x2'] = arange(0, 3, 0.5)
         data['y1'] = sin(cosh(data['x1']))
         data['y2'] = cos(cosh(data['x2']) + 1e-3*data['x2'])
         maxdiff = 0.973
         norm2diff = 0.5924
         return (data, maxdiff, norm2diff)
# 
# This data has wrong dimensions for 'y1':
# 
    def get_invalid_test_data():
         data = {}
         data['x1'] = arange(0, 1, 0.01)
         data['x2'] = arange(0, 3, 0.5)
         data['y1'] = 0
         data['y2'] = data['x1']
         return data
# 
# Then, perform some tests: require that the calculated maximum and
# norm2 differences values should agree with precomputed ones:
# 
    good = get_good_test_data()
    bad = get_bad_test_data()
    invalid = get_invalid_test_data()
        
    good_max_diff = curve_max_difference(**good[0])
    bad_max_diff = curve_max_difference(**bad[0])

    good_norm2_diff = curve_norm2_difference(**good[0])
    bad_norm2_diff = curve_norm2_difference(**bad[0])

    assert is_zero((good_max_diff - good[1])/good[1], tolerance=1e-3)
    assert is_zero((bad_max_diff - bad[1])/bad[1], tolerance=1e-2)
    assert is_zero((good_norm2_diff - good[2])/good[2], tolerance=1e-3)
    assert is_zero((bad_norm2_diff - bad[2])/bad[2], tolerance=2e-2)
# 
# Invalid data should cause errors:
#         
    assert raises(lambda: curve_max_difference(**invalid), ValueError)
    assert raises(lambda: curve_norm2_difference(**invalid), ValueError)

#
# Loading arrays from text files
# ------------------------------
# 
def test_array_loading():
    x = load_ascii_array('test_threeprobe.1.dat')
    y = load_ascii_array('test_threeprobe.2.dat')
    x2 = load_ascii_array('test_threeprobe.1.dat')
    assert is_zero(max(abs(ravel(x - x2))), tolerance=1e-10)

# 
# Rounding
# --------
# 
def test_rounding():
    assert array_abs_round(array([12345.6789101112]), digits=2)[0] == 12000
    assert array_abs_round(array([12505.6789101112]), digits=2)[0] == 13000
    x = array([123.4567891011, 12.505678910, 1.234567890, 0.123456789])
    assert arrays_equal(array_abs_round(x, digits=3),
                        [123.0, 13.0, 1.0, 0.0])

if __name__ == "__main__":
    run_tests()
