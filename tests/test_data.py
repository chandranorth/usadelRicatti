"""
Test HDF5 datafile utilities

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *
from numpy import *
import tables, scipy, sys, os, random
from usadel1 import data

__revision__ = "$Id: test_data.py 3184 2006-10-02 07:05:45Z pauli $"

#
# First, introduce simple serializable classes:
#
class SimpleSer(object):
     def __init__(self):
         self.a = scipy.rand(10, 20, 30)
         self.b = scipy.rand()
         self.c = scipy.rand()
         self.x = float(random.random())
#
class SimpleESer(object):
     def __init__(self):
         self.a = scipy.rand(2, 3, 4)
         self.b = scipy.rand()
         self.c = scipy.rand()
         self.d = scipy.rand(2, 4, 3)
#
# Serialization testing function: save and load data from file
#

def do_serialize(d, path='obj'):
     hdf5 = tables.openFile('test.h5', 'w')
     data.save(hdf5.root, path, d)
     hdf5.close()
     d = d.__class__()
     hdf5 = tables.openFile('test.h5', 'r+')
     d = data.load(hdf5.root, path)
     data.save(hdf5.root, path, d)
     d = d.__class__()
     d = data.load(hdf5.root, path)
     return (hdf5, d)
#
#
# Simple test
# -----------
#
# First, a simple test of serialization
#

@in_tempdir
def test_1():
    d = SimpleSer()
    hdf5, d2 = do_serialize(d)
    assert arrays_equal(d.a, d2.a)
    assert arrays_equal(d.b, d2.b)
    assert arrays_equal(d.c, d2.c)
    hdf5.close()
#
# More complicated test
# ---------------------
#
# A more complicated test
#

@in_tempdir
def test_2():
    d = SimpleSer()
#
    f = tables.openFile('test.h5', 'w')
    data.save(f.root, 'fubar', d)
    d2 = data.load(f.root, 'fubar')
    assert arrays_equal(d.a, d2.a)
    assert arrays_equal(d.b, d2.b)
    assert arrays_equal(d.c, d2.c)
    f.close()
#
# Loading partial data:
#
    f = tables.openFile('test.h5', 'r')
    da = data.load(f.root, 'fubar/a')
    dx = data.load(f.root, 'fubar/x')
    f.close()
    assert arrays_equal(d.a, da)
    assert dx == d.x
    assert type(dx) == float
#
#
# Another test
# ------------
#

@in_tempdir
def test_3():
    d = SimpleESer()
#
# Check that we can save the whole small piece
#
    hdf5, d2 = do_serialize(d)
    assert arrays_equal(d.a, d2.a)
    assert arrays_equal(d.b, d2.b)
    assert arrays_equal(d.c, d2.c)
    assert arrays_equal(d.d, d2.d)
    hdf5.close()
#
#
# More testing: saving into subgroups
# -----------------------------------
#

@in_tempdir
def test_save_subgroups():
    a = scipy.rand(1,2,3,4,5)

    hdf5 = tables.openFile('test.h5', 'w')
    data.save(hdf5.root, 'foo/bar/quux', a)
    hdf5.close()

    hdf5 = tables.openFile('test.h5', 'r')
    b = data.load(hdf5.root, 'foo/bar/quux')
    hdf5.close()

    assert arrays_equal(a, b)

#
# Type mapping test
# -----------------
#

@in_tempdir
def test_type_mapping():
    x = [1,2,3,4,5]
    data._type_map[int] = float32
    hdf5 = tables.openFile('test.h5', 'w')
    data.save(hdf5.root, 'fubar', x)
    assert hdf5.root.fubar.atom.dtype == float32
    hdf5.close()
    del data._type_map[int]
#
# 64-bitness test
# ---------------
#
# We do not want to save 64 bit ints, these cause problems for Matlab
# 14R3 (nevertheless, R2006a seems to be fixed).
#

@in_tempdir
def test_64_bitness():
    x = [1,2,3,4,5]
    y = array([1,2,3,4,5,6], dtype=int64)
    data._type_map[int] = int32
    hdf5 = tables.openFile('test.h5', 'w')
    data.save(hdf5.root, 'fubar', x)
    data.save(hdf5.root, 'baz', y)
    assert hdf5.root.fubar.atom.dtype == int32
    assert hdf5.root.baz.atom.dtype == int32
    hdf5.close()

if __name__ == "__main__":
    run_tests()
