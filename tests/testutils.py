from __future__ import division
import shutil
import tempfile
import subprocess
import tables
import scipy
import scipy.io
import scipy.interpolate
import scipy.integrate
import numpy as _n
import sys, os, errno
import pdb, distutils.util

__revision__ = "$Id$"

from numpy.testing import assert_allclose

class KnownFailure(Exception):
    pass

# Prepare sys.path
def prepare_sys_path():
    dirs = [d for d in os.listdir(os.path.abspath(os.path.join(__file__,
                                                               os.path.pardir,
                                                               os.path.pardir,
                                                               'build')))
            if d == 'lib.%s-%d.%d' % (distutils.util.get_platform(),
                                      sys.version_info[0],
                                      sys.version_info[1])]
    for d in dirs:
        sys.path.insert(0, os.path.join(os.path.pardir, 'build', d))

def run_tests():
    subprocess.call([sys.executable,
                     os.path.join(os.path.dirname(sys.argv[0]), 'test.py'),
                     os.path.basename(sys.argv[0])] + sys.argv[1:])

def curve_max_difference(x1, y1, x2, y2, axis=-1):
    yy2 = scipy.interpolate.interp1d(x2, y2, axis=axis)
    diff = y1 - yy2(x1)
    return max(abs(diff))

def curve_norm2_difference(x1, y1, x2, y2, axis=-1):
    yy2 = scipy.interpolate.interp1d(x2, y2, axis=axis)
    diff = scipy.integrate.trapz((y1 - yy2(x1))**2, x1, axis=axis)
    return _n.sqrt(diff)

def curve_max_relative_difference(x1, y1, x2, y2, axis=-1):
    diff = curve_max_difference(x1, y1, x2, y2, axis=axis)
    norm = max(max(abs(y1)), max(abs(y2)))
    return diff / norm

def curve_norm2_relative_difference(x1, y1, x2, y2, axis=-1):
    diff = curve_max_difference(x1, y1, x2, y2, axis=axis)
    norm = sqrt(scipy.integrate.trapz(y1**2, x1, axis=axis)
                + scipy.integrate.trapz(y2**2, x2, axis=axis))
    return diff / norm

def array_abs_round(x, digits):
    m = array_normmax(x)
    if m == 0: return x
    ndigits = _n.floor(_n.log(m)/_n.log(10))
    divisor = 10**(ndigits - digits + 1)
    x = _n.around(x / divisor)*divisor
    return x

def array_norm2(x):
    return _n.sqrt(sum(abs(_n.ravel(x))**2))

def array_normmax(x):
    return max(abs(_n.ravel(x)))

def arrays_equal(x, y, tolerance=None, relative=None):
    """Test whether two arrays are equal. Return ``True`` if yes, \
    ``False`` if no."""
    # numarray's alltrue returns 1 instead of True: work around this
    if tolerance is None and relative is None:
        tolerance = 0
    diff = array_normmax(x - y)
    if tolerance is not None and diff > tolerance:
        return False
    if relative is not None and diff > relative * max(array_normmax(x),
                                                      array_normmax(y)):
        return False
    return True

def load_ascii_array(fn):
    return _n.loadtxt(fn, comments='%')

def is_zero(value, tolerance):
    if _n.all(abs(value) < tolerance):
        return True
    else:
        raise ValueError(value, 'is not zero within tolerance %g'%tolerance)

def raises(obj, excclasses, *args, **kwargs):
    if not hasattr(excclasses, '__len__'):
        excclasses = (excclasses,)

    try:
        obj(*args, **kwargs)
    except tuple(excclasses):
        return True
    else:
        return False

def rm_f(filename):
    try:
        os.unlink(filename)
    except OSError, e:
        if e.errno != errno.ENOENT:
            raise

def in_tempdir(func):
    def wrapper(*a, **kw):
        cwd = os.getcwd()
        tmpdir = os.path.abspath(tempfile.mkdtemp(prefix='usadel1-test-'))
        try:
            os.chdir(tmpdir)
            return func(*a, **kw)
        finally:
            shutil.rmtree(tmpdir)
            os.chdir(cwd)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper.__module__ = func.__module__
    return wrapper

def slow(func):
    func.__dict__.setdefault('__tags__', []).append('slow')
    return func

_old_tables_openfile = None
if _old_tables_openfile is None:
    _old_tables_openfile = tables.openFile
    def _openFile(fn, *a, **kw):
        # Work around hdf5 file name caching in temporary dirs?
        fn = os.path.abspath(fn)
        return _old_tables_openfile(fn, *a, **kw)
    tables.openFile = _old_tables_openfile
