"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

Utility functions
=================

Common utility functions.

"""

from __future__ import division
import datetime as datetime
import scipy, scipy.integrate, scipy.interpolate
import numpy as _n

__all__ = ['linspace', 'linspace_weighed', 'complex_to_polar',
           'polar_to_complex', 'estimate_time', 'continuous_phase',
           'divide_vector_to_chunks', 'stretch_shape',
           'norm2', 'Sliceable', 'ArrayProperty']
__docformat__ = "restructuredtext en"

def linspace(a,b,n):
    """Returns an array of n points in the range [a,b].

       First point is always at a, and last point at b."""
    if n > 0:
        return a + ((float(b) - float(a))/(n-1)) * _n.arange(n)
    else:
        return _n.array([a])

def linspace_weighed(a, b, n, points):
    """Positions 'n' points in range ['a', 'b'] such that space around
       points[:][0] has an additional weight points[:][1] and
       half-width points[:][2].

       The shape of the weight is ``(|x - x0|/w + 1)**(-2)``, so that if the
       range is infinite, w is indeed the half-width of the distribution.
       points[:][1] describes the relative weights of such peaks."""

    x = linspace(a, b, 5*n)
    density = _n.zeros([len(x)], _n.float_)

    for point in points:
        point = list(point)
        
        if point[0] < a:
            point[0] = a
        if point[0] > b:
            point[0] = b

        density_shape = 1 / (abs(x - point[0])/point[2] + 1)**2
        base_weight = scipy.integrate.trapz(density_shape, x)

        density += (point[1] / base_weight) * density_shape

    if len(points) == 0:
        density[:] = 1

    cumdensity = scipy.cumsum(density) - density[0]
    cumdensity /= cumdensity[-1]
    interpolant = scipy.interpolate.interp1d(cumdensity, x)

    y = linspace(0, 1, n)
    return interpolant(y)

def complex_to_polar(z, continuous_axis=0, center=False):
    """Return the polar representation (abs, phase) of complex numbers.

    The phase can be made to vary continuously along a given axis and
    also centered around zero. See function continuous_phase.
    """
    return (abs(z), continuous_phase(_n.angle(z),
                                     axis=continuous_axis,
                                     center=center))

def polar_to_complex(r, phi):
    """Return the complex number corresponding to its polar representation."""
    return r * _n.exp(1j * phi)

def dump_data(r):
    for ie in range(r.ne):
        for ix in range(r.nx):
            print "%15.7g %15.7g %15.7g %15.7g" % (
                r.E[ie], r.x[ix],
                r.theta.real[ix,ie,0],
                r.theta.imag[ix,ie,0] )

def estimate_time(start, k, kmax):
    """Form a string describing the elapsed time and estimate completion.

    'start' is the initial time as returned by datetime.datetime.now()
    'k' and 'kmax' are positive scalars describing current progress.
    """
    now = datetime.datetime.now()
    diff = now - start
    secs = diff.seconds + diff.days * 24 * 3600

    etasecs = secs*(abs(kmax)+1)/(abs(k)+1)
    
    return "Elapsed: %02d:%02d:%02d. Estimated: %02d:%02d:%02d" % (
        int(secs/60/60), int(secs/60) % 60, int(secs) % 60,
        int(etasecs/60/60), int(etasecs/60) % 60, int(etasecs) % 60,
        )

def continuous_phase(phase, axis=0, center=False):
    """Add and subtract 2 pi such that the phase in the array is
       as continuous as possible, along first or given axis. Optionally,
       it also centers the phase data so that the average is smallest."""

    phase = _n.array(phase, copy=0)

    rowshape = list(phase.shape)
    
    if len(rowshape) > 0:
        rowshape[axis] = 1

        slip = _n.concatenate([ _n.zeros(rowshape),
                                scipy.diff(phase, axis=axis) ],
                              axis=axis)
        slip = _n.around(slip/(2*_n.pi))
        cumslip = scipy.cumsum(slip, axis=axis)

        phase = phase - 2*_n.pi*cumslip
    else:
        pass

    if center:
        offset = _n.around(scipy.average(phase, axis=axis)/(2*_n.pi))
        offset = _n.reshape(offset, rowshape)
        offset = _n.repeat(offset, cumslip.shape[axis], axis=axis)
        phase = phase - 2*_n.pi*offset
    
    return phase

def divide_vector_to_chunks(v, p):
    """Sort a vector and divides split it to chunks at given points.

    Input arguments are a vector 'v', [ v_1, ..., v_n ] and a vector 'p'
    of split points [ p_1, ..., p_m ].

    The result is a list of vectors of type
    
        [ v_{k}, v_{k+1}, ..., v_{k+j} ]
        
    where
    
        p_{t}  <=  v_{k}, ..., v_{k+j} < p_{t+1}
        
    for some t.

    Precondition: Both 'v' and 'p' must be sorted in ascending order!
    """

    vectors = []

    last_i = 0
    for p_low in p:
        i = scipy.searchsorted(v, p_low)

        if i > last_i:
            vectors.append(v[last_i:i])
            last_i = i
     
    if last_i < len(v):
        vectors.append(v[last_i:])
   
    return vectors


def stretch_shape(arr, newshape):
    """Reshape the array, adding new axes with repeated data.

    'newshape' should be the new shape of the array, where negative
    entry -j can label the element j-1 of the old shape.
    """
    nshape = list(newshape)
    arr = _n.array(arr)

    for i, length in enumerate(nshape):
        if length < 0:
            nshape[i] = arr.shape[- length - 1]
        else:
            nshape[i] = 1

    arr = _n.reshape(arr, nshape)

    for i, length in enumerate(newshape):
        if length > 0:
            arr = _n.repeat(arr, length, i)

    return arr

def norm2(v):
    """Evaluate the 2-norm of the given vector."""
    return _n.sqrt(_n.sum(abs(_n.ravel(v)**2)))


class Sliceable(object):
    """
    Add slicing to an object.

    That is, instances of this class have a syntax ``object[slices]``
    that returns a new item of the same type, such that
    ``newitem.attr = item.attr[xslices]`` for all specified attributes
    -- where ``xslices`` is a permutation of ``slices`` specific to ``attr``.

    New items are constructed via self.__class__.__new__ and then
    updating item.__dict__. Note that data is not necessarily copied:
    for example numarray arrays return "Views" when sliced.
    """
    
    def __init__(self, axes):
        """
        Initialize sliceable behavior.

        :param axes:
            A mapping of axes of attributes to virtual axes, given in form::
            
                { 'attr_1': (a_1_1, ..., a_m_1),
                  'attr_2': (a_1_2, ..., a_m_2),
                  ... }

            Integers ``a_n_m`` in the tuples map indices of the object
            to indices of each attribute. The above indicates that
            ``obj[i_1, ..., i_m]`` corresponds to
            ``obj.attr_j[i_{a_1_j}, ..., i_{a_m_j}]``.

            If ``a_n_m`` is ``None``, then the value of the corresponding
            attribute is broadcasted along the axis. For example,
            ``a_2_1 == None`` corresponds to
            ``obj.attr_1[i_{a_1_1}, i_{a_3_1}, ..., i_{a_m_1}]``,
            for any value of ``i_2``.
        
        :type axes: dict
        """
        self.__axes = {}
        self.__naxes = 0
        for key, value in axes.iteritems():
            self.__axes[key] = self.__tuplify(value)
            self.__naxes = max(len(self.__axes[key]), self.__naxes)

    def __tuplify(value):
        try:
            return tuple(value)
        except TypeError:
            return (value,)
    __tuplify = staticmethod(__tuplify)

    def __get_slices(self, in_slices):
        in_slices = self.__tuplify(in_slices)
        slices = {}
        for key, axis in self.__axes.iteritems():
            value = getattr(self, key)
            cur_slices = [slice(None)] * len(value.shape)
            for i, sl in zip(axis, in_slices):
                if i is None:
                    continue # broadcast along this axis
                cur_slices[int(i)] = sl
            slices[key] = tuple(cur_slices)
        return slices

    def __getitem__(self, in_slices):
        """Slice the item into pieces."""
        in_slices = self.__tuplify(in_slices)
        if len(in_slices) != self.__naxes:
            raise KeyError("Invalid number of slices")
        
        item = self.__class__.__new__(self.__class__)
        item.__dict__.update(self.__dict__)
        for key, sl in self.__get_slices(in_slices).iteritems():
            setattr(item, key, getattr(self, key)[sl])
        return item

    def __setitem__(self, in_slices, other):
        """Assign to slices."""
        in_slices = self.__tuplify(in_slices)
        if len(in_slices) != self.__naxes:
            raise KeyError("Invalid number of slices")
        
        for key, sl in self.__get_slices(in_slices).iteritems():
            getattr(self, key)[sl] = getattr(other, key)

class ArrayProperty(property):
    """
    Property specifying an array attribute that will slice-assign when set,
    except for the first time.
    """
    
    def __init__(self, name, doc=None, settable=False):
        property.__init__(self, fset=self.__fset, fget=self.__fget, doc=doc)
        self.name = name
        self.settable = settable
        self.set_name = '__ArrayProperty_%s_set' % self.name
        self.__doc__ = doc

    def __fset(self, obj, value):
        if not (self.settable or self.set_name in obj.__dict__):
            obj.__dict__[self.set_name] = True
            obj.__dict__[self.name] = value
        else:
            obj.__dict__[self.name][...] = value

    def __fget(self, obj):
        return obj.__dict__[self.name]
