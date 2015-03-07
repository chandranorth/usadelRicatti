"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

Serializing data to and from HDF5 files.
"""

from __future__ import division

import numpy as _n
import sys as _sys
import re as _re
import os as _os
import scipy as _scipy
import tables as _tbl
import shutil as _shutil
import datetime as _datetime

import hdf5pickle

from util import *

__docformat__ = "restructuredtext en"

__all__ = ['save', 'load']

######################################################################

def _typeconvert(value):
    """Work around bugs in old matlab versions on 64-bit architectures..."""
    value = _n.asarray(value)
    savetype = value.dtype
    if savetype == _n.int64:
        savetype = _n.int32
    elif savetype == _n.uint64:
        savetype = _n.uint32
    value = _n.asarray(value, dtype=savetype)
    return value

_type_map = { int: _n.int32 }
"""Again, matlab bug workaround"""

def parent_paths(path, offdepth=0):
    x = path.split('/')
    parent = '/'
    for i in range(1, len(x)-offdepth):
        child = x[i]
        parent = '/'.join(x[:i])
        full = '/'.join([parent,child])
        if parent == '':
            parent = '/'
        yield (parent, child, full)

def create_space(where, name, overwrite=True):
    """
    Create parent groups, and optionally remove an existing node from
    a HDF5 file.
    """
    file = where._v_file
    path = '/'.join([where._v_pathname, name]).replace('//', '/')

    for (parent, child, full) in parent_paths(path):
        try:
            file.getNode(full, classname='Group')
        except LookupError:
            if overwrite:
                try:
                    file.removeNode(full)
                except LookupError:
                    pass
            file.createGroup(parent, child)

    if overwrite:
        file.removeNode(path, recursive=True)

def create_group(where, name, overwrite=True):
    if isinstance(where, _tbl.File):
        where = where.root

    file = where._v_file
    path = '/'.join([where._v_pathname, name]).replace('//', '/')

    create_space(where, name, overwrite=overwrite)

    parent, realname = path.rsplit('/', 1)
    if parent == '':
        parent = '/'
    return file.createGroup(parent, realname)

def save(where, name, obj, overwrite=True):
    """
    hdf5pickle an object to given node in a HDF5 file.

    :Parameters:
      where : tables.Group : The parent node where to save
      name : str : Sub-path where to save
      obj : anything : The object to save

      overwrite : bool: Whether to overwrite the object, if it already is there

    :raises: tables.NodeError, if overwrite==False and node exists
    """
    if isinstance(where, _tbl.File):
        where = where.root

    file = where._v_file
    path = '/'.join([where._v_pathname, name]).replace('//', '/')

    create_space(where, name, overwrite=overwrite)

    if isinstance(obj, _n.ndarray):
        obj = _typeconvert(obj)
    
    hdf5pickle.dump(obj, file, path, type_map=_type_map)

def load(where, path=None):
    """
    hdf5-unpickle an object from a given node in a HDF5 file.
    """
    if isinstance(where, _tbl.File):
        where = where.root
    
    p = where._v_pathname
    if path is not None:
        p = '/'.join([p, path]).replace('//', '/')
    return hdf5pickle.load(where._v_file, p)

def find_node_in_file(filename, path=None):
    if path is None:
        components = filename.split(_os.path.sep)
        for k in reversed(xrange(len(components))):
            filename = _os.path.sep.join(components[:k])
            path = _os.path.sep.join(components[k:])
            if _os.path.isfile(filename):
                break

    if not filename:
        filename = path
        path = ''

    if not path.startswith('/'):
        path = '/' + path

    f = _tbl.openFile(filename, 'r')
    return f, f.getNode(path)

