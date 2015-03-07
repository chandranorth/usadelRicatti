#!/usr/bin/env python
import os, sys

# compiler & distutils setup
# --------------------------

try:
    import setuptools
except ImportError:
    pass

from numpy.distutils.misc_util import Configuration
from distutils.util import get_platform

if not hasattr(sys, 'version_info') or sys.version_info < (2,2,0,'alpha',0):
    raise SystemExit, "Python 2.2 or later required to build this module."

# usadel1
# -------

def configuration(parent_package='', top_path=None):
    ext_sources = """
    src/usadel1/_solvercore.pyf

    src/third-party/blas/daxpy.f
    src/third-party/blas/dcopy.f
    src/third-party/blas/ddot.f
    src/third-party/blas/dscal.f
    src/third-party/blas/dswap.f
    src/third-party/blas/idamax.f
    src/third-party/linpack/dgefa.f
    src/third-party/linpack/dgesl.f
    src/third-party/slatec/dintrv.f
    src/third-party/slatec/d1mach.f
    src/third-party/bvp/colnew.f90
    src/third-party/bvp/twpbvpc.f90

    src/usadel1/lazy_alloc.f90
    src/usadel1/miscmath.f90
    src/usadel1/interpolate.f90
    src/usadel1/params.f90
    src/usadel1/error.f90
    src/usadel1/kin_equations.f90
    src/usadel1/kin_solve.f90
    src/usadel1/kin_solve2.f90
    src/usadel1/kin_solve3.f90
    src/usadel1/sp_equations.f90
    src/usadel1/sp_solve.f90
    src/usadel1/sp_solve2.f90
    
    src/usadel1/solvercore.f90
    """.split()

    version = "0.2.4"
    
    generate_version(version, os.path.dirname(__file__))

    config = Configuration('usadel1', parent_package, top_path,
                           package_path='src/usadel1',
                           install_requires=["numpy >= 1.0",
                                             "scipy >= 0.5",],
                           version=version,
                           )
    config.add_extension('_solvercore', sources=ext_sources)
    config.add_subpackage('hdf5pickle', 'src/hdf5pickle')
    config.add_data_dir('tests')

    return config

def generate_version(version, top_path):
    version_file = os.path.join(top_path, 'src', 'usadel1', 'version.py')
    hg_path = os.path.join(top_path, '.hg')
    git_path = os.path.join(top_path, '.git')
    if os.path.isdir(git_path):
        import subprocess
        try:
            p = subprocess.Popen(['git', 'rev-parse', 'HEAD'],
                                 stdout=subprocess.PIPE)
            rev = p.communicate()[0].strip()[:6]
            if '.dev' in version:
                version += "-" + rev
        except:
            pass
    elif os.path.isdir(hg_path):
        import subprocess
        try:
            p = subprocess.Popen(['hg', 'id'], stdout=subprocess.PIPE)
            hgid = p.communicate()[0].split()
            p = subprocess.Popen(['hg', 'id', '-n'], stdout=subprocess.PIPE)
            hgnum = p.communicate()[0].strip().replace('+', '')
            if len(hgid) == 2:
                if hgid[1] != version and hgid[1][0] in '0123456789':
                    print "HG version does not match setup.py, overriding!"
                    version = hgid[1]
            rev = hgid[0].replace('+', 'x')
            if rev.endswith('x') and '.dev' not in version:
                print "HG version is modified, adding .dev"
                version += ".dev"
            if '.dev' in version:
                version += hgnum + "." + rev
        except:
            pass
    f = open(version_file, 'w')
    f.write('__version__ = "%s"\n' % version)
    f.close()


# setup
# -----

if __name__ == "__main__":
    from numpy.distutils.core import setup
    cfg = configuration(top_path='').todict()
    setup(**cfg)
