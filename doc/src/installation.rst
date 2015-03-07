.. index:: installing

.. _installation:

Installation
============

The following are required to be installed before this code will work:

- Python interpreter (version >= 2.4)
- Numpy library (version >= 1.0.4)
- Scipy library (version >= 0.5.0)
- HDF5 and Pytables libraries (version >= 1.3.0)
- Fortran 90 compiler and a C compiler

How to install everything is discussed below.

After the installation, you should be able to open a Python interpreter,
and type:

    .. code-block:: python
  
       import usadel1

If it doesn't report errors (``Traceback ...``), installation was
successful, and you can proceed to :ref:`scripting`.

If you want to run tests for :mod:`usadel1`, go to ``tests/``
directory and run::

    test.py

Note that this will exercise most of the functionality, compare
results to literature, and solve several "large" problems, and can
take several hours to complete.

.. index:: 
   pair: installing; Windows

Windows
-------

The most straightforward way to get all the prerequisites is to install
either (both are free for academic use)

- `Python(x,y) <http://www.pythonxy.com/>`_
- `Enthought Python distribution <http://www.enthought.com/products/epd.php>`_

The former is bundled with C and Fortran 77/90 compilers, that is, all that you
need. The latter is not, so you need to get compilers separately.  I'd
recommend using Python(x,y).

After the prerequisites are installed, edit ``install.bat`` and change
the ``--fcompiler`` switch to match your Fortran compiler, and
``--compiler`` to match your C compiler (Python(x,y) comes with
``mingw32``). Supported compiler choices are listed below.

Finally, run ``install.bat`` to compile and install the Usadel library.


.. rubric:: Test if it works

Double-click on the ``sns_minigap_vs_phi.py`` script in the
``scripts`` subdirectory. It should start calculating the size of the
minigap in an SNS junction.


.. rubric:: Supported compilers

============  ===================================================
F90 Compiler
============  ===================================================
``gnu95``     GNU Fortran **recommended**
``absoft``    Absoft Corp Fortran Compiler
``compaq``    Compaq Fortran Compiler
``g95``       G95 Fortran Compiler
``intel``     Intel Fortran Compiler for 32-bit apps
``intele``    Intel Fortran Compiler for Itanium apps
``intelem``   Intel Fortran Compiler for EM64T-based apps
``lahey``     Lahey/Fujitsu Fortran 95 Compiler
``nag``       NAGWare Fortran 95 Compiler
``pg``        Portland Group Fortran Compiler
``vast``      Pacific-Sierra Research Fortran 90 Compiler
``hpux``      HP Fortran 90 Compiler
``ibm``       IBM XL Fortran Compiler
``intelev``   Intel Visual Fortran Compiler for Itanium apps
``intelv``    Intel Visual Fortran Compiler for 32-bit apps
``mips``      MIPSpro Fortran Compiler
``sun``       Sun or Forte Fortran 95 Compiler
============  ===================================================

===========  =======================================================
C Compiler
===========  =======================================================
``bcpp``     Borland C++ Compiler
``cygwin``   Cygwin port of GNU C Compiler for Win32
``emx``      EMX port of GNU C Compiler for OS/2
``intel``    Intel C Compiler for 32-bit applications
``intele``   Intel C Itanium Compiler for Itanium-based applications
``mingw32``  Mingw32 port of GNU C Compiler for Win32
``msvc``     Microsoft Visual C++
``mwerks``   MetroWerks CodeWarrior
``unix``     standard UNIX-style compiler
===========  =======================================================


.. index::
   pair: installing; Linux

Linux 
-----

Recent Linux distributions may have the necessary packages available.



Debian-based

    Starting from Debian 4.0, including most versions of Ubuntu, you
    can easily install everything necessary::

        apt-get install python-numpy python-scipy python-tables gfortran python-dev python-numpy-dev

Fedora/Redhat

    .. code-block:: sh

        yum install numpy scipy python-devel hdf5-devel

    You'll need to compile Pytables yourself:
    get http://www.pytables.org/download/stable/tables-2.1.1.tar.gz
    and do::

        tar xzf tables-2.1.1.tar.gz
        cd tables-2.1.1
        python setup.py build
        su -c 'python setup.py install --skip-build --prefix=/usr/local'


.. rubric:: This package

To build this Usadel package, do

    .. code-block:: sh

       cd usadel1
       python setup.py config_fc --fcompiler=gnu95 --noarch build

and install with

    .. code-block:: sh

       sudo python setup.py install --skip-build --prefix=/usr/local

or:

    .. code-block:: sh

       su -c 'python setup.py install --skip-build --prefix=/usr/local'


.. rubric:: Test if it works

Try to run:

    .. code-block:: sh

       cd usadel1/scripts
       python sns_minigap_vs_phi.py

It should start calculating the size of the minigap in an SNS junction.


.. index::
   pair: installing; Unix

Unixes
------

In case prebuilt packages are not available on your platform (eg. on a
cluster where the administrators don't want to install additional
software system-wide, or the installed software is too old), you need
to build the prerequisites yourself. Unfortunately, this may take some
more effort.

SPD
^^^

The easiest way is to install SPD:

    http://code.google.com/p/spdproject/

which bundles most necessary parts in one piece that can be easily
built from sources. Also binaries are available for many platforms.

SPD does not ship with HDF5 or PyTables libraries, so you'll have
to install them separately.

.. rubric:: HDF5

Get ftp://ftp.hdfgroup.org/HDF5/current/src/hdf5-1.8.3.tar.gz

Do: (replace ``$SPD_PATH`` with the unpacked SPD directory)

    .. code-block:: sh

       cd $SPD_PATH
       tar xzf hdf5-1.8.3.tar.gz
       cd hdf5-1.8.3/
       ./configure --prefix=$SPD_PATH/local --disable-parallel
       make
       make install

.. rubric:: Pytables

Get http://www.pytables.org/download/pytables-2.1.1/tables-2.1.1.tar.gz

Install HDF5 first.

Do:

    .. code-block:: sh

       cd $SPD_PATH
       tar xzf tables-2.1.1.tar.gz
       cd tables-2.1.1
       export HDF5_DIR=$SPD_PATH/local
       $SPD_PATH/local/bin/python setup.py install

.. rubric:: This package

Do:

    .. code-block:: sh

       cd usadel1
       $SPD_PATH/local/bin/python setup.py config_fc --fcompiler=gnu95 --noarch build
       $SPD_PATH/local/bin/python setup.py install --skip-build
   

.. rubric:: Test if it works

Try to run:

    .. code-block:: sh

       cd usadel1/scripts
       export LD_LIBRARY_PATH=$SPD_PATH/local/lib
       $SPD_PATH/local/bin/python sns_minigap_vs_phi.py

It should start calculating the size of the minigap in an SNS junction.


Manual build
^^^^^^^^^^^^

If you do not want to use SPD, you can build everything manually, as
instructed below.

Note that you also need BLAS, LAPACK, and Zlib libraries. It is very
likely that these are already installed on a cluster used for
scientific computations.

When everything is installed, you can run:

    .. code-block:: sh

       $HOME/local/bin/python some-script.py

to run scripts.

In case your Unix shell is CSH, instead of:

    .. code-block:: sh

       export VARIABLE=value

write:

    .. code-block:: sh

       setenv VARIABLE value

In general, before installing or running codes, do:

    .. code-block:: sh

       export PATH=$HOME/local/bin:$PATH
       export LD_LIBRARY_PATH=$HOME/local/lib
       export LIBRARY_PATH=$HOME/local/lib
       export CPATH=$HOME/local/include


.. rubric:: zlib

Typically this is already installed (zlib.h is present), but if not,
get ftp://ftp.hdfgroup.org/lib-external/zlib/1.2/src/zlib-1.2.3.tar.gz

Do:

    .. code-block:: sh

       cd $HOME
       mkdir $HOME/local
       tar xzf zlib-1.2.3.tar.gz
       cd zlib-1.2.3
       ./configure --prefix=$HOME/local --shared
       make
       make install
   

.. rubric:: Python

Get http://python.org/ftp/python/2.5.2/Python-2.5.2.tgz

Do:

    .. code-block:: sh

       cd $HOME
       mkdir $HOME/local
       tar xzf Python-2.5.2.tgz
       cd Python-2.5.2
       ./configure --prefix=$HOME/local
       make
       make install
   

.. rubric:: HDF5

Get  ftp://ftp.hdfgroup.org/HDF5/current/src/hdf5-1.8.3.tar.gz

Do:

    .. code-block:: sh

       cd $HOME
       tar xzf hdf5-1.8.3.tar.gz
       cd hdf5-1.8.3/
       ./configure --prefix=$HOME/local --disable-parallel
       make
       make install
   

.. rubric:: Numpy

Get http://downloads.sourceforge.net/numpy/numpy-1.1.1.tar.gz?use_mirror=osdn

Do:

    .. code-block:: sh

       cd $HOME
       tar xzf numpy-1.1.1.tar.gz
       cd numpy-1.1.1
       $HOME/local/bin/python setup.py install
   

.. rubric:: Scipy

Get http://prdownloads.sourceforge.net/scipy/scipy-0.6.0.tar.gz?download

Do:
    
    .. code-block:: sh

       cd $HOME
       tar xzf scipy-0.6.0.tar.gz
       cd scipy-0.6.0
       $HOME/local/bin/python setup.py install
   
If it complains about missing BLAS or LAPACK, try to point it at your
BLAS/LAPACK libraries:

    .. code-block:: sh

       export BLAS=/usr/lib/libblas.so
       export LAPACK=/usr/lib/liblapack.so
   
the actual file names may be different on your system.


.. rubric:: Pytables

Get http://www.pytables.org/download/pytables-2.1.1/tables-2.1.1.tar.gz

Install HDF5 first.

Do:

    .. code-block:: sh

       cd $HOME
       tar xzf tables-2.1.1.tar.gz
       cd tables-2.1.1
       export HDF5_DIR=$HOME/local
       $HOME/local/bin/python setup.py install
   

.. rubric:: This package

Do:

    .. code-block:: sh

       cd usadel1
       $HOME/local/bin/python setup.py config_fc --fcompiler=gnu95 --noarch build
       $HOME/local/bin/python setup.py install --skip-build
   

.. rubric:: Test if it works

Try to run:

    .. code-block:: sh

       cd usadel1/scripts
       $HOME/local/bin/python sns_minigap_vs_phi.py
   
It should start calculating the size of the minigap in an SNS junction.
