#!/usr/bin/env python

import unittest
import os
import sys
import traceback
import hotshot
import shutil

from testutils import *

__revision__ = "$Id: do-profile.py 3147 2006-09-29 08:07:21Z pauli $"

modules = [ os.path.splitext(os.path.basename(x))[0] for x in os.listdir('.')
            if x.endswith('.py') and x.startswith('test-') ]

if __name__ == "__main__":

    profilefilename = 'profile.prof'
    qprofilefilename = 'profile.qprof'

    mainmodule = __import__('__main__')
    for modulename in modules:
        module = __import__(modulename)
        mainmodule.__dict__.update(module.__dict__)

    if (len(sys.argv) == 2 and sys.argv[1]
          and not sys.argv[1].startswith('-')):

        try:
            os.unlink(qprofilefilename)
        except OSError:
            pass

        os.spawnlp(os.P_WAIT, 'qprof', 'qprof', '-o', qprofilefilename,
                   sys.argv[0], sys.argv[1], 'running-in-qprof')

        os.spawnlp(os.P_WAIT, 'sort', 'sort', '-n', '-k2', '-r',
                   '-o', qprofilefilename+'.tmp', qprofilefilename)
        shutil.move(qprofilefilename+'.tmp', qprofilefilename)
        
        print "=" * 79
        print ""
        sys.exit(0)
        
    elif (len(sys.argv) == 3 and sys.argv[2] == 'running-in-qprof'
          and sys.argv[1]):
        sys.argv = sys.argv[:2]

        testname = sys.argv[1]
        try:
            os.unlink(profilefilename)
        except OSError:
            pass

        print "=" * 79
        print "Profiling %s" % testname
        print "=" * 79

        def cmd():
            try:
                os.nice(10)
                main(defaultTest=testname)
            except SystemExit:
                pass
            except:
                traceback.print_exc()

        prof = hotshot.Profile(profilefilename)
        prof.runcall(cmd)

        print "=" * 79
        print "Profiling finished."
        print "  Dumping to %s ..." % profilefilename
        prof.close()
        print "  Dumping to %s ..." % qprofilefilename

    else:
        sys.argv = [sys.argv[0], "--help"]
        main()
        raise RuntimeError('Supply a test name')

