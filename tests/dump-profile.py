#!/usr/bin/env python

import os
import sys
import hotshot, hotshot.stats

__revision__ = "$Id: dump-profile.py 3147 2006-09-29 08:07:21Z pauli $"

def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'dump-them':
        x = hotshot.stats.load('profile.prof')
        x.strip_dirs()
        x.sort_stats('time', 'calls')
        x.print_stats()
    else:
        print "Dumping..."
        os.system('%s dump-them > profile.txt 2>&1' % sys.argv[0])
        os.system('echo ============================= >> profile.txt')
        os.system('echo qprof >> profile.txt')
        os.system('echo ============================= >> profile.txt')
        os.system('cat profile.qprof >> profile.txt')

if __name__ == "__main__":
    main()

