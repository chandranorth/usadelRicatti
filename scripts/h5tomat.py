#!/usr/bin/env python
import tables, scipy, optparse, os, sys
scipy.pkgload('io')

__revision__ = "$Id: h5tomat.py 3152 2006-09-29 10:26:22Z pauli $"

def main():
    parser = optparse.OptionParser(usage="%prog infile.h5 [outdir]")
    (options, args) = parser.parse_args()

    if len(args) < 1 or len(args) > 2:
        parser.error("Wrong number of arguments")

    infile = args[0]
    if len(args) > 1:
        outdir = args[1]
    else:
        outdir = os.path.splitext(infile)[0] + '.mat'

    f = None
    try:
        # Open
        
        try:
            f = tables.openFile(infile, 'r')
        except IOError, err:
            print err
            raise SystemExit(1)

        # Convert

        contents = {}

        def get_name(node):
            n = node._v_pathname.replace('/', '_')
            while n.startswith('_'): n = n[1:]
            return n

        def walk_tree(node):
            if isinstance(node, tables.Group):
                for name, item in node._v_children.iteritems():
                    if name.startswith('_'): continue
                    walk_tree(item)
            elif isinstance(node, tables.Array):
                contents[get_name(node)] = scipy.asarray(node.read())
            elif isinstance(node, tables.Table):
                name = get_name(node)
                for colname in node.colnames:
                    contents[name + '_' + colname] = scipy.asarray(
                        node.col(colname))

        walk_tree(f.root)

        # Fix array scalars
        
        for n in contents.iterkeys():
            if len(contents[n].shape) == 0:
                contents[n] = scipy.array([contents[n]])

        # Save

        scipy.io.savemat(outdir, contents)
    finally:
        if f: f.close()

    raise SystemExit(0)

if __name__ == "__main__":
    main()
