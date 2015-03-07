#!/usr/bin/env python
from __future__ import division
import unittest, optparse
import os, sys, re, imp, inspect, glob, datetime
import testutils
import pdb

__revision__ = "$Id$"

class ModuleFunctionTestCase(unittest.FunctionTestCase):
    _teardown_count = {}
    _setup_done = {}

    def __init__(self, testFunc, setUpModule=None, tearDownModule=None,
                 *a, **kw):
        self.tags = kw.pop('tags', [])
        self.tags.extend(testFunc.__dict__.get('__tags__', []))
        self.module = testFunc.__module__
        self._setUpModule = setUpModule
        self._tearDownModule = tearDownModule

        unittest.FunctionTestCase.__init__(self, testFunc, *a, **kw)
        self._name = "%s.%s" % (testFunc.__module__, testFunc.__name__)

    def setUp(self):
        if not ModuleFunctionTestCase._setup_done[self.module]:
            ModuleFunctionTestCase._setup_done[self.module] = True
            if self._setUpModule:
                self._setUpModule()
        unittest.FunctionTestCase.setUp(self)

    def tearDown(self):
        unittest.FunctionTestCase.tearDown(self)
        ModuleFunctionTestCase._teardown_count[self.module] -= 1
        if ModuleFunctionTestCase._teardown_count[self.module] == 0:
            if self._tearDownModule:
                self._tearDownModule()

    @classmethod
    def prepare(cls, tests):
        for test in tests:
            cls._teardown_count.setdefault(test.module, 0)
            cls._teardown_count[test.module] += 1
            cls._setup_done.setdefault(test.module, False)

    def id(self):
        return self._name

    def __str__(self):
        return self.id()

class MyTestResult(unittest._TextTestResult):
    def startTest(self, test):
        if self.showAll:
            self.stream.write(
                datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S] '))
        return unittest._TextTestResult.startTest(self, test)

class MyTestRunner(unittest.TextTestRunner):
    def _makeResult(self):
        return MyTestResult(self.stream, self.descriptions, self.verbosity)

def load_tests(filenames):
    """
    Load all Python modules listed in `filenames`, and
    return a test suite of all TestCases and test functions in the
    module.
    """
    suite = unittest.TestSuite()

    for filename in filenames:
        if filename.endswith('.py'): filename = filename[:-3]
        try:
            m = __import__(filename, {}, {}, [])
        except ImportError, e:
            print "Error importing %s: %s" % (filename, e)
            continue

        setup_module = m.__dict__.get('setup_module', None)
        if not inspect.isfunction(setup_module):
            setup_module = None

        teardown_module = m.__dict__.get('teardown_module', None)
        if not inspect.isfunction(teardown_module):
            teardown_module = None

        setup = m.__dict__.get('setup', None)
        if not inspect.isfunction(setup):
            setup = None

        teardown = m.__dict__.get('teardown', None)
        if not inspect.isfunction(teardown):
            teardown = None

        for name, value in m.__dict__.iteritems():
            if not hasattr(value, '__module__'):
                continue
            if value.__module__ != m.__name__:
                continue

            if inspect.isfunction(value) and name.startswith('test'):
                suite.addTest(ModuleFunctionTestCase(
                    value,
                    setUp=setup,
                    tearDown=teardown,
                    setUpModule=setup_module,
                    tearDownModule=teardown_module,
                    description=value.__doc__))

    # Unpack the test suite
    tests = []
    stack = [suite]
    while stack:
        item = stack[0]
        del stack[0]
        if isinstance(item, unittest.TestSuite):
            stack.extend(item._tests)
        else:
            tests.append(item)

    tests.sort(lambda a, b: cmp(a.id(), b.id()))

    # Return
    return tests

def get_testfiles(prefixes):
    """
    Find all test files (i.e. ``test_*.py``)
    """
    prefixes = [x.split('.')[0] for x in prefixes]
    if len(prefixes) == 1 and prefixes[0].startswith('test_'):
        return glob.glob('%s*.py' % prefixes[0])
    else:
        return glob.glob('test_*.py')

def get_matching_tests(prefixes, no_slow=False):
    tests = load_tests(get_testfiles(prefixes))
    ok_tests = []
    for test in tests:
        if prefixes:
            ok = False
            for prefix in prefixes:
                if test.id().startswith(prefix):
                    ok = True
                    break
        else:
            ok = True
        if ok and (no_slow and 'slow' in test.tags):
            print "SKIP (slow):", test
            ok = False
        if ok:
            ok_tests.append(test)
    return ok_tests

def list_tests(tests):

    print "Available tests:"

    for test in tests:
        if test.shortDescription():
            print " - %s : %s" % (test.id(), test.shortDescription())
        else:
            print " - %s" % (test.id(),)

def run_tests(tests, descriptions=False, verbosity=2, debug=False):
    suite = unittest.TestSuite()
    suite.addTests(tests)

    ModuleFunctionTestCase.prepare(tests)

    runner = MyTestRunner(stream=sys.stdout,
                          descriptions=descriptions,
                          verbosity=verbosity)
    if debug:
        pdb.set_trace()
    return runner.run(suite)

def main():
    testutils.prepare_sys_path()

    parser = optparse.OptionParser(usage="test.py [options] [test1]...")
    parser.add_option("-l", "--list", help="list available tests",
                      action="store_true", dest="list", default=False)
    parser.add_option("-v", "--verbose", help="verbose output",
                      action="count", dest="verbosity", default=0)
    parser.add_option("-d","--describe",help="show descriptions while running",
                      action="store_true", dest="descriptions", default=False)
    parser.add_option("-g","--debug",help="run tests in debugger",
                      action="store_true", dest="debug", default=False)
    parser.add_option("--no-slow",help="do not run extremely slow tests",
                      action="store_true", dest="no_slow", default=False)
    (options, args) = parser.parse_args()

    for i, arg in enumerate(args):
        if arg.endswith('.py'):
            args[i] = arg[:-2]

    tests = get_matching_tests(args, no_slow=options.no_slow)

    if not tests:
        parser.error('No tests like %s found!' % args)

    if options.list:
        list_tests(tests)
    else:
        r = run_tests(tests, options.descriptions, options.verbosity,
                      options.debug)
        if r.errors > 0 or r.failures > 0:
            sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
    os.nice(10)
    main()
