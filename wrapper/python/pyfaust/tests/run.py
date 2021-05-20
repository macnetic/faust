import unittest
from pyfaust.tests.TestFaust import TestFaust
import sys

if __name__ == "__main__":
    nargs = len(sys.argv)
    if(nargs > 1):
        dev = sys.argv[1]
        if dev != 'cpu' and not dev.startswith('gpu'):
            raise ValueError("dev argument must be cpu or gpu.")
        if(nargs > 2):
            field = sys.argv[2]
            if field not in ['complex', 'real']:
                raise ValueError("field must be complex or float")
        del sys.argv[2]  # deleted to avoid interfering with unittest
        del sys.argv[1]
    if(len(sys.argv) > 1):
        # ENOTE: test only a single test if name passed on command line
        singleton = unittest.TestSuite()
        singleton.addTest(TestFaust(sys.argv[1]))
        unittest.TextTestRunner().run(singleton)
    else:
        unittest.main()
