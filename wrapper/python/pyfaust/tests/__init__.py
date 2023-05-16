import unittest
from pyfaust.tests.TestFaust import TestFaust
from pyfaust.tests.TestPoly import TestPoly


def run_tests(dev, dtype):
    """
    Runs all available tests using device dev ('cpu' or 'gpu') and scalar type
    dtype ('real' or 'complex') when it applies.
    """
    runner = unittest.TextTestRunner()
    suite = unittest.TestSuite()
    for class_name in ['TestFaust', 'TestPoly']:
        testloader = unittest.TestLoader()
        test_names = eval("testloader.getTestCaseNames("+class_name+")")
        for meth_name in test_names:
            test = eval(""+class_name+"('"+meth_name+"', dev=dev, dtype=dtype)")
            suite.addTest(test)
    runner.run(suite)
