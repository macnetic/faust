import unittest
from pyfaust.tests.TestFaust import TestFaust

dev = 'cpu'
field = 'real'

def run_tests(_dev, _field):
    global dev, field
    dev = _dev
    field = _field
    suite = unittest.makeSuite(TestFaust, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

