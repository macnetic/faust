import unittest
from pyfaust.factparams import ConstraintName


class TestFactParams(unittest.TestCase):



    def __init__(self, methodName='runTest', dev='cpu', dtype='double'):
        super(TestFactParams, self).__init__(methodName)


    def test_int2str_str2int(self):
        print("ConstraintName.name_int2str & name_str2int")
        max_int_name = 0
        for attr in ConstraintName.__dict__:
            if isinstance(ConstraintName.__dict__[attr], int):
                str_name = attr.lower().replace('_', '').replace('blkdiag',
                                                                 'blockdiag')
                self.assertEqual(ConstraintName.name_int2str(ConstraintName.__dict__[attr]),
                                 str_name)
                self.assertEqual(ConstraintName.__dict__[attr],
                                 ConstraintName.str2name_int(str_name))
                if ConstraintName.__dict__[attr] > max_int_name:
                    max_int_name = ConstraintName.__dict__[attr]
        err_msg = "Invalid argument to designate a ConstraintName."
        self.assertRaisesRegex(ValueError, err_msg,
                               ConstraintName.str2name_int, "notaconstraint")
        self.assertRaisesRegex(ValueError, err_msg,
                               ConstraintName.name_int2str, max_int_name+1)
