import unittest
from pyfaust.factparams import ConstraintName
from pyfaust.factparams import ParamsFact, ParamsPalm4MSAFGFT
import numpy as np


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
                self.assertEqual(ConstraintName.name_int2str(
                    ConstraintName.__dict__[attr]),
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

    def test_factor_format_str2int_int2str(self):
        print("ParamsFact.factor_format_int2str/str2int")
        formats = ['dense', 'sparse', 'dynamic']
        for i, s in enumerate(formats):
            self.assertEqual(ParamsFact.factor_format_int2str(i), s)
            self.assertEqual(ParamsFact.factor_format_int2str(s), s)
            self.assertEqual(ParamsFact.factor_format_str2int(s), i)
            self.assertEqual(ParamsFact.factor_format_str2int(i), i)
        type_error_msg = 'factor_format must be int or str'
        self.assertRaisesRegex(TypeError, type_error_msg,
                               ParamsFact.factor_format_int2str, None)
        self.assertRaisesRegex(TypeError, type_error_msg,
                               ParamsFact.factor_format_str2int, None)
        int_range_error_msg = r"factor_format as int must be in \[0, 1, 2\]"
        self.assertRaisesRegex(ValueError, int_range_error_msg,
                               ParamsFact.factor_format_str2int, 3)
        self.assertRaisesRegex(ValueError, int_range_error_msg,
                               ParamsFact.factor_format_int2str, 3)
        str_range_error_msg = "factor_format as str must be in " + \
            repr(formats).replace('[', r'\[').replace(']', r'\]')
        self.assertRaisesRegex(ValueError, str_range_error_msg,
                               ParamsFact.factor_format_str2int, 'anyformat')
        self.assertRaisesRegex(ValueError, str_range_error_msg,
                               ParamsFact.factor_format_int2str,  'anyformat')

    def test_ParamsPalm4MSAFGFT_set_init_D(self):
        print("ParamsPalm4MSAFGFT._set_init_D")
        # check errors are detected
        not_nparr_err = "init_D must be a numpy ndarray"
        not_vec_err = "init_D must be a vector."
        invalid_sz_err = "init_D must have the same size as first "
        "constraint number of rows"
        self.assertRaisesRegex(TypeError, not_nparr_err,
                               ParamsPalm4MSAFGFT._set_init_D, 'anything', 10)
        self.assertRaisesRegex(ValueError, not_vec_err,
                               ParamsPalm4MSAFGFT._set_init_D, np.ones((10,
                                                                        10)),
                               10)
        self.assertRaisesRegex(ValueError, invalid_sz_err,
                               ParamsPalm4MSAFGFT._set_init_D, np.ones((11)),
                               10)
        # check None is accepted and give default vector of ones
        self.assertTrue(np.allclose(ParamsPalm4MSAFGFT._set_init_D(None, 10),
                                    np.ones((10))))
        # check a proper vector is accepted (i.e. returned)
        v = np.random.rand(10)
        self.assertTrue(np.allclose(ParamsPalm4MSAFGFT._set_init_D(v, 10),
                                    v))
