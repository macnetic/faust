# -*- coding: utf-8 -*-
from pyfaust import *

"""
    This module provides all the classes that represent the input parameters needed

    by factorization algorithms FaustFactory.fact_palm4msa()
    FaustFactory.fact_hierarchical()
"""

class ConstraintGeneric(object):
    """
    This is the parent class for representing a factor constraint in FAµST factorization algorithms.

    This class shouldn't be instantiated, rather rely on sub-classes.

    <b/> See also: ConstraintInt, ConstraintReal, ConstraintMat, FaustFactory.fact_palm4msa, FaustFactory.fact_hierarchical.

    Attributes:
        _name: The name of the constraint applied to the factor (ConstraintName instance).
        _num_rows: The number of rows to constrain the factor to.
        _num_cols: The number of columns to constrain the factor to.
        _cons_value: The parameter value of the constraint.

    """

    def __init__(self, name, num_rows, num_cols, cons_value):
        """
        Constructs a generic constraint.

        Args:
            name: the name of the constraint applied to the factor (ConstraintName instance).
            num_rows: the number of rows to constrain the factor to.
            num_cols: the number of columns to constrain the factor to.
            cons_value: the parameter value of the constraint.

        """
        if(isinstance(name, str)):
            name = ConstraintName(name)
        self._name = name
        self._num_rows = num_rows
        self._num_cols = num_cols
        self._cons_value = cons_value

    @property
    def name(self):
        """
            Property to access the ConstraintName of the constraint.
        """
        return self._name.name

    def is_int_constraint(self):
        """
            Returns True if this constraint is a ConstraintInt, False otherwise.
        """
        return self._name.is_int_constraint()

    def is_real_constraint(self):
        """
            Returns True if this constraint is a ConstraintReal, False otherwise.
        """
        return self._name.is_real_constraint()

    def is_mat_constraint(self):
        """
            Returns True if this constraint is a ConstraintMat, False otherwise.
        """
        return self._name.is_mat_constraint()

class ConstraintInt(ConstraintGeneric):
    """
        This class represents an integer constraint on a matrix-factor.

        It constrains a matrix by its column/row-vectors sparsity or 0-norm
        (ConstraintName.SPLIN, ConstraintName.SPCOL, ConstraintName.SPLINCOL).

    """
    def __init__(self, name, num_rows, num_cols, cons_value):
        super(ConstraintInt, self).__init__(name, num_rows, num_cols, cons_value)
        if(not isinstance(cons_value, np.int)):
            raise TypeError('ConstraintInt must receive a int as cons_value '
                            'argument.')
        if(not isinstance(self._name, ConstraintName) or not self._name.is_int_constraint()):
            raise TypeError('ConstraintInt first argument must be a '
                            'ConstraintName with a int type name '
                            '(name.is_int_constraint() must return True).')

class ConstraintMat(ConstraintGeneric):

    def __init__(self, name, num_rows, num_cols, cons_value):
        super(ConstraintMat, self).__init__(name, num_rows, num_cols, cons_value)
        if(not isinstance(cons_value, np.matrix) and not isinstance(cons_value,
                                                               np.ndarray)):
            raise TypeError('ConstraintMat must receive a numpy matrix as cons_value '
                            'argument.')
        self.cons_value = float(self.cons_value)
        if(not isinstance(self._name, ConstraintName) or not self._name.is_mat_constraint()):
            raise TypeError('ConstraintMat first argument must be a '
                            'ConstraintName with a matrix type name '
                            '(name.is_mat_constraint() must return True)')

class ConstraintReal(ConstraintGeneric):
    """
        This class represents a real constraint on a matrix-factor.

        It constrains a matrix by a column/row-vector 2-norm
        (ConstraintName.NORMCOL, ConstraintName.NORMLIN).

    """
    def __init__(self, name, num_rows, num_cols, cons_value):
        """
        Constructs a real type constraint.

        Args:
            name: the name of the constraint applied to the factor. It must be
            a ConstraintName instance.
            num_rows: the number of rows to constrain the factor to.
            num_cols: the number of columns to constrain the factor to.
            cons_value: the parameter value of the constraint, it must be a
            float number that designates the 2-norm imposed to all columns (if
            name.name == ConstraintName.NORMCOL) or rows (if
            name.name == ConstraintName.NORMLIN).
        """
        super(ConstraintReal, self).__init__(name, num_rows, num_cols, cons_value)
        if(not isinstance(cons_value, np.float) and not isinstance(cons_value, np.int)):
            raise TypeError('ConstraintReal must receive a float as cons_value '
                            'argument.')
        self._cons_value = float(self._cons_value)
        if(not isinstance(self._name, ConstraintName) or not self._name.is_real_constraint()):
            raise TypeError('ConstraintReal first argument must be a '
                            'ConstraintName with a real type name '
                            '(name.is_real_constraint() must return True).')

class ConstraintName:
    """
    Attributes:
        SP: Designates a constraint on the sparsity/0-norm of a matrix.
        SPCOL: Designates a sparsity/0-norm constraint on the columns of a
        matrix.
        SPLIN: Designates a sparsity/0-norm constraint on the lines of a
        matrix.
        SPLINCOL: Designates a constraint that imposes both SPLIN and SPCOL
        constraints.
        SP_POS: Designates a constraint that imposes a SP constraints and
        besides erase the negative coefficients (it doesn't apply to complex
        matrices).
        NORMCOL: Designates a 2-norm constraint on the columns of a matrix.
        NORMLIN: Designates a 2-norm constraint on the lines of a matrix.
        CONST: Designates a constraint imposing to a matrix to be constant.

        
    """
    SP = 0 # Int Constraint
    SPCOL = 1 # Int Constraint
    SPLIN=2 # Int Constraint
    NORMCOL = 3 # Real Constraint
    SPLINCOL = 4 # Int Constraint
    CONST = 5 # Mat Constraint
    SP_POS = 6 # Int Constraint
    #BLKDIAG = 7 # ?? Constraint #TODO
    SUPP = 8 # Mat Constraint
    NORMLIN = 9 # Real Constraint

    def __init__(self, name):
        if(isinstance(name,str)):
            name = ConstraintName.str2name_int(name)
        if(not isinstance(name, np.int) or name < ConstraintName.SP or name > ConstraintName.NORMLIN):
            raise ValueError("name must be an integer among ConstraintName.SP,"
                             "ConstraintName.SPCOL, ConstraintName.NORMCOL,"
                             "ConstraintName.SPLINCOL, ConstraintName.CONST,"
                             "ConstraintName.SP_POS," # ConstraintName.BLKDIAG,
                             "ConstraintName.SUPP, ConstraintName.NORMLIN")
        self.name = name

    def is_int_constraint(self):
        return self.name in [ ConstraintName.SP, ConstraintName.SPCOL,
                             ConstraintName.SPLIN, ConstraintName.SPLINCOL,
                             ConstraintName.SP_POS ]

    def is_real_constraint(self):
        return self.name in [ ConstraintName.NORMCOL, ConstraintName.NORMLIN ]

    def is_mat_constraint(self):
        return self.name in [ConstraintName.SUPP, ConstraintName.CONST ]

    @staticmethod
    def str2name_int(_str):
        err_msg = "Invalid argument to designate a ConstraintName."
        if(not isinstance(_str, str)):
            raise ValueError(err_msg)
        if(_str == 'sp'):
            id = ConstraintName.SP
        elif(_str == 'splin'):
            id = ConstraintName.SPLIN
        elif(_str == 'spcol'):
            id = ConstraintName.SPCOL
        elif(_str == 'splincol'):
            id = ConstraintName.SPLINCOL
        elif(_str == 'sppos'):
            id = ConstraintName.SP_POS
        elif(_str == 'normcol'):
            id = ConstraintName.NORMCOL
        elif(_str == 'normlin'):
            id = ConstraintName.NORMLIN
        elif(_str == 'supp'):
            id = ConstraintName.SUPP
        elif(_str == 'const'):
            id = ConstraintName.CONST
        else:
            raise ValueError(err_msg)
        return id

class ParamsFact(object):

    def __init__(self, num_facts, is_update_way_R2L, init_lambda,
                 constraints, step_size, constant_step_size=False, is_verbose=False,
                 ):
        self.num_facts = num_facts
        self.is_update_way_R2L = is_update_way_R2L
        self.init_lambda = init_lambda
        self.step_size = step_size
        self.constraints = constraints
        self.is_verbose = is_verbose
        self.constant_step_size = constant_step_size

class ParamsHierarchicalFact(ParamsFact):

    def __init__(self, num_facts, is_update_way_R2L, init_lambda,
                 fact_constraints, res_constraints, data_num_rows,
                 data_num_cols, stop_crits,
                 step_size=10.0**-16, constant_step_size=False,
                 is_verbose=False,
                 is_fact_side_left = False):
        constraints = fact_constraints + res_constraints
        super(ParamsHierarchicalFact, self).__init__(num_facts,
                                                     is_update_way_R2L,
                                                     init_lambda,
                                                     constraints, step_size,
                                                     constant_step_size,
                                                     is_verbose)
        self.data_num_rows = data_num_rows
        self.data_num_cols = data_num_cols
        self.stop_crits = stop_crits
        self.is_fact_side_left = is_fact_side_left
        #TODO: verify number of constraints is consistent with num_facts in
        if((not isinstance(stop_crits, list) and not isinstance(stop_crits,
                                                                tuple)) or
           len(stop_crits) != 2 or
           not isinstance(stop_crits[0],StoppingCriterion) or not
           isinstance(stop_crits[1],StoppingCriterion)):
            raise TypeError('ParamsHierarchicalFact stop_crits argument must be a list/tuple of two '
                            'StoppingCriterion objects')
        if((not isinstance(constraints, list) and not isinstance(constraints,
                                                                tuple)) or
           np.array([not isinstance(constraints[i],ConstraintGeneric) for i in
                    range(0,len(constraints))]).any()):
            raise TypeError('constraints argument must be a list/tuple of '
                            'ConstraintGeneric (or subclasses) objects')

class ParamsPalm4MSA(ParamsFact):

    def __init__(self, num_facts, is_update_way_R2L, init_lambda,
                 constraints, stop_crit, init_facts=None, step_size=10.0**-16,
                 constant_step_size=False,
                 is_verbose=False):
        super(ParamsPalm4MSA, self).__init__(num_facts, is_update_way_R2L,
                                             init_lambda,
                                             constraints, step_size,
                                             constant_step_size,
                                             is_verbose)
        if(init_facts != None and (not isinstance(init_facts, list) and not isinstance(init_facts,
                                                               tuple) or
           len(init_facts) != num_facts)):
            raise ValueError('ParamsPalm4MSA init_facts argument must be a '
                             'list/tuple of '+str(num_facts)+" (num_facts) arguments.")
        else:
            self.init_facts = init_facts
        if(not isinstance(stop_crit, StoppingCriterion)):
           raise TypeError('ParamsPalm4MSA stop_crit argument must be a StoppingCriterion '
                           'object')
        self.stop_crit = stop_crit
        #TODO: verify number of constraints is consistent with num_facts

class StoppingCriterion(object):

    def __init__(self, is_criterion_error = False , error_treshold = 0.3,
                 num_its = 500,
                 max_num_its = 1000):
        self.is_criterion_error = is_criterion_error
        self.error_treshold = error_treshold
        self.num_its = num_its
        self.max_num_its = max_num_its
        #TODO: check_validity() like C++ code does
        if(is_criterion_error and num_its != 500):
            raise ValueError("It's forbidden to set a number of iterations as stopping"
                             " criterion when is_criterion_error == True.")
        elif(not is_criterion_error and (max_num_its != 1000 or error_treshold
                                         != 0.3)):
            raise ValueError("When is_criterion_error == True it's forbidden to use"
                             " other arguments than num_its argument to define "
                             "the stopping criterion.")

