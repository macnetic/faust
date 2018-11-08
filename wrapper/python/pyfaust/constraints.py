# -*- coding: utf-8 -*-
from pyfaust import *

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
        if(not isinstance(name, ConstraintName) or not name.is_int_constraint()):
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
        if(not isinstance(name, ConstraintName) or not name.is_mat_constraint()):
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
        if(not isinstance(name, ConstraintName) or not name.is_real_constraint()):
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
