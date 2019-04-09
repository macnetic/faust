# -*- coding: utf-8 -*-
# @PYFAUST_LICENSE_HEADER@
from pyfaust import *
import FaustCorePy
import sys
if sys.version_info > (3,0):
    from abc import ABC, abstractmethod
else:
    from abc import abstractmethod
    ABC = object # trick to handle py2 missing ABC
                 # but not using abstract class in py2.7

## @package pyfaust.factparams @brief The module for the parametrization of FAuST's algorithms (Palm4MSA and Hierarchical Factorization).
## <b/> See also: FaustFactory.fact_hierarchical, FaustFactory.fact_palm4msa

"""
    This module provides all the classes that represent the input parameters needed
    by factorization algorithms FaustFactory.fact_palm4msa()
    FaustFactory.fact_hierarchical()
"""

class ConstraintGeneric(ABC):
    """
    This is the parent class for representing a factor constraint in FAuST factorization algorithms.

    This class shouldn't be instantiated, rather rely on sub-classes.
    Typically, a constraint finds its place into a ParamsFact or sub-class
    instance (as a container for the factorization parameters).
    It's also possible to set a list of constraints with the ConstraintList class.

    <b/> See also: ConstraintInt, ConstraintReal, ConstraintMat,
    FaustFactory.fact_palm4msa, FaustFactory.fact_hierarchical, ParamsPalm4MSA,
    ParamsHierarchicalFact.

    Attributes:
        _name: The name of the constraint applied to the factor (ConstraintName instance).
        _num_rows: the number of columns of the constrained matrix.
        _num_cols: the number of columns of the constrained matrix.
        _cons_value: The value of the constraint.

    """

    def __init__(self, name, num_rows, num_cols, cons_value):
        """
        Constructs a generic constraint.

        Warning: This constructor shouldn't be called directly as the class is
        abstract.

        Args:
            name: The name of the constraint applied to the factor (ConstraintName instance).
            num_rows: the number of columns of the constrained matrix.
            num_cols: the number of columns of the constrained matrix.
            cons_value: The value of the constraint.

        Raises:
            TypeError: Can't instantiate abstract class ConstraintGeneric with
            abstract methods project. This exception is python 3 only, but
            this class shouldn't be instantiated in python 2.7 either.

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


    @abstractmethod
    def project(self, M):
        """
            Applies the constraint to the matrix M.

            NOTE: The project function is also called a proximal operator.

            Args:
                M: a numpy array, it must be of the same size as set in object attributes self._num_rows, self._num_cols.

            Raises:
                ValueError: if M.shape and self._num_rows, self._num_cols don't agree.
                TypeError: if M is not a numpy.ndarray

            Returns:
                The proximal operator result as a numpy array.
        """
        if(not isinstance(M, np.ndarray)):
           raise TypeError("M must be a numpy array.")
        if(M.shape[0] != self._num_rows or M.shape[1] != self._num_cols):
            raise ValueError("The dimensions must agree.")

class ConstraintInt(ConstraintGeneric):
    """
        This class represents an integer constraint on a matrix.

        It constrains a matrix by its column/row-vectors sparsity or also
        called 0-norm: ConstraintName.SPLIN, ConstraintName.SPCOL, ConstraintName.SPLINCOL.

        The other constraint names of this type are: ConstraintName.SP,
        ConstraintName.SP_POS, which both designate the 0-norm based constraint of the
        matrix with besides for SP_POS the zero replacement of all
        negative values.


    """
    def __init__(self, name, num_rows, num_cols, cons_value):
        """
            Args:
                name: must be a ConstraintName instance set with a value among SP_POS, SP, SPLIN, SPCOL, SPLINCOL (cf. ConstraintName) or it can also be one of the more handy str aliases which are respectively: 'sppos', 'sp', 'splin', 'spcol', 'splincol'.
                num_rows: the number of rows of the constrained matrix.
                num_cols: the number of columns of the constrained matrix.
                cons_value: the integer value of the constraint (the 0-norm as
                sparsity).

            Example:
                >>> from pyfaust.factparams import ConstraintInt
                >>> from numpy.random import rand
                >>> import numpy as np
                >>> cons = ConstraintInt('sppos', 10, 10, 2) # a short for ConstraintInt(ConstraintName(ConstraintName.SP_POS), 10, 10, 2)
                >>> # cons_value == 2 here and this is the 0-norm we want to constrain M to
                >>> M = rand(10,10)
                >>> np.count_nonzero(M)
                100
                >>> np.count_nonzero(cons.project(M))
                2

            <b/> See also: ConstraintGeneric.__init__
        """
        super(ConstraintInt, self).__init__(name, num_rows, num_cols, cons_value)
        if(not isinstance(cons_value, np.int)):
            raise TypeError('ConstraintInt must receive a int as cons_value '
                            'argument.')
        if(not isinstance(self._name, ConstraintName) or not self._name.is_int_constraint()):
            raise TypeError('ConstraintInt first argument must be a '
                            'ConstraintName with a int type name '
                            '(name.is_int_constraint() must return True).')

    def project(self, M):
        """
            <b/> See: ConstraintGeneric.project
        """
        super(ConstraintInt, self).project(M)
        return FaustCorePy.ConstraintIntCore.project(M, self._name.name, self._num_rows,
                                              self._num_cols, self._cons_value)


class ConstraintMat(ConstraintGeneric):
    """
        This class represents a matrix-based constraint to apply on a matrix.

    """
    def __init__(self, name, num_rows, num_cols, cons_value):
        """
        Constructs a matrix type constraint.

        Args:
            name: must be a ConstraintName instance set with a value among
            SUPP or CONST (cf. ConstraintName) or it can also be one of the
            more handy str aliases which are respectively: 'supp' and 'const'.
            num_rows: the number of rows of the constrainted matrix.
            num_cols: the number of columns of the constrained matrix.
            cons_value: the value of the constraint, it must be a numpy.array
            that defines the constraint (the matrix support for SUPP and the
            constant matrix for CONST).

        Example:
            >>> from pyfaust.factparams import ConstraintMat
            >>> from numpy.random import rand
            >>> from numpy import eye
            >>> from numpy.linalg import norm
            >>> cons = ConstraintMat('supp', 10, 10, eye(10))
            >>> M = rand(10,10)
            >>> from numpy import count_nonzero
            >>> count_nonzero(M)
            100
            >>> count_nonzero(cons.project(M))
            10
            >>> from numpy import diag
            >>> diag(M)
            array([ 0.77861201,  0.88512726,  0.56276019,  0.89159211,  0.85893333,
                           0.59919467,  0.02603014,  0.16725741,  0.20578577,  0.40803648])
            >>> diag(cons.project(M))
            array([ 0.39756071,  0.4519476 ,  0.28734638,  0.45524856,  0.43857293,
                           0.30594989,  0.01329104,  0.08540194,  0.10507459,  0.20834417])
            >>> norm(cons.project(M))
            0.99999999999999989

        """
        super(ConstraintMat, self).__init__(name, num_rows, num_cols, cons_value)
        if(not isinstance(cons_value, np.matrix) and not isinstance(cons_value,
                                                               np.ndarray)):
            raise TypeError('ConstraintMat must receive a numpy matrix as cons_value '
                            'argument.')
        self.cons_value = self._cons_value
        if(not isinstance(self._name, ConstraintName) or not self._name.is_mat_constraint()):
            raise TypeError('ConstraintMat first argument must be a '
                            'ConstraintName with a matrix type name '
                            '(name.is_mat_constraint() must return True)')

    def project(self, M):
        """
            <b/> See: ConstraintGeneric.project
        """
        super(ConstraintMat, self).project(M)
        return FaustCorePy.ConstraintMatCore.project(M, self._name.name, self._num_rows,
                                              self._num_cols, self._cons_value)


class ConstraintReal(ConstraintGeneric):
    """
        This class represents a real constraint on a matrix.

        It constrains a matrix by a column/row-vector 2-norm
        (ConstraintName.NORMCOL, ConstraintName.NORMLIN).

    """
    def __init__(self, name, num_rows, num_cols, cons_value):
        """
        Constructs a real type constraint.

        Args:
            name: must be a ConstraintName instance set with a value among NORMCOL, NORMLIN (cf. ConstraintName) or it can also be one of the more handy str aliases which are respectively: 'normcol', 'normlin'.
            num_rows: the number of columns of the constrained matrix.
            num_cols: the number of columns of the constrained matrix.
            cons_value: the parameter value of the constraint, it must be a
            float number that designates the 2-norm imposed to all columns (if
            name is ConstraintName.NORMCOL) or rows (if
            name is ConstraintName.NORMLIN).

        Example:
            >>> from pyfaust.factparams import ConstraintReal
            >>> from numpy.random import rand
            >>> from numpy.linalg import norm
            >>> cons = ConstraintReal('normcol', 10, 10, 2.) # a short for ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 10, 10, 2.)
            >>> M = rand(10,10)*10
            >>> norm(M[:,2],2)
            17.041462424512272
            >>> norm(cons.project(M)[:,2])
            2.0


            <b/> See also: ConstraintGeneric.__init__
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

    def project(self, M):
        """
        <b/> See: ConstraintGeneric.project
        """
        super(ConstraintReal, self).project(M)
        return FaustCorePy.ConstraintRealCore.project(M, self._name.name, self._num_rows,
                                              self._num_cols, self._cons_value)


class ConstraintName:
    """
    This class defines the names for the sub-types of constraints into the ConstraintGeneric hierarchy of classes.

    The table <a href="constraint.png">here</a> is a summary of the available
    constraints.

    Attributes:
        SP: Designates a constraint on the sparsity/0-norm of a matrix.
        SPCOL: Designates a sparsity/0-norm constraint on the columns of a matrix.
        SPLIN: Designates a sparsity/0-norm constraint on the rows of a matrix.
        SPLINCOL: Designates a constraint that imposes both SPLIN and SPCOL constraints (see example above for clarification).
        SP_POS: Designates a constraint that imposes a SP constraints and besides set to zero the negative coefficients (it doesn't apply to complex matrices).
        NORMCOL: Designates a 2-norm constraint on each column of a matrix.
        NORMLIN: Designates a 2-norm constraint on each row of a matrix.
        CONST: Designates a constraint imposing to a matrix to be constant.
        SUPP: Designates a constraint by a support matrix S (element-wisely multiplying the matrix to constrain to obtain a matrix for which the 2-norm equals 1, see: ConstraintMat.project()).
        name: The name of the constraint (actually an integer among the valid constants).

    Example:
        >>> # SPLINCOL Comprehensive Example
        >>> # This constraint doesn't necessarily
        >>> # lead to a  projected matrix with asked sparsity respected
        >>> # both for columns and rows
        >>> from numpy.random import rand
        >>> from numpy.linalg import norm
        >>> n = 10; m = 10; v = 2;
        >>> M = rand(10,10)
        >>> Mspcol = ConstraintInt('spcol', n, m, v).project(M)
        >>> Msplin = ConstraintInt('splin', n, m, v).project(M)
        >>> Mp = ConstraintInt('splincol', n, m, v).project(M)
        >>> Mp_ = Mspcol + np.where(Mspcol != 0, 0, Msplin) # the sum of Mspcol and Msplin minus their non-zero matrix intersection
        >>> Mp_/= norm(Mp_)
        >>> # Mp is approximately equal to Mp_
        >>> print(norm(Mp-Mp_,2)/norm(Mp_), 2)
        0.00769281206496
        >>> from numpy import count_nonzero
        >>> count_nonzero(Mp[:,1])
        3
        >>> # sparsity value v is not respected
        >>> count_nonzero(Mp_[:,1])
        3
        >>> count_nonzero(Mp_[1,:])
        2
        >>> count_nonzero(Mp[1,:])
        2
        >>> # v is respected for this row

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
        """
            Constructor of the ConstraintName object.

            Args:
                name: must be a valid constraint name (integer among the
                static constants defined in the class: ConstraintName.SP, ...).
        """
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
        """
            A delegate for ConstraintGeneric.is_int_constraint.
        """
        return self.name in [ ConstraintName.SP, ConstraintName.SPCOL,
                             ConstraintName.SPLIN, ConstraintName.SPLINCOL,
                             ConstraintName.SP_POS ]

    def is_real_constraint(self):
        """
            A delegate for ConstraintGeneric.is_real_constraint.
        """
        return self.name in [ ConstraintName.NORMCOL, ConstraintName.NORMLIN ]

    def is_mat_constraint(self):
        """
            A delegate for ConstraintGeneric.is_mat_constraint.
        """
        return self.name in [ConstraintName.SUPP, ConstraintName.CONST ]

    @staticmethod
    def str2name_int(_str):
        """
            Converts a str constraint shortname to its integer constant name equivalent.

            For example, str2name_int('sp') returns ConstraintName.SP.
        """
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

class ConstraintList(object):
    """
        An helper class for constructing a list of consistent ConstraintGeneric objects.
    """
    def __init__(self, *args):
        # constraint definition tuple
        tuple_len = 4 # name, value, nrows, ncols
        i = 0
        j = 1 # the number of processed constraints
        self.clist = []
        while(i < len(args)):
              cname = ConstraintName(args[i])
              if(i+1 > len(args)):
                raise ValueError("No value/parameter given to define the "
                                 +str(j)+"-th constraint.")
              cval = args[i+1]
              if(i+2 > len(args)):
                raise ValueError("No number of rows given to define the "
                                 +str(j)+"-th constraint.")
              nrows = args[i+2]
              if(i+3 > len(args)):
                raise ValueError("No number of columns given to define the "
                                 +str(j)+"-th constraint.")
              ncols = args[i+3]
              if(cname.is_int_constraint()):
                cons = ConstraintInt(cname, nrows, ncols, cval)
              elif(cname.is_real_constraint()):
                cons = ConstraintReal(cname, nrows, ncols, cval)
              elif(cname.is_mat_constraint()):
                cons = ConstraintMat(cname, nrows, ncols, cval)
              else:
                raise Exception(cname +" is not a valid name for a "
                                "ConstraintGeneric object")
              self.clist += [ cons ]
              i += tuple_len

    def __len__(self):
        return len(self.clist)

    def __add__(self, other):
        """
        Warning: it  returns a python list not a ConstraintList.
        """
        if(not isinstance(other, ConstraintList)):
           raise TypeError("Can't concatenate a ConstraintList with something"
                           " else.")
        return self.clist + other.clist

class ParamsFact(ABC):
    """
        The parent abstract class to represent the general factorization parameters.

        The class is the base parameters for Palm4MSA and Hierarchical
        factorization but as an abstract class it's not for direct use.

    <b/> See also  ParamsHierarchicalFact, ParamsPalm4MSA
    """
    def __init__(self, num_facts, is_update_way_R2L, init_lambda,
                 constraints, step_size, constant_step_size,
                 is_verbose):
        self.num_facts = num_facts
        self.is_update_way_R2L = is_update_way_R2L
        self.init_lambda = init_lambda
        self.step_size = step_size
        if(isinstance(constraints, ConstraintList)):
            self.constraints = constraints.clist
        else:
            self.constraints = constraints
        self.is_verbose = is_verbose
        self.constant_step_size = constant_step_size

    @abstractmethod
    def is_mat_consistent(self, M):
        if(not isinstance(M, np.ndarray)):
            raise ValueError("M must be a numpy ndarray")
        #print("M.shape=", M.shape)
        return M.shape[0] == self.constraints[0]._num_rows and \
                M.shape[1] == self.constraints[-1]._num_cols

class ParamsHierarchicalFact(ParamsFact):
    """
        The parent class to set input parameters for the hierarchical factorization algorithm.

        The class' goal is to instantiate a fully defined set of parameters
        for the algorithm. But it exists simplified parametrizations for the
        same algorithm as child classes.

        <b/> See also ParamsHierarchicalFactSquareMat,
        ParamsHierarchicalFactRectMat, FaustFactory.fact_hierarchical
    """
    def __init__(self, fact_constraints, res_constraints, stop_crit1,
                 stop_crit2, is_update_way_R2L=False, init_lambda=1.0,
                 step_size=10.0**-16, constant_step_size=False,
                 is_fact_side_left=False,
                 is_verbose=False):
        if(not isinstance(fact_constraints, list) and not
           isinstance(fact_constraints, ConstraintList)):
            raise TypeError('fact_constraints must be a list or a'
                            ' ConstraintList.')
        if(not isinstance(res_constraints, list) and not
           isinstance(res_constraints, ConstraintList)):
            raise TypeError('res_constraints must be a list or a ConstraintList.')
        if(len(fact_constraints) != len(res_constraints)):
            raise ValueError('fact_constraints and res_constraints must have'
                             ' same length.')
        num_facts = len(fact_constraints)+1
        if(is_fact_side_left):
            constraints = res_constraints + fact_constraints
        else:
            constraints = fact_constraints + res_constraints
        stop_crits = [ stop_crit1, stop_crit2 ]
        super(ParamsHierarchicalFact, self).__init__(num_facts,
                                                     is_update_way_R2L,
                                                     init_lambda,
                                                     constraints, step_size,
                                                     constant_step_size,
                                                     is_verbose)
        self.stop_crits = stop_crits
        self.is_fact_side_left = is_fact_side_left
        if((not isinstance(stop_crits, list) and not isinstance(stop_crits,
                                                                tuple)) or
           len(stop_crits) != 2 or
           not isinstance(stop_crits[0],StoppingCriterion) or not
           isinstance(stop_crits[1],StoppingCriterion)):
            raise TypeError('ParamsHierarchicalFact stop_crits argument must be'
                            ' a list/tuple of two StoppingCriterion objects')
        if((not isinstance(constraints, list) and not isinstance(constraints,
                                                                tuple)) or
           np.array([not isinstance(constraints[i],ConstraintGeneric) for i in
                    range(0,len(constraints))]).any()):
            raise TypeError('constraints argument must be a list/tuple of '
                            'ConstraintGeneric (or subclasses) objects')
        # auto-infer matrix dimension sizes according to the constraints
        if(is_fact_side_left):
            self.data_num_rows = res_constraints[-1]._num_rows
            self.data_num_cols = fact_constraints[0]._num_cols
        else:
            self.data_num_rows = constraints[0]._num_rows
            self.data_num_cols = constraints[-1]._num_cols

    def is_mat_consistent(self, M):
        if(not isinstance(M, np.ndarray)):
            raise ValueError("M must be a numpy ndarray")
        return M.shape[0] == self.data_num_rows and \
                M.shape[1] == self.data_num_cols

class ParamsHierarchicalFactSquareMat(ParamsHierarchicalFact):
    """
    The simplified parameterization class for factorizing a square matrix (of order a power of two) with the hierarchical factorization algorithm.

    This type of parameters is typically used for Hadamard matrix
    factorization.

    <b/> See also FaustFactory.fact_hierarchical, pyfaust.demo.hadamard
    """
    def __init__(self, n):
        d = 2**int(n)
        stop_crit = StoppingCriterion(num_its=30)
        super(ParamsHierarchicalFactSquareMat,
              self).__init__([ConstraintInt(ConstraintName(ConstraintName.SPLINCOL),d,d,2)
                                        for i in range(0,n-1)],
                                        [ConstraintInt(ConstraintName(ConstraintName.SPLINCOL),d,d,int(d/2.**(i+1)))
                                         for i in range(0,n-1)],
                                        stop_crit, stop_crit,
                                        is_update_way_R2L=True)

    @staticmethod
    def createParams(M, p):
        pot = np.log2(M.shape[0])
        if(pot > int(pot) or M.shape[0] != M.shape[1]):
            raise ValueError('M must be a '
                             'square matrix of order a power of '
                             'two.')
        pot = int(pot)
        return ParamsHierarchicalFactSquareMat(pot)



class ParamsHierarchicalFactRectMat(ParamsHierarchicalFact):
    """
    The simplified parameterization class for factorizing a rectangular matrix with the hierarchical factorization algorithm.

    <b/> See also FaustFactory.fact_hierarchical, pyfaust.demo.bsl
    """

    DEFAULT_P_CONST_FACT = 1.4

    def __init__(self, m, n, j, k, s, rho=0.8, P=None):
        from math import ceil
        #test args
        for arg,aname in zip([m, n, j, k, s],["m","n","j","k","s"]):
            if(not isinstance(m, int) and not isinstance(m, np.integer)):
                raise TypeError(aname+" must be an integer.")
        if(not isinstance(rho, float)):
            raise TypeError('rho must be a float')
        if(not isinstance(rho, float)):
            raise TypeError('p must be a float')
        if(not P):
            P=ParamsHierarchicalFactRectMat.DEFAULT_P_CONST_FACT*m**2
        S1_cons = ConstraintInt('spcol', m, n, k)
        S_cons = [S1_cons];
        for i in range(j-2):
            S_cons += [ ConstraintInt('sp', m, m, s*m) ]

        R_cons = []
        for i in range(j-1):
            R_cons += [ConstraintInt('sp', m, m, int(ceil(P*rho**i)))]

        stop_crit = StoppingCriterion(num_its=30)

        super(ParamsHierarchicalFactRectMat, self).__init__(S_cons, R_cons,
                                                            stop_crit,
                                                            stop_crit,
                                                            is_update_way_R2L=True,
                                                            is_fact_side_left=True)

    @staticmethod
    def createParams(M, p):
        # caller is responsible to check if name in p is really 'rectmat'
        def parse_p(p):
            # p = ('rectmat', j, k, s)
            # or p is (['rectmat', j, k, s ],{'rho':rho, P: P})
            if(isinstance(p, tuple) or isinstance(p, list)):
                if(len(p) == 2 and (isinstance(p[0], list) or isinstance(p[0],
                                                                        tuple))
                  and len(p[0]) == 4 and isinstance(p[1], dict) and 'rho' in
                   p[1].keys() and 'P' in p[1].keys()):
                    # ENOTE: concatenation instead of unpacking into list
                    # because of py2 (it would be ok for py3)
                    p = list(p[0][:])+[p[1]['rho'], p[1]['P']]
                elif(len(p) == 4 and (isinstance(p, list) or isinstance(p,
                                                                        tuple))):
                    pass #nothing to do
                else:
                    raise ValueError('The valid formats for p are: '
                                     '("rectmat",j,k,s) or '
                                     '[("rectmat",j,k,s),{"rho": rho, "P": P}]'
                                     ' with j, k, s being integers and rho and'
                                     ' P being floats')
            return p
        p = parse_p(p)
        if(not isinstance(M, np.ndarray)):
            raise TypeError('M must be a numpy.ndarray.')
        p = ParamsHierarchicalFactRectMat(M.shape[0], M.shape[1], *p[1:])
        return p

class ParamsPalm4MSA(ParamsFact):
    """
        The class is to set input parameters for the Palm4MSA algorithm.

        <b/> See also  FaustFactory.fact_palm4msa
    """

    def __init__(self, constraints, stop_crit, init_facts=None,
                 is_update_way_R2L=False, init_lambda=1.0,
                 step_size=10.0**-16,
                 constant_step_size=False,
                 is_verbose=False):
        if(not isinstance(constraints, list) and not
           isinstance(constraints, ConstraintList)):
            raise TypeError('constraints argument must be a list or a'
                            ' ConstraintList.')
        num_facts = len(constraints)
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

    def is_mat_consistent(self, M):
        return super(ParamsPalm4MSA, self).is_mat_consistent(M)

class ParamsPalm4MSAFGFT(ParamsPalm4MSA):
    """
    """
    def __init__(self, constraints, stop_crit, init_facts=None,
                 init_D=None,
                 is_update_way_R2L=False, init_lambda=1.0,
                 step_size=10.0**-16,
                 is_verbose=False):
        super(ParamsPalm4MSAFGFT, self).__init__(constraints, stop_crit,
                                                 init_facts, is_update_way_R2L,
                                                 init_lambda, step_size,
                                                 True, is_verbose)
        self.init_D = _init_init_D(init_D, self.constraints[0]._num_rows)

def _init_init_D(init_D, dim_sz):
    """
        Utility function for ParamsHierarchicalFactFGFT, ParamsPalm4MSAFGFT
    """
    def _check_init_D_is_consistent(init_D):
        if(not isinstance(init_D, np.ndarray)):
            raise ValueError("init_D must be a numpy ndarray")
        if(init_D.ndim != 1):
            raise ValueError("init_D must be a vector.")
        if(init_D.shape[0] != dim_sz):
            raise ValueError("init_D must have the same size as first "
                             "constraint's number of rows")

    if(not isinstance(init_D, np.ndarray) and init_D == None):
        init_D = np.ones(dim_sz)
    _check_init_D_is_consistent(init_D)
    return init_D



class StoppingCriterion(object):
    """
        This class defines a StoppingCriterion for the FAuST's algorithms.

        A stopping criterion can be of two kinds:
            - number of iterations,
            - error treshold for the approximation of the matrix.
    """
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
            raise ValueError("When is_criterion_error == False it's forbidden to use"
                             " other arguments than num_its argument to define "
                             "the stopping criterion.")

class ParamsFactFactory:
    """
        The factory for creating simplified FAuST hierarchical algorithm parameters (ParamsHierarchicalFact).

        Note: this factory is not related to ParamsPalm4MSA, it only creates
        ParamsHierarchicalFact instances.

        <b/> See also  ParamsHierarchicalFactRectMat,
        ParamsHierarchicalFactSquareMat, FaustFactory.fact_hierarchical()
    """
    SIMPLIFIED_PARAM_NAMES = [
        [ "squaremat", "hadamard"],
        ["rectmat", "meg"]
    ]
    SQRMAT_ID = 0
    RECTMAT_ID = 1

    @staticmethod
    def createParams(M, p):
        """
        """
        from pyfaust.factparams import \
        (ParamsHierarchicalFactSquareMat,
        ParamsHierarchicalFactRectMat)
        param_id = None
        c = ParamsFactFactory # class alias
        if(not c.is_a_valid_simplification(p)):
            raise TypeError('Invalid p to represent a simplified '
                            'parametrization.')
        param_id = c.get_simplification_name(p)
        if(param_id.lower() in c.SIMPLIFIED_PARAM_NAMES[c.SQRMAT_ID]):
            return ParamsHierarchicalFactSquareMat.createParams(M, p)
        elif(param_id.lower() in c.SIMPLIFIED_PARAM_NAMES[c.RECTMAT_ID]):
            return ParamsHierarchicalFactRectMat.createParams(M, p)
        else:
            raise ValueError("p is not a known simplified parametrization.")

    @staticmethod
    def get_simplification_name(p):
        # to be a valid simplification form
        # p must be something among:
        # 1. a str
        # 2. a list/tuple with l[0] being a str
        # 3. a list/tuple with first elt a list/tuple such that l[0][0] is a str
        max_depth=3
        l = [p]
        for i in range(max_depth):
            if((isinstance(l, list) or isinstance(l, tuple))):
                if(isinstance(l[0], str)):
                    return l[0]
                else:
                    l = l[0]
        return None

    @staticmethod
    def is_a_valid_simplification(p):
        return ParamsFactFactory.get_simplification_name(p) != None
