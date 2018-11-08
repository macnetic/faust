# -*- coding: utf-8 -*-
from pyfaust import *
from pyfaust.constraints import *

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
