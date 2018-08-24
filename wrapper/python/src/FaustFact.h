#ifndef FAUST_FACT_H
#define FAUST_FACT_H


class PyxConstraintGeneric
{

    public:
        int name;
        unsigned long num_rows;
        unsigned long num_cols;

        bool is_int_constraint();
        bool is_real_constraint();
        bool is_mat_constraint();
};

class PyxConstraintInt : public PyxConstraintGeneric
{
    public:
        unsigned long parameter;
};

template<typename FPP>
class PyxConstraintScalar : public PyxConstraintGeneric
{
    public:
        FPP parameter;
};

template<typename FPP>
class PyxConstraintMat : public PyxConstraintGeneric
{
    public:
        FPP* parameter;
};

template<typename FPP>
class PyxStoppingCriterion
{
    public:
        bool is_criterion_error;
        int num_its;
        FPP error_treshold;
        unsigned long max_num_its;
};

template<typename FPP>
class PyxParamsFact
{
    public:
        int num_facts;
        bool is_update_way_R2L;
        FPP init_lambda;
        FPP** init_facts;// num_facts elts
        unsigned long* init_fact_sizes;
        FPP step_size;
        PyxConstraintGeneric** constraints; // size: (num_facts-1)*2
        unsigned int num_constraints;
        bool is_verbose;
        bool constant_step_size;
};

template<typename FPP>
class PyxParamsFactPalm4MSA : public PyxParamsFact<FPP>
{
    public:
        PyxStoppingCriterion<FPP> stop_crit;
};

template<typename FPP>
class PyxParamsHierarchicalFact : public PyxParamsFact<FPP>
{
    public:
        unsigned int num_rows;
        unsigned int num_cols;
        PyxStoppingCriterion<FPP>* stop_crits; // must be of size 2
        bool is_fact_side_left;
};

template<typename FPP>
FaustCoreCpp<FPP>* fact_palm4MSA(FPP*, unsigned int, unsigned int, PyxParamsFactPalm4MSA<FPP>*, bool);

template<typename FPP>
FaustCoreCpp<FPP>* fact_hierarchical(FPP*, unsigned int, unsigned int, PyxParamsHierarchicalFact<FPP>*, bool);

#include "FaustFact.hpp"

#endif
