#ifndef FAUST_FACT_H
#define FAUST_FACT_H

#include <vector>

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
        unsigned long int parameter_sz;
};

template<typename FPP>
void prox_blockdiag(FPP* mat_data,  unsigned long mat_nrows, unsigned long mat_ncols, unsigned long *m_ptr, unsigned long *n_ptr, unsigned int vec_size, const bool normalized, const bool pos, FPP* mat_out);

template<typename FPP>
void prox_int(unsigned int cons_type, unsigned long cons_param, FPP* mat_in, unsigned long num_rows,
        unsigned long num_cols, FPP* mat_out, const bool normalized=true, const bool pos=false);

template<typename FPP, typename FPP2>
void prox_real(unsigned int cons_type, FPP2 cons_param, FPP* mat_in, unsigned long num_rows,
        unsigned long num_cols, FPP* mat_out, const bool normalized=false, const bool pos=false);

template<typename FPP>
void prox_mat(unsigned int cons_type, FPP* cons_param, FPP* mat_in, unsigned long num_rows, unsigned long num_cols, FPP* mat_out, const bool normalized=false, const bool pos=false);

template<typename FPP>
class PyxStoppingCriterion
{
    public:
        bool is_criterion_error;
        int num_its;
        FPP error_threshold;
        unsigned long max_num_its;
};

template<typename FPP, typename FPP2 = double>
class PyxParamsFact
{
    public:
        int num_facts;
        bool is_update_way_R2L;
        FPP2 init_lambda;
        FPP2 step_size;
        PyxConstraintGeneric** constraints; // size: (num_facts-1)*2
        unsigned int num_constraints;
        bool is_verbose;
        bool constant_step_size;
        unsigned int grad_calc_opt_mode;
        int norm2_max_iter;
        double norm2_threshold;
};

template<typename FPP, typename FPP2 = double>
class PyxParamsFactPalm4MSA : public PyxParamsFact<FPP,FPP2>
{
    public:
        FPP** init_facts;// num_facts elts
        unsigned long* init_fact_sizes;
        PyxStoppingCriterion<FPP2> stop_crit;
        virtual PyxStoppingCriterion<FPP2>& get_stop_crit() {return stop_crit;};

};

template<typename FPP, typename FPP2 = double>
class PyxParamsFactPalm4MSAFFT : public PyxParamsFactPalm4MSA<FPP,FPP2>
{
    public:
        FPP* init_D; // a vector for the diagonal
};

template<typename FPP, typename FPP2 = double>
class PyxParamsHierarchicalFact : public PyxParamsFact<FPP,FPP2>
{
    public:
        unsigned int num_rows;
        unsigned int num_cols;
        PyxStoppingCriterion<FPP2>* stop_crits; // must be of size 2
        bool is_fact_side_left;
        virtual const PyxStoppingCriterion<FPP2>* get_stop_crits() {return stop_crits;};
};

template<typename FPP, typename FPP2 = double>
class PyxParamsHierarchicalFactFFT : public PyxParamsHierarchicalFact<FPP,FPP2>
{
    public:
        FPP* init_D;
};

/**
 * Generic function acting on behalf of fact_palm4MSAFFT() and fact_palm4MSA().
 */
template<typename FPP, typename FPP2 = double>
FaustCoreCpp<FPP>* fact_palm4MSA(FPP*, unsigned int, unsigned int, PyxParamsFactPalm4MSA<FPP,FPP2>*, bool, FPP*);

template<typename FPP, typename FPP2 = double>
FaustCoreCpp<FPP>* fact_palm4MSAFFT(FPP*, unsigned int, unsigned int, PyxParamsFactPalm4MSAFFT<FPP,FPP2>*, bool, FPP*);

/** this function is here to factorize the code of the two functions above */
template<typename FPP, typename FPP2 = double>
FaustCoreCpp<FPP>* fact_palm4MSA_gen(FPP*, unsigned int, unsigned int, PyxParamsFactPalm4MSA<FPP,FPP2>*, bool, FPP*);

template<typename FPP, typename FPP2 = double>
FaustCoreCpp<FPP>* fact_hierarchical_gen(FPP*, FPP*, unsigned int, unsigned int, PyxParamsHierarchicalFact<FPP,FPP2>*, bool, FPP*);

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_hierarchical(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsHierarchicalFact<FPP, FPP2>* p, FPP* out_lambda);

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_hierarchical_fft(FPP* U, FPP* L, unsigned int num_rows, unsigned int num_cols, PyxParamsHierarchicalFactFFT<FPP, FPP2>* p, FPP* out_lambda_D);

template<typename FPP>
FaustCoreCpp<FPP>* hierarchical2020(FPP* mat, unsigned int num_rows, unsigned int num_cols, /*unsigned int nites*/ PyxStoppingCriterion<double>* sc, PyxConstraintGeneric** constraints, unsigned int num_cons, unsigned int num_facts, double* out_buf, bool is_update_way_R2L, bool is_fact_side_left, bool use_csr, bool packing_RL, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size, const bool on_gpu=false);

template<typename FPP>
FaustCoreCpp<FPP>* palm4msa2020(FPP* mat, unsigned int num_rows, unsigned int num_cols,  PyxConstraintGeneric** constraints, unsigned int num_cons, double* out_buf, PyxStoppingCriterion<double> sc, bool is_update_way_R2L, bool use_csr, bool packing_RL, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size, const bool on_gpu=false);

#include "FaustFact.hpp"

#endif
