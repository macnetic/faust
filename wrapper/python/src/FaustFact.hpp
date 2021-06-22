#include <iostream>
#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include "faust_Params.h"
#include "faust_ParamsFGFT.h"
#include "faust_ParamsPalm.h"
#include "faust_StoppingCriterion.h"
#include "faust_Palm4MSA.h"
#include "faust_Palm4MSAFGFT.h"
#include "faust_HierarchicalFact.h"
#include "faust_HierarchicalFactFGFT.h"
#include "faust_BlasHandle.h"
#include "faust_ConstraintGeneric.h"
#include "faust_hierarchical.h" // 2020
#include "faust_palm4msa2020.h"
#ifdef USE_GPU_MOD
#include "faust_gpu_mod_utils.h"
#include "faust_Transform_gpu.h"
#include "faust_TransformHelper_gpu.h"
#include "faust_MatDense_gpu.h"
#endif
#include "faust_butterfly.h"



bool PyxConstraintGeneric::is_int_constraint()
{
    switch(static_cast<faust_constraint_name>(this->name))
	{
		case CONSTRAINT_NAME_SP:
		case CONSTRAINT_NAME_SPCOL:
		case CONSTRAINT_NAME_SPLIN:
		case CONSTRAINT_NAME_SPLINCOL:
		case CONSTRAINT_NAME_SP_POS:
		case CONSTRAINT_NAME_SKPERM:
			return true;
		default:
			return false;
	}
}

bool PyxConstraintGeneric::is_real_constraint()
{
    switch(static_cast<faust_constraint_name>(this->name))
    {
        case CONSTRAINT_NAME_NORMCOL:
        case CONSTRAINT_NAME_NORMLIN:
            return true;
        default:
            return false;
    }
}

bool PyxConstraintGeneric::is_mat_constraint()
{
    switch(static_cast<faust_constraint_name>(this->name))
    {
        case CONSTRAINT_NAME_CONST:
        case CONSTRAINT_NAME_SUPP:
        case CONSTRAINT_NAME_BLKDIAG:
        case CONSTRAINT_NAME_TOEPLITZ:
        case CONSTRAINT_NAME_CIRC:
        case CONSTRAINT_NAME_HANKEL:
            return true;
        default:
            return false;
    }
}

template<typename FPP>
MHTPParams<Real<FPP>> convPyxMHTPParams2FaustMHTPParams(const PyxMHTPParams<Real<FPP>>& MHTP_params)
{
    Faust::MHTPParams<Real<FPP>> _MHTP_params;
    _MHTP_params.used = true;
    _MHTP_params.constant_step_size = MHTP_params.constant_step_size;
    _MHTP_params.step_size = MHTP_params.step_size;
    _MHTP_params.updating_lambda = MHTP_params.updating_lambda;
    _MHTP_params.palm4msa_period = MHTP_params.palm4msa_period;
    auto sc = MHTP_params.stop_crit;
    Faust::StoppingCriterion<Real<FPP>> sc_(sc.num_its, sc.is_criterion_error, sc.error_threshold, sc.max_num_its);
    _MHTP_params.sc = sc_;
    return _MHTP_params;
}

template<typename FPP>
int prox_blockdiag(FPP* mat_data,  unsigned long mat_nrows, unsigned long mat_ncols, unsigned long *m_ptr, unsigned long *n_ptr, unsigned int vec_size, const bool normalized, const bool pos, FPP* mat_out)
{
    int ret = 0;
    Faust::MatDense<FPP,Cpu> mat(mat_data, mat_nrows, mat_ncols);
    std::vector<unsigned long> m_vec, n_vec;
    for(int i = 0; i < vec_size; i++)
    {
        m_vec.push_back(m_ptr[i]);
        n_vec.push_back(n_ptr[i]);
    }
    try
    {
        Faust::prox_blockdiag(mat, m_vec, n_vec, normalized, pos);
        memcpy(mat_out, mat.getData(), sizeof(FPP) * mat_ncols * mat_nrows);
    }
    catch(std::domain_error & e)
    {
        ret = -1;
    }
    return ret;
}

template<typename FPP>
int prox_mat(unsigned int cons_type, FPP* cons_param, unsigned long cons_param_sz, FPP* mat_in, unsigned long num_rows, unsigned long num_cols, FPP* mat_out, const bool normalized /* deft to false */, const bool pos /* = false*/)
{
    int ret = 0;
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
    try
    {
        switch(static_cast<faust_constraint_name>(cons_type))
        {
            case CONSTRAINT_NAME_CONST: /**< Matrix equal to A ; MAT */
                // nothing to do, same mat returned
                Faust::prox_const(fmat, Faust::MatDense<FPP,Cpu>(cons_param, num_rows, num_cols), normalized, pos);
                break;
            case CONSTRAINT_NAME_BLKDIAG:
                Faust::prox_blockdiag(fmat, Faust::MatDense<FPP,Cpu>(cons_param, cons_param_sz/2, 2), normalized, pos);
                break;
            case CONSTRAINT_NAME_SUPP: /**< Matrix which support is equal to A ; MAT ; (frobenius norm 1)*/
                Faust::prox_supp(fmat, Faust::MatDense<FPP,Cpu>(cons_param, num_rows, num_cols), normalized, pos);
                break;
            case CONSTRAINT_NAME_TOEPLITZ:
                Faust::prox_toeplitz(fmat, normalized, pos); //, Faust::MatDense<FPP,Cpu>(cons_param, num_rows, num_cols), normalized, pos);
                memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
                return ret;
            case CONSTRAINT_NAME_CIRC:
                Faust::prox_circ(fmat, normalized, pos);//, Faust::MatDense<FPP,Cpu>(cons_param, num_rows, num_cols), normalized, pos);
                memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
                return ret;
            case CONSTRAINT_NAME_HANKEL:
                Faust::prox_hankel(fmat, normalized, pos);//, Faust::MatDense<FPP,Cpu>(cons_param, num_rows, num_cols), normalized, pos);
                memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
                return ret;
            default:
                throw invalid_argument("PyxConstraintMat::project() inconsistent constraint name");
        }
        memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
    }
    catch(std::domain_error & e)
    {
        ret = -1;
    }
    catch(std::invalid_argument & e)
    {
        ret = -2;
    }
    return ret;
}

template<typename FPP>
int prox_int(unsigned int cons_type, unsigned long cons_param, FPP* mat_in, unsigned long num_rows,
        unsigned long num_cols, FPP* mat_out, const bool normalized /* default to true */, const bool pos /* = false*/)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
    int ret = 0;
    try
    {
        switch(static_cast<faust_constraint_name>(cons_type))
        {
            case CONSTRAINT_NAME_SPCOL: /*!< fixed number of non zero elements per column INT (frobenius norm 1) */
                Faust::prox_spcol(fmat, (faust_unsigned_int) cons_param, normalized, pos);
                break;
            case CONSTRAINT_NAME_SPLIN: /*!< fixed number of non zero elements per line INT (frobenius norm 1) */
                Faust::prox_splin(fmat, (faust_unsigned_int) cons_param, normalized, pos);
                break;
            case CONSTRAINT_NAME_SPLINCOL:
                Faust::prox_splincol(fmat, (faust_unsigned_int) cons_param, normalized, pos);
                break;
            case CONSTRAINT_NAME_SP_POS:/**< fixed number of non zeros coefficients: INT (frobenius norm 1) */
                Faust::prox_sp_pos(fmat, (faust_unsigned_int) cons_param, normalized, pos);
                break;
            case CONSTRAINT_NAME_SP:
                Faust::prox_sp(fmat, (faust_unsigned_int) cons_param, normalized, pos);
                break;
            case CONSTRAINT_NAME_SKPERM:
                Faust::prox_skperm(fmat, (faust_unsigned_int) cons_param, normalized, pos);
                break;
            default:
                throw invalid_argument("PyxConstraintInt::project() inconsistent constraint name");
        }
        memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
    }
    catch(std::domain_error & e)
    {
        ret = -1;
    }
    catch(std::invalid_argument & e)
    {
        ret = -2;
    }
    return ret;
}

template<typename FPP, typename FPP2>
int prox_real(unsigned int cons_type, FPP2 cons_param, FPP* mat_in, unsigned long num_rows, unsigned long num_cols, FPP* mat_out, const bool normalized /* default to false */, const bool pos /* = false*/)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
    int ret = 0;
    try
    {
        switch(static_cast<faust_constraint_name>(cons_type))
        {
            case CONSTRAINT_NAME_NORMLIN:/**< 2nd norm of the lines of matrix A ; REAL  */
                Faust::prox_normlin(fmat, cons_param, normalized, pos);
                break;
            case CONSTRAINT_NAME_NORMCOL:/*!< 2nd norm of the columns of A REAL */
                Faust::prox_normcol(fmat, cons_param, normalized, pos);
                break;
            default:
                throw invalid_argument("PyxConstraintScalar::project() inconsistent constraint name");
        }
        memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
    }
    catch(std::domain_error & e)
    {
        // normatalization error
        ret = -1;
    }
    catch(std::invalid_argument & e)
    {
        ret = -2;
    }
    return ret;
}


template<typename FPP, typename FPP2>
void prepare_fact(const FPP* mat, const unsigned int num_rows, const unsigned int num_cols, const PyxParamsFact<FPP,FPP2>* p,
        /* out args */vector<const Faust::ConstraintGeneric*>& cons)
{
    Faust::ConstraintGeneric* tmp_cons;
    // the C++ core is responsible to print info in verbose mode
//    if(p->is_verbose)
//    {
//        cout << "fact_palm4MSA() mat[0]= " << mat[0] << endl;
//        cout << "p->num_facts: " << p->num_facts << endl;
//        cout << "p->is_update_way_R2L: " << p->is_update_way_R2L << endl;
//        cout << "p->init_lambda: " << p->init_lambda << endl;
//        cout << "p->step_size: " << p->step_size << endl;
//        cout << "p->is_verbose: " << p->is_verbose << endl;
//        cout << "p->constant_step_size: " << p->constant_step_size << endl;
//        cout << "p->grad_calc_opt_mode: " << p->grad_calc_opt_mode << endl; 
//        cout << "p->norm2_max_iter:" << p->norm2_max_iter << "(0 <=> default val.)" << endl;
//        cout << "p->norm2_threshold:" << p->norm2_threshold <<"(0 <=> default val.)"<<  endl;
//    }
    PyxConstraintInt* cons_int;
    PyxConstraintScalar<FPP2>* cons_real;
    PyxConstraintMat<FPP>* cons_mat;
    for(int i=0; i < p->num_constraints;i++)
    {
//        if(p->is_verbose) {
//            cout << "constraint[" << i << "]->name: " << p->constraints[i]->name << endl;
//            cout << "constraint[" << i << "]->num_rows: " << p->constraints[i]->num_rows << endl;
//            cout << "constraint[" << i << "]->num_cols: " << p->constraints[i]->num_cols << endl;
//        }
        //TODO: make ConstraintGeneric virtual and add a display() function to
        //avoid this mess and also a function to convert to Faust::Constraint*
        // corresponding object
        if(p->constraints[i]->is_int_constraint())
        {
            cons_int = static_cast<PyxConstraintInt*>(p->constraints[i]);
//            if(p->is_verbose)
//                cout << "constraint[" << i << "]->parameter: " << cons_int->parameter << endl;
            tmp_cons = new Faust::ConstraintInt<FPP,Cpu>(static_cast<faust_constraint_name>(p->constraints[i]->name), cons_int->parameter, p->constraints[i]->num_rows, p->constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else if(p->constraints[i]->is_real_constraint())
        {
            cons_real = static_cast<PyxConstraintScalar<FPP2>*>(p->constraints[i]);
//            if(p->is_verbose)
//                cout << "constraint[" << i << "]->parameter: " << cons_real->parameter << endl;
            tmp_cons = new Faust::ConstraintFPP<FPP,Cpu,FPP2>(static_cast<faust_constraint_name>(p->constraints[i]->name), cons_real->parameter, p->constraints[i]->num_rows, p->constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else if(p->constraints[i]->is_mat_constraint())
        {
            cons_mat = static_cast<PyxConstraintMat<FPP>*>(p->constraints[i]);
//            if(p->is_verbose)
//                cout << "constraint[" << i << "]->parameter: " << cons_mat->parameter[0] << endl;
            Faust::MatDense<FPP, Cpu> P;
            if(p->constraints[i]->num_rows * p->constraints[i]->num_cols == cons_mat->parameter_sz)
                P = Faust::MatDense<FPP, Cpu>(cons_mat->parameter, p->constraints[i]->num_rows, p->constraints[i]->num_cols);
            else
                P = Faust::MatDense<FPP, Cpu>(cons_mat->parameter, cons_mat->parameter_sz/2, 2);
            tmp_cons = new Faust::ConstraintMat<FPP,Cpu>(static_cast<faust_constraint_name>(p->constraints[i]->name), P, p->constraints[i]->num_rows, p->constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else
            handleError("FaustFact", "Invalid constraint.");
    }

}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_palm4MSA(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsFactPalm4MSA<FPP,FPP2>* p, FPP* out_lambda_D)
{
    return fact_palm4MSA_gen(mat, num_rows, num_cols, p, out_lambda_D);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_palm4MSAFFT(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsFactPalm4MSAFFT<FPP,FPP2>* p, FPP* out_lambda)
{
    return fact_palm4MSA_gen(mat, num_rows, num_cols, p, out_lambda);
}


template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_palm4MSA_gen(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsFactPalm4MSA<FPP,FPP2>* p, FPP* out_buf)
{
    FaustCoreCpp<FPP>* core = nullptr;
    Faust::ParamsPalm<FPP,Cpu,FPP2>* params;
    Faust::MatDense<FPP,Cpu> inMat(mat, num_rows, num_cols);
    vector<const Faust::ConstraintGeneric*> cons;
    vector<Faust::MatDense<FPP,Cpu> > initFacts;
    Faust::Palm4MSA<FPP,Cpu,FPP2>* palm;
    PyxParamsFactPalm4MSAFFT<FPP,FPP2> *p_fft = nullptr;
    Faust::BlasHandle<Cpu> blasHandle;

//    if(p->is_verbose) {
//        cout << "stop_crit.is_criterion_error: " << p->stop_crit.is_criterion_error << endl;
//        cout << "stop_crit.error_threshold: " << p->stop_crit.error_threshold << endl;
//        cout << "stop_crit.num_its: " << p->stop_crit.num_its << endl;
//        cout << "stop_crit.max_num_its: " << p->stop_crit.max_num_its << endl;
//    }
    prepare_fact(mat, num_rows, num_cols, p, cons);
    for(int i=0; i < p->num_facts;i++) {
//        if(p->is_verbose)
//        {
//            cout << "init_facts[" << i << "] ele[0][0]: " << p->init_facts[i][0] << " size: " << p->init_fact_sizes[i*2] << "x"<< p->init_fact_sizes[i*2+1] << endl;
//        }
        initFacts.push_back(Faust::MatDense<FPP, Cpu>(p->init_facts[i], p->init_fact_sizes[i*2], p->init_fact_sizes[i*2+1]));
    }
    // set all constructor arguments because they could be at non-default
    // values
    Faust::StoppingCriterion<FPP2> crit(p->stop_crit.num_its, p->stop_crit.is_criterion_error, p->stop_crit.error_threshold, p->stop_crit.max_num_its);

    if(p_fft = dynamic_cast<PyxParamsFactPalm4MSAFFT<FPP,FPP2>*>(p))
    {
        params = new ParamsPalmFGFT<FPP,Cpu,FPP2>(inMat, p->num_facts, cons, initFacts, p_fft->init_D, crit, p->is_verbose, p->is_update_way_R2L, p->init_lambda, p->step_size, static_cast<Faust::GradientCalcOptMode>(p->grad_calc_opt_mode));

        if(p->norm2_max_iter != 0)
            params->norm2_max_iter = p->norm2_max_iter;
        if(p->norm2_threshold != 0)
            params->norm2_threshold = p->norm2_threshold;
        palm = new Palm4MSAFGFT<FPP,Cpu,FPP2>(*static_cast<ParamsPalmFGFT<FPP,Cpu,FPP2>*>(params),blasHandle,true);
    }
    else {
        params = new ParamsPalm<FPP,Cpu,FPP2>(inMat, p->num_facts, cons, initFacts, crit, p->is_verbose, p->is_update_way_R2L, p->init_lambda, p->constant_step_size, p->step_size, static_cast<Faust::GradientCalcOptMode>(p->grad_calc_opt_mode));
        if(p->norm2_max_iter != 0)
            params->norm2_max_iter = p->norm2_max_iter;
        if(p->norm2_threshold != 0)
            params->norm2_threshold = p->norm2_threshold;
        palm = new Palm4MSA<FPP,Cpu,FPP2>(*params,blasHandle,true);
    }

    if(p->is_verbose) params->Display();

    try {
        palm->compute_facts();
    }
    catch(std::logic_error& e)
    {
        //intercept error like what():  Faust::Palm4MSA : compute_lambda :
        //Xhatt_Xhat_tr is too small or Xt_Xhat.trace is too big so lambda is
        //infinite
        cerr << e.what() << endl;
        return core; // core is nullptr (way to detect error in caller code)
    }

    FPP2 lambda = palm->get_lambda();

    std::vector<Faust::MatDense<FPP,Cpu> > facts;
    std::vector<Faust::MatGeneric<FPP, Cpu>*> rf_facts;
    facts=palm->get_facts();

	int sparse_weight;
	for(typename std::vector<Faust::MatDense<FPP, Cpu>>::iterator it = facts.begin(); it != facts.end(); it++)
	{
		sparse_weight = 2*it->getNonZeros()+it->getNbRow()+1;
		Faust::MatGeneric<FPP, Cpu>* M;
		if(sparse_weight < it->getNbRow()*it->getNbCol())
			M = new Faust::MatSparse<FPP, Cpu>(*it);
		else
			M = new Faust::MatDense<FPP, Cpu>(*it);
		rf_facts.push_back(M);
	}

    Faust::TransformHelper<FPP, Cpu> *th = new Faust::TransformHelper<FPP,Cpu>(rf_facts, FPP(lambda), false, false, /* internal call */ true);

    if(p->is_verbose) th->display();

    core = new FaustCoreCpp<FPP>(th);

    blasHandle.Destroy();

    for (typename std::vector<const Faust::ConstraintGeneric*>::iterator it = cons.begin() ; it != cons.end(); ++it)
        delete *it;

    if(p_fft == nullptr)
        // palm4MSA basis case, just get lambda
        *out_buf = FPP(lambda);
    else
    {
        // retrieve D matrix from Palm4MSAFFT
        // out buffer must have been allocated from outside
        dynamic_cast<Palm4MSAFGFT<FPP,Cpu,FPP2>*>(palm)->get_D(out_buf+1);
        // add lambda at the first position
        out_buf[0] = FPP(lambda);
    }

    delete params;
    delete palm;
    return core;

}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_hierarchical(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsHierarchicalFact<FPP, FPP2>* p, FPP* out_lambda)
{
    return fact_hierarchical_gen(mat, static_cast<FPP*>(nullptr), num_rows, num_cols, p, out_lambda);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_hierarchical_fft(FPP* U, FPP* L, unsigned int num_rows, unsigned int num_cols, PyxParamsHierarchicalFactFFT<FPP, FPP2>* p, FPP* out_lambda_D)
{
    return fact_hierarchical_gen(U, L, num_rows, num_cols, p, out_lambda_D);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_hierarchical_gen(FPP* mat, FPP* mat2, unsigned int num_rows, unsigned int num_cols, PyxParamsHierarchicalFact<FPP, FPP2>* p, FPP* out_buf)
{
    FaustCoreCpp<FPP>* core = nullptr;
    Faust::MatDense<FPP,Cpu> inMat(mat, num_rows, num_cols);
    vector<const Faust::ConstraintGeneric*> cons;
    vector<std::vector<const Faust::ConstraintGeneric*>> cons_;
    vector<const Faust::ConstraintGeneric*> fact_cons;
    vector<const Faust::ConstraintGeneric*> residuum_cons;
    vector<Faust::MatDense<FPP,Cpu> > initFacts_deft;
    Faust::HierarchicalFact<FPP,Cpu,FPP2>* hierFact;
    Faust::Params<FPP,Cpu,FPP2>* params;
    PyxParamsHierarchicalFactFFT<FPP,FPP2>* p_fft = nullptr;
    Faust::BlasHandle<Cpu> blasHandle;
    Faust::SpBlasHandle<Cpu> spblasHandle;
    Faust::MatDense<FPP,Cpu>* inMat2;

//    if(p->is_verbose)
//    {
//        cout << "p->num_rows: " << p->num_rows << endl;
//        cout << "p->num_cols: " << p->num_cols << endl;
//        cout << "p->is_fact_side_left: " << p->is_fact_side_left;
//        cout << "stop_crits[0].is_criterion_error: " << p->stop_crits[0].is_criterion_error << endl;
//        cout << "stop_crits[0].error_threshold: " << p->stop_crits[0].error_threshold << endl;
//        cout << "stop_crits[0].num_its: " << p->stop_crits[0].num_its << endl;
//        cout << "stop_crits[0].max_num_its: " << p->stop_crits[0].max_num_its << endl;
//        cout << "stop_crits[1].is_criterion_error: " << p->stop_crits[1].is_criterion_error << endl;
//        cout << "stop_crits[1].error_threshold: " << p->stop_crits[1].error_threshold << endl;
//        cout << "stop_crits[1].num_its: " << p->stop_crits[1].num_its << endl;
//        cout << "stop_crits[1].max_num_its: " << p->stop_crits[1].max_num_its << endl;
//        cout << "p->grad_calc_opt_mode: " << p->grad_calc_opt_mode << endl;
//
//    }
    prepare_fact(mat, num_rows, num_cols, p, cons);

    // set all constructor arguments because they could be at non-default
    // values
    Faust::StoppingCriterion<FPP2> crit0(p->stop_crits[0].num_its, p->stop_crits[0].is_criterion_error,p->stop_crits[0].error_threshold, p->stop_crits[0].max_num_its); //2 facts
    Faust::StoppingCriterion<FPP2> crit1(p->stop_crits[1].num_its, p->stop_crits[1].is_criterion_error, p->stop_crits[1].error_threshold, p->stop_crits[1].max_num_its); //global

    for(int i=0;i<p->num_facts-1;i++)
        fact_cons.push_back(cons[i]);

    for(int i=0;i<p->num_facts-1;i++)
        residuum_cons.push_back(cons[i+p->num_facts-1]);

    cons_.push_back(fact_cons);
    cons_.push_back(residuum_cons);

    if(p_fft = dynamic_cast<PyxParamsHierarchicalFactFFT<FPP,FPP2>*>(p))
    {
        inMat2 = new Faust::MatDense<FPP,Cpu>(mat2, num_rows, num_cols);
        params = new Faust::ParamsFGFT<FPP,Cpu,FPP2>(p->num_rows, p->num_cols, p->num_facts, cons_, initFacts_deft, p_fft->init_D, crit0, crit1, p->is_verbose, p->is_update_way_R2L, p->is_fact_side_left, p->init_lambda, p->constant_step_size, p->step_size, static_cast<Faust::GradientCalcOptMode>(p->grad_calc_opt_mode));
        if(p->norm2_max_iter != 0)
            params->norm2_max_iter = p->norm2_max_iter;
        if(p->norm2_threshold != 0)
            params->norm2_threshold = p->norm2_threshold;
        hierFact = new HierarchicalFactFGFT<FPP,Cpu,FPP2>(inMat, *inMat2, *(static_cast<ParamsFGFT<FPP,Cpu,FPP2>*>(params)), blasHandle, spblasHandle);
    }
    else
    {
        params = new Params<FPP,Cpu,FPP2>(p->num_rows, p->num_cols, p->num_facts, cons_, initFacts_deft, crit0, crit1, p->is_verbose, p->is_update_way_R2L, p->is_fact_side_left, p->init_lambda, p->constant_step_size, p->step_size, static_cast<Faust::GradientCalcOptMode>(p->grad_calc_opt_mode));
        if(p->norm2_max_iter != 0)
            params->norm2_max_iter = p->norm2_max_iter;
        if(p->norm2_threshold != 0)
            params->norm2_threshold = p->norm2_threshold;
        hierFact = new HierarchicalFact<FPP,Cpu,FPP2>(inMat, *params, blasHandle, spblasHandle);
    }

    if(p->is_verbose) params->Display();

    try {

        hierFact->compute_facts();
    }
    catch(std::logic_error& e)
    {
        //intercept error like what():  Faust::Palm4MSA : compute_lambda :
        //Xhatt_Xhat_tr is too small or Xt_Xhat.trace is too big so lambda is
        //infinite
        cerr << e.what() << endl;
        return core; // core is nullptr (way to detect error in caller code)
    }

    vector<Faust::MatSparse<FPP,Cpu> > facts;
    hierFact->get_facts(facts);
    FPP2 lambda = hierFact->get_lambda();
    *out_buf =  FPP(lambda);
    (facts[0]) *= FPP(lambda);
    // transform the sparse matrix into matrix pointers
    std::vector<Faust::MatGeneric<FPP,Cpu> *> list_fact_generic;
    list_fact_generic.resize(facts.size());
    for (int i=0;i<list_fact_generic.size();i++)
        list_fact_generic[i]=facts[i].Clone();

    //don't delete list_fact_generic because they are directly used in
    //TransformHelper (without copy)

    Faust::TransformHelper<FPP, Cpu>* th = new Faust::TransformHelper<FPP,Cpu>(list_fact_generic, 1.0, false, true, /* internal call */ true);

    if(p->is_verbose) th->display();
    core = new FaustCoreCpp<FPP>(th);

    blasHandle.Destroy();
    spblasHandle.Destroy();

    for (typename std::vector<const Faust::ConstraintGeneric*>::iterator it = cons.begin() ; it != cons.end(); ++it)
        delete *it;


    if(p_fft == nullptr)
        // palm4MSA basis case, just get lambda
        *out_buf = FPP(lambda);
    else
    {
        // retrieve D matrix from HierarchicalFactFFT
        // out buffer must have been allocated from outside
        dynamic_cast<HierarchicalFactFGFT<FPP,Cpu,FPP2>*>(hierFact)->get_D(out_buf+1);
        // add lambda at the first position
        out_buf[0] = FPP(lambda);
        delete inMat2;
    }

    delete params;
    delete hierFact;
    return core;

}

template<typename FPP, FDevice DEV>
TransformHelper<FPP, DEV>* palm4msa2020_gen(FPP* mat, unsigned int num_rows, unsigned int num_cols,  PyxConstraintGeneric** constraints, unsigned int num_cons, double* out_buf, PyxStoppingCriterion<Real<FPP>> sc, bool is_update_way_R2L, const int factor_format, bool packing_RL, PyxMHTPParams<Real<FPP>> &MHTP_params, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size, FaustCoreCpp<FPP, Cpu>* init_facts/*=nullptr*/)
{
    std::vector<Faust::ConstraintGeneric*> fact_cons;
    Faust::MatDense<FPP,DEV> inMat(num_rows, num_cols, mat);
    Faust::ConstraintGeneric* tmp_cons;
    PyxConstraintInt* cons_int;
    PyxConstraintScalar<double>* cons_real;
    PyxConstraintMat<FPP>* cons_mat;
    MHTPParams<Real<FPP>> _MHTP_params;
    TransformHelper<FPP,DEV> * th = nullptr, *th_times_lambda = nullptr;
    for(int i=0; i < num_cons;i++)
    {
        // corresponding object
        if(constraints[i]->is_int_constraint())
        {
            cons_int = static_cast<PyxConstraintInt*>(constraints[i]);
            //            if(is_verbose)
            //                cout << "constraint[" << i << "]->parameter: " << cons_int->parameter << endl;
            tmp_cons = new Faust::ConstraintInt<FPP,DEV>(static_cast<faust_constraint_name>(constraints[i]->name), cons_int->parameter, constraints[i]->num_rows, constraints[i]->num_cols);
            fact_cons.push_back(tmp_cons);
        }
        else if(constraints[i]->is_real_constraint())
        {
            cons_real = static_cast<PyxConstraintScalar<double>*>(constraints[i]);
            //            if(is_verbose)
            //                cout << "constraint[" << i << "]->parameter: " << cons_real->parameter << endl;
            tmp_cons = new Faust::ConstraintFPP<FPP,DEV,double>(static_cast<faust_constraint_name>(constraints[i]->name), cons_real->parameter, constraints[i]->num_rows, constraints[i]->num_cols);
            fact_cons.push_back(tmp_cons);
        }
        else if(constraints[i]->is_mat_constraint())
        {
            cons_mat = static_cast<PyxConstraintMat<FPP>*>(constraints[i]);
            //            if(is_verbose)
            //                cout << "constraint[" << i << "]->parameter: " << cons_mat->parameter[0] << endl;
            Faust::MatDense<FPP, DEV> P;
			faust_unsigned_int nrows, ncols;
            if(constraints[i]->num_rows * constraints[i]->num_cols == cons_mat->parameter_sz)
			{
				nrows = constraints[i]->num_rows;
				ncols = constraints[i]->num_cols;
			}
            else
			{
				nrows = cons_mat->parameter_sz/2;
				ncols = 2;
			}
			P = Faust::MatDense<FPP, DEV>(nrows, ncols, cons_mat->parameter);
            tmp_cons = new Faust::ConstraintMat<FPP,DEV>(static_cast<faust_constraint_name>(constraints[i]->name), P, constraints[i]->num_rows, constraints[i]->num_cols);
            fact_cons.push_back(tmp_cons);
        }
        else
            handleError("FaustFact", "Invalid constraint.");
    }

    if(MHTP_params.used)
        _MHTP_params = convPyxMHTPParams2FaustMHTPParams<FPP>(MHTP_params);

    bool free_cth = false;
    if(init_facts == nullptr)
    {
        th = new TransformHelper<FPP,DEV>();
    }
    else // init_facts has been set from the wrapper
    {
        free_cth = init_facts->make_transform(&th);
    }


    Faust::StoppingCriterion<Real<FPP>> sc0(sc.num_its, sc.is_criterion_error, sc.error_threshold, sc.max_num_its);
    Faust::palm4msa2(inMat, fact_cons, *th, out_buf[0], sc0, is_update_way_R2L,
            static_cast<FactorsFormat>(factor_format),
            packing_RL,
            _MHTP_params,
            false,
            norm2_threshold,
            norm2_max_iter,
            constant_step_size,
            step_size,
            false, /* on_gpu useless argument TODO: delete */
            is_verbose);

    for (typename std::vector</*const */Faust::ConstraintGeneric*>::iterator it = fact_cons.begin() ; it != fact_cons.end(); ++it)
        delete *it;
    
    if(free_cth)
        delete init_facts;

    return th;
}


template<typename FPP>
FaustCoreCpp<FPP>* palm4msa2020_cpu(FPP* mat, unsigned int num_rows, unsigned int num_cols,  PyxConstraintGeneric** constraints, unsigned int num_cons, double* out_buf, PyxStoppingCriterion<Real<FPP>> sc, bool is_update_way_R2L, const int factor_format, bool packing_RL, PyxMHTPParams<Real<FPP>> &MHTP_params, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size, FaustCoreCpp<FPP>* cth/*=nullptr*/)
{
    FaustCoreCpp<FPP> *core = nullptr;
    TransformHelper<FPP,Cpu> * th = nullptr, *th_times_lambda = nullptr;
    try
    {
        th = palm4msa2020_gen<FPP,Cpu>(mat, num_rows, num_cols, constraints, num_cons, out_buf, sc, is_update_way_R2L, factor_format, packing_RL, MHTP_params, norm2_max_iter, norm2_threshold, is_verbose, constant_step_size, step_size, cth);
        FPP _lambda = FPP(out_buf[0]);
        th_times_lambda = th->multiply(_lambda);
        delete th;
        th = th_times_lambda;
    }
    catch(std::logic_error& e)
    {
        cerr << e.what() << endl;
        return core;
    }

    if(is_verbose) th->display();
    core = new FaustCoreCpp<FPP>(th);

    return core;
}

#ifdef USE_GPU_MOD
template<typename FPP>
FaustCoreCpp<FPP>* palm4msa2020_gpu2(FPP* mat, unsigned int num_rows, unsigned int num_cols,  PyxConstraintGeneric** constraints, unsigned int num_cons, double* out_buf, PyxStoppingCriterion<Real<FPP>> sc, bool is_update_way_R2L, const int factor_format, bool packing_RL, PyxMHTPParams<Real<FPP>> &MHTP_params, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size, FaustCoreCpp<FPP>* cth/*=nullptr*/)
{
    FaustCoreCpp<FPP> *core = nullptr;
    TransformHelper<FPP,Cpu> * th = nullptr;
    try
    {
        auto gpu_th = palm4msa2020_gen<FPP,GPU2>(mat, num_rows, num_cols, constraints, num_cons, out_buf, sc, is_update_way_R2L, factor_format, packing_RL, MHTP_params, norm2_max_iter, norm2_threshold, is_verbose, constant_step_size, step_size, cth);
        FPP lambda = FPP(out_buf[0]);
        gpu_th->multiply(lambda);
        if(is_verbose) gpu_th->display();
		th = new Faust::TransformHelper<FPP,Cpu>();
        gpu_th->tocpu(*th);
    }
    catch(std::logic_error& e)
    {
        cerr << e.what() << endl;
        return core;
    }
    if(is_verbose) th->display();
    core = new FaustCoreCpp<FPP>(th);

    return core;
}
#endif

template<typename FPP>
FaustCoreCpp<FPP>* palm4msa2020(FPP* mat, unsigned int num_rows, unsigned int num_cols,  PyxConstraintGeneric** constraints, unsigned int num_cons, double* out_buf, PyxStoppingCriterion<Real<FPP>> sc, bool is_update_way_R2L, const int factor_format, bool packing_RL, PyxMHTPParams<Real<FPP>> &MHTP_params, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size, const bool full_gpu /*= false*/, FaustCoreCpp<FPP>* cth/*=nullptr*/)
{
#ifdef USE_GPU_MOD
    if(full_gpu)
        return palm4msa2020_gpu2(mat, num_rows, num_cols, constraints, num_cons, out_buf, sc, is_update_way_R2L, factor_format, packing_RL, MHTP_params, norm2_max_iter, norm2_threshold, is_verbose, constant_step_size, step_size, cth);
    else
#endif
        return palm4msa2020_cpu(mat, num_rows, num_cols, constraints, num_cons, out_buf, sc, is_update_way_R2L, factor_format, packing_RL, MHTP_params, norm2_max_iter, norm2_threshold, is_verbose, constant_step_size, step_size, cth);
}

template<typename FPP, FDevice DEV>
Faust::TransformHelper<FPP,DEV>* hierarchical2020_gen(FPP* mat, unsigned int num_rows, unsigned int num_cols, /* unsigned int nites*/PyxStoppingCriterion<double>* sc, PyxConstraintGeneric** constraints, unsigned int num_cons, unsigned int num_facts, double* inout_lambda, bool is_update_way_R2L, bool is_fact_side_left, const int factor_format, bool packing_RL, PyxMHTPParams<Real<FPP>>& MHTP_params, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size)
{
    Faust::MatDense<FPP,DEV> inMat(num_rows, num_cols, mat);
    vector<const Faust::ConstraintGeneric*> cons;
    vector<std::vector<const Faust::ConstraintGeneric*>> cons_;
    vector<const Faust::ConstraintGeneric*> fact_cons;
    vector<const Faust::ConstraintGeneric*> residuum_cons;
    Faust::ConstraintGeneric* tmp_cons;
    Faust::MHTPParams<Real<FPP>> _MHTP_params;

    PyxConstraintInt* cons_int;
    PyxConstraintScalar<double>* cons_real;
    PyxConstraintMat<FPP>* cons_mat;
    for(int i=0; i < num_cons;i++)
    {
        // corresponding object
        if(constraints[i]->is_int_constraint())
        {
            cons_int = static_cast<PyxConstraintInt*>(constraints[i]);
            //            if(is_verbose)
            //                cout << "constraint[" << i << "]->parameter: " << cons_int->parameter << endl;
            tmp_cons = new Faust::ConstraintInt<FPP,DEV>(static_cast<faust_constraint_name>(constraints[i]->name), cons_int->parameter, constraints[i]->num_rows, constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else if(constraints[i]->is_real_constraint())
        {
            cons_real = static_cast<PyxConstraintScalar<double>*>(constraints[i]);
            //            if(is_verbose)
            //                cout << "constraint[" << i << "]->parameter: " << cons_real->parameter << endl;
            tmp_cons = new Faust::ConstraintFPP<FPP,DEV,double>(static_cast<faust_constraint_name>(constraints[i]->name), cons_real->parameter, constraints[i]->num_rows, constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else if(constraints[i]->is_mat_constraint())
        {
            cons_mat = static_cast<PyxConstraintMat<FPP>*>(constraints[i]);
			faust_unsigned_int nrows, ncols;
            //            if(is_verbose)
            //                cout << "constraint[" << i << "]->parameter: " << cons_mat->parameter[0] << endl;
            Faust::MatDense<FPP, DEV> P;
            if(constraints[i]->num_rows * constraints[i]->num_cols == cons_mat->parameter_sz)
			{
				nrows = constraints[i]->num_rows;
				ncols = constraints[i]->num_cols;
			}
            else
			{
				nrows = cons_mat->parameter_sz/2;
				ncols = 2;
			}
			P = Faust::MatDense<FPP, DEV>(nrows, ncols, cons_mat->parameter);
            tmp_cons = new Faust::ConstraintMat<FPP, DEV>(static_cast<faust_constraint_name>(constraints[i]->name), P, constraints[i]->num_rows, constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else
            handleError("FaustFact", "Invalid constraint.");
    }

    for(int i=0;i<num_facts-1;i++)
        fact_cons.push_back(cons[i]);

    for(int i=0;i<num_facts-1;i++)
        residuum_cons.push_back(cons[i+num_facts-1]);

    if(MHTP_params.used)
        _MHTP_params = convPyxMHTPParams2FaustMHTPParams<FPP>(MHTP_params);

    if(norm2_threshold == 0)
        norm2_threshold = FAUST_PRECISION; //faust_constant.h
    if(norm2_max_iter == 0)
        norm2_max_iter = FAUST_NORM2_MAX_ITER;

    Faust::StoppingCriterion<Real<FPP>> sc0(sc[0].num_its, sc[0].is_criterion_error, sc[0].error_threshold, sc[0].max_num_its);
    Faust::StoppingCriterion<Real<FPP>> sc1(sc[1].num_its, sc[1].is_criterion_error, sc[1].error_threshold, sc[1].max_num_its);
    //        vector<Faust::StoppingCriterion<Real<FPP>>> sc_ = {sc[0].num_its, sc[1].num_its};
    vector<Faust::StoppingCriterion<Real<FPP>>> sc_ = {sc0, sc1};
    auto th = Faust::hierarchical(inMat, sc_, fact_cons, residuum_cons, inout_lambda[0], is_update_way_R2L,
            is_fact_side_left, static_cast<FactorsFormat>(factor_format),
            packing_RL,
            _MHTP_params,
            /* compute_2norm_on_array*/ false,
            norm2_threshold,
            norm2_max_iter, is_verbose, constant_step_size, step_size);
    for (typename std::vector<const Faust::ConstraintGeneric*>::iterator it = cons.begin() ; it != cons.end(); ++it)
        delete *it;
    return th;
}

template<typename FPP>
FaustCoreCpp<FPP>* hierarchical2020(FPP* mat, unsigned int num_rows, unsigned int num_cols, /* unsigned int nites*/PyxStoppingCriterion<double>* sc, PyxConstraintGeneric** constraints, unsigned int num_cons, unsigned int num_facts, double* inout_lambda, bool is_update_way_R2L, bool is_fact_side_left, const int factor_format, bool packing_RL, PyxMHTPParams<Real<FPP>> &MHTP_params, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size, const bool full_gpu /* = false*/)
{
#ifdef USE_GPU_MOD
    if(full_gpu)
        return hierarchical2020_gpu2(mat, num_rows, num_cols, sc, constraints, num_cons, num_facts, inout_lambda, is_update_way_R2L, is_fact_side_left, factor_format, packing_RL, MHTP_params, norm2_max_iter, norm2_threshold, is_verbose, constant_step_size, step_size);
    else
#endif
        return hierarchical2020_cpu(mat, num_rows, num_cols, sc, constraints, num_cons, num_facts, inout_lambda, is_update_way_R2L, is_fact_side_left, factor_format, packing_RL, MHTP_params, norm2_max_iter, norm2_threshold, is_verbose, constant_step_size, step_size);
}

template<typename FPP>
FaustCoreCpp<FPP>* hierarchical2020_cpu(FPP* mat, unsigned int num_rows, unsigned int num_cols, /* unsigned int nites*/PyxStoppingCriterion<double>* sc, PyxConstraintGeneric** constraints, unsigned int num_cons, unsigned int num_facts, double* inout_lambda, bool is_update_way_R2L, bool is_fact_side_left, const int factor_format, bool packing_RL, PyxMHTPParams<Real<FPP>> &MHTP_params, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size)
{
    FaustCoreCpp<FPP>* core = nullptr;
    Faust::TransformHelper<FPP,Cpu>* th, *th_times_lambda;
    try
    {
        th = hierarchical2020_gen<FPP,Cpu>(mat,num_rows,num_cols, sc, constraints, num_cons, num_facts, inout_lambda, is_update_way_R2L, is_fact_side_left, factor_format, packing_RL, MHTP_params, norm2_max_iter, norm2_threshold, is_verbose, constant_step_size, step_size);
        FPP lambda = FPP(inout_lambda[0]);
        th_times_lambda = th->multiply(lambda);
        delete th;
        th = th_times_lambda;
    }
    catch(std::logic_error& e)
    {
        cerr << e.what() << endl;
        return core;
    }

    if(is_verbose) th->display();
    core = new FaustCoreCpp<FPP>(th);
    return core;
}

#ifdef USE_GPU_MOD
template<typename FPP>
FaustCoreCpp<FPP>* hierarchical2020_gpu2(FPP* mat, unsigned int num_rows, unsigned int num_cols, /* unsigned int nites*/PyxStoppingCriterion<double>* sc, PyxConstraintGeneric** constraints, unsigned int num_cons, unsigned int num_facts, double* inout_lambda, bool is_update_way_R2L, bool is_fact_side_left, const int factor_format, bool packing_RL, PyxMHTPParams<Real<FPP>> &MHTP_params, unsigned int norm2_max_iter, double norm2_threshold, bool is_verbose, bool constant_step_size, double step_size)
{
    FaustCoreCpp<FPP>* core = nullptr;
    Faust::TransformHelper<FPP,Cpu>* cpu_th = nullptr, *cpu_th1 = nullptr;
    try
    {
        auto th = hierarchical2020_gen<FPP,GPU2>(mat,num_rows,num_cols, sc, constraints, num_cons, num_facts, inout_lambda, is_update_way_R2L, is_fact_side_left, factor_format, packing_RL, MHTP_params, norm2_max_iter, norm2_threshold, is_verbose, constant_step_size, step_size);
        if(is_verbose) th->display();
        Faust::TransformHelper<FPP,GPU2>* th_times_lambda = th->multiply(inout_lambda[0]);
        if(is_verbose) th->display();
        cpu_th = new Faust::TransformHelper<FPP,Cpu>();
        th_times_lambda->tocpu(*cpu_th);
        delete th;
        delete th_times_lambda;
    }
    catch(std::logic_error& e)
    {
        cerr << e.what() << endl;
        return core;
    }

    if(is_verbose) cpu_th->display();
    core = new FaustCoreCpp<FPP>(cpu_th);

    return core;
}
#endif

template<typename FPP>
FaustCoreCpp<FPP>* butterfly_hierarchical(FPP* mat, unsigned int num_rows, unsigned int num_cols, int dir)
{

    FaustCoreCpp<FPP>* core = nullptr;
    TransformHelper<FPP, Cpu> * th = nullptr;
    Faust::MatDense<FPP, Cpu> inMat(num_rows, num_cols, mat);
    th = Faust::butterfly_hierarchical(inMat, static_cast<Faust::ButterflyFactDir>(dir));
    core = new FaustCoreCpp<FPP>(th);
    return core;
}
