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

using namespace std;
using namespace Faust;

bool PyxConstraintGeneric::is_int_constraint()
{
    switch(static_cast<faust_constraint_name>(this->name))
    {
        case CONSTRAINT_NAME_SP:
        case CONSTRAINT_NAME_SPCOL:
        case CONSTRAINT_NAME_SPLIN:
        case CONSTRAINT_NAME_SPLINCOL:
        case CONSTRAINT_NAME_SP_POS:
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
void prox_blockdiag(FPP* mat_data,  unsigned long mat_nrows, unsigned long mat_ncols, unsigned long *m_ptr, unsigned long *n_ptr, unsigned int vec_size, const bool normalized, const bool pos, FPP* mat_out)
{
    Faust::MatDense<FPP,Cpu> mat(mat_data, mat_nrows, mat_ncols);
    std::vector<unsigned long> m_vec, n_vec;
    for(int i = 0; i < vec_size; i++)
    {
        m_vec.push_back(m_ptr[i]);
        n_vec.push_back(n_ptr[i]);
    }
    Faust::prox_blockdiag(mat, m_vec, n_vec, normalized, pos);
    memcpy(mat_out, mat.getData(), sizeof(FPP) * mat_ncols * mat_nrows);
}

template<typename FPP>
void prox_mat(unsigned int cons_type, FPP* cons_param, unsigned long cons_param_sz, FPP* mat_in, unsigned long num_rows, unsigned long num_cols, FPP* mat_out, const bool normalized /* deft to false */, const bool pos /* = false*/)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
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
            return;
            break;
        case CONSTRAINT_NAME_CIRC:
            Faust::prox_circ(fmat, normalized, pos);//, Faust::MatDense<FPP,Cpu>(cons_param, num_rows, num_cols), normalized, pos);
            memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
            return;
            break;
        case CONSTRAINT_NAME_HANKEL:
            Faust::prox_hankel(fmat, normalized, pos);//, Faust::MatDense<FPP,Cpu>(cons_param, num_rows, num_cols), normalized, pos);
            memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
            return;
            break;
		default:
			throw invalid_argument("PyxConstraintMat::project() inconsistent constraint name");
	}
    memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
}

template<typename FPP>
void prox_int(unsigned int cons_type, unsigned long cons_param, FPP* mat_in, unsigned long num_rows,
        unsigned long num_cols, FPP* mat_out, const bool normalized /* default to true */, const bool pos /* = false*/)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
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
        default:
            throw invalid_argument("PyxConstraintInt::project() inconsistent constraint name");
    }
    memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
}

template<typename FPP, typename FPP2>
void prox_real(unsigned int cons_type, FPP2 cons_param, FPP* mat_in, unsigned long num_rows, unsigned long num_cols, FPP* mat_out, const bool normalized /* default to false */, const bool pos /* = false*/)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
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

template<typename FPP, typename FPP2>
void prepare_fact(const FPP* mat, const unsigned int num_rows, const unsigned int num_cols, const PyxParamsFact<FPP,FPP2>* p,
        /* out args */vector<const Faust::ConstraintGeneric*>& cons)
{
    Faust::ConstraintGeneric* tmp_cons;
    if(p->is_verbose)
    {
        cout << "fact_palm4MSA() mat[0]= " << mat[0] << endl;
        cout << "p->num_facts: " << p->num_facts << endl;
        cout << "p->is_update_way_R2L: " << p->is_update_way_R2L << endl;
        cout << "p->init_lambda: " << p->init_lambda << endl;
        cout << "p->step_size: " << p->step_size << endl;
        cout << "p->is_verbose: " << p->is_verbose << endl;
        cout << "p->constant_step_size: " << p->constant_step_size << endl;
        cout << "p->grad_calc_opt_mode: " << p->grad_calc_opt_mode << endl; 
        cout << "p->norm2_max_iter:" << p->norm2_max_iter << "(0 <=> default val.)" << endl;
        cout << "p->norm2_threshold:" << p->norm2_threshold <<"(0 <=> default val.)"<<  endl;
    }
    PyxConstraintInt* cons_int;
    PyxConstraintScalar<FPP2>* cons_real;
    PyxConstraintMat<FPP>* cons_mat;
    for(int i=0; i < p->num_constraints;i++)
    {
        if(p->is_verbose) {
            cout << "constraint[" << i << "]->name: " << p->constraints[i]->name << endl;
            cout << "constraint[" << i << "]->num_rows: " << p->constraints[i]->num_rows << endl;
            cout << "constraint[" << i << "]->num_cols: " << p->constraints[i]->num_cols << endl;
        }
        //TODO: make ConstraintGeneric virtual and add a display() function to
        //avoid this mess and also a function to convert to Faust::Constraint*
        // corresponding object
        if(p->constraints[i]->is_int_constraint())
        {
            cons_int = static_cast<PyxConstraintInt*>(p->constraints[i]);
            if(p->is_verbose)
                cout << "constraint[" << i << "]->parameter: " << cons_int->parameter << endl;
            tmp_cons = new Faust::ConstraintInt<FPP,Cpu>(static_cast<faust_constraint_name>(p->constraints[i]->name), cons_int->parameter, p->constraints[i]->num_rows, p->constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else if(p->constraints[i]->is_real_constraint())
        {
            cons_real = static_cast<PyxConstraintScalar<FPP2>*>(p->constraints[i]);
            if(p->is_verbose)
                cout << "constraint[" << i << "]->parameter: " << cons_real->parameter << endl;
            tmp_cons = new Faust::ConstraintFPP<FPP,Cpu,FPP2>(static_cast<faust_constraint_name>(p->constraints[i]->name), cons_real->parameter, p->constraints[i]->num_rows, p->constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else if(p->constraints[i]->is_mat_constraint())
        {
            cons_mat = static_cast<PyxConstraintMat<FPP>*>(p->constraints[i]);
            if(p->is_verbose)
                cout << "constraint[" << i << "]->parameter: " << cons_mat->parameter[0] << endl;
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

    if(p->is_verbose) {
        cout << "stop_crit.is_criterion_error: " << p->stop_crit.is_criterion_error << endl;
        cout << "stop_crit.error_treshold: " << p->stop_crit.error_treshold << endl;
        cout << "stop_crit.num_its: " << p->stop_crit.num_its << endl;
        cout << "stop_crit.max_num_its: " << p->stop_crit.max_num_its << endl;
    }
    prepare_fact(mat, num_rows, num_cols, p, cons);
    for(int i=0; i < p->num_facts;i++) {
        if(p->is_verbose)
        {
            cout << "init_facts[" << i << "] ele[0][0]: " << p->init_facts[i][0] << " size: " << p->init_fact_sizes[i*2] << "x"<< p->init_fact_sizes[i*2+1] << endl;
        }
        initFacts.push_back(Faust::MatDense<FPP, Cpu>(p->init_facts[i], p->init_fact_sizes[i*2], p->init_fact_sizes[i*2+1]));
    }
    // set all constructor arguments because they could be at non-default
    // values
    Faust::StoppingCriterion<FPP2> crit(p->stop_crit.num_its, p->stop_crit.is_criterion_error, p->stop_crit.error_treshold, p->stop_crit.max_num_its);

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
    std::vector<Faust::MatGeneric<FPP, Cpu>*> sp_facts;
    facts=palm->get_facts();

    for(typename std::vector<Faust::MatDense<FPP, Cpu>>::iterator it = facts.begin(); it != facts.end(); it++)
    {
        Faust::MatSparse<FPP, Cpu> * M = new Faust::MatSparse<FPP, Cpu>(*it);
        sp_facts.push_back(M);
    }

    Faust::TransformHelper<FPP, Cpu> *th = new Faust::TransformHelper<FPP,Cpu>(sp_facts, FPP(lambda), false, true, /* internal call */ true);

    for(typename std::vector<Faust::MatGeneric<FPP,Cpu>*>::iterator it = sp_facts.begin(); it != sp_facts.end(); it++)
    {
        delete *it;
    }
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

    if(p->is_verbose)
    {
        cout << "p->num_rows: " << p->num_rows << endl;
        cout << "p->num_cols: " << p->num_cols << endl;
        cout << "p->is_fact_side_left: " << p->is_fact_side_left;
        cout << "stop_crits[0].is_criterion_error: " << p->stop_crits[0].is_criterion_error << endl;
        cout << "stop_crits[0].error_treshold: " << p->stop_crits[0].error_treshold << endl;
        cout << "stop_crits[0].num_its: " << p->stop_crits[0].num_its << endl;
        cout << "stop_crits[0].max_num_its: " << p->stop_crits[0].max_num_its << endl;
        cout << "stop_crits[1].is_criterion_error: " << p->stop_crits[1].is_criterion_error << endl;
        cout << "stop_crits[1].error_treshold: " << p->stop_crits[1].error_treshold << endl;
        cout << "stop_crits[1].num_its: " << p->stop_crits[1].num_its << endl;
        cout << "stop_crits[1].max_num_its: " << p->stop_crits[1].max_num_its << endl;
        cout << "p->grad_calc_opt_mode: " << p->grad_calc_opt_mode << endl;

    }
    prepare_fact(mat, num_rows, num_cols, p, cons);

    // set all constructor arguments because they could be at non-default
    // values
    Faust::StoppingCriterion<FPP2> crit0(p->stop_crits[0].num_its, p->stop_crits[0].is_criterion_error,p->stop_crits[0].error_treshold, p->stop_crits[0].max_num_its); //2 facts
    Faust::StoppingCriterion<FPP2> crit1(p->stop_crits[1].num_its, p->stop_crits[1].is_criterion_error, p->stop_crits[1].error_treshold, p->stop_crits[1].max_num_its); //global

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

template<typename FPP>
FaustCoreCpp<FPP>* hierarchical2020(FPP* mat, unsigned int num_rows, unsigned int num_cols, /* unsigned int nites*/PyxStoppingCriterion<double>* sc, PyxConstraintGeneric** constraints, unsigned int num_cons, unsigned int num_facts, double* inout_lambda, bool is_update_way_R2L, bool is_fact_side_left, bool use_csr, bool packing_RL, unsigned int norm2_max_iter, double norm2_threshold)
{
    FaustCoreCpp<FPP>* core = nullptr;
    Faust::MatDense<FPP,Cpu> inMat(mat, num_rows, num_cols);
    vector<const Faust::ConstraintGeneric*> cons;
    vector<std::vector<const Faust::ConstraintGeneric*>> cons_;
    vector<const Faust::ConstraintGeneric*> fact_cons;
    vector<const Faust::ConstraintGeneric*> residuum_cons;
    vector<Faust::MatDense<FPP,Cpu> > initFacts_deft;
    Faust::ConstraintGeneric* tmp_cons;

    PyxConstraintInt* cons_int;
    PyxConstraintScalar<double>* cons_real;
    PyxConstraintMat<FPP>* cons_mat;
    bool is_verbose = false; // TODO: should be an argument
    for(int i=0; i < num_cons;i++)
    {
        // corresponding object
        if(constraints[i]->is_int_constraint())
        {
            cons_int = static_cast<PyxConstraintInt*>(constraints[i]);
            if(is_verbose)
                cout << "constraint[" << i << "]->parameter: " << cons_int->parameter << endl;
            tmp_cons = new Faust::ConstraintInt<FPP,Cpu>(static_cast<faust_constraint_name>(constraints[i]->name), cons_int->parameter, constraints[i]->num_rows, constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else if(constraints[i]->is_real_constraint())
        {
            cons_real = static_cast<PyxConstraintScalar<double>*>(constraints[i]);
            if(is_verbose)
                cout << "constraint[" << i << "]->parameter: " << cons_real->parameter << endl;
            tmp_cons = new Faust::ConstraintFPP<FPP,Cpu,double>(static_cast<faust_constraint_name>(constraints[i]->name), cons_real->parameter, constraints[i]->num_rows, constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else if(constraints[i]->is_mat_constraint())
        {
            cons_mat = static_cast<PyxConstraintMat<FPP>*>(constraints[i]);
            if(is_verbose)
                cout << "constraint[" << i << "]->parameter: " << cons_mat->parameter[0] << endl;
            Faust::MatDense<FPP, Cpu> P;
            if(constraints[i]->num_rows * constraints[i]->num_cols == cons_mat->parameter_sz)
                P = Faust::MatDense<FPP, Cpu>(cons_mat->parameter, constraints[i]->num_rows, constraints[i]->num_cols);
            else
                P = Faust::MatDense<FPP, Cpu>(cons_mat->parameter, cons_mat->parameter_sz/2, 2);
            tmp_cons = new Faust::ConstraintMat<FPP,Cpu>(static_cast<faust_constraint_name>(constraints[i]->name), P, constraints[i]->num_rows, constraints[i]->num_cols);
            cons.push_back(tmp_cons);
        }
        else
            handleError("FaustFact", "Invalid constraint.");
    }

    for(int i=0;i<num_facts-1;i++)
        fact_cons.push_back(cons[i]);

    for(int i=0;i<num_facts-1;i++)
        residuum_cons.push_back(cons[i+num_facts-1]);

    if(norm2_threshold == 0)
        norm2_threshold = FAUST_PRECISION; //faust_constant.h
    if(norm2_max_iter == 0)
        norm2_max_iter = FAUST_NORM2_MAX_ITER;

    Faust::TransformHelper<FPP,Cpu>* th, *th_times_lambda;
    try
    {
        vector<Faust::StoppingCriterion<Real<FPP>>> sc_ = {sc[0].num_its, sc[1].num_its};
        th = Faust::hierarchical(inMat, sc_, fact_cons, residuum_cons, inout_lambda[0], is_update_way_R2L,
                is_fact_side_left, use_csr,
                packing_RL,
                /* compute_2norm_on_array*/ false,
                norm2_threshold,
                norm2_max_iter);
        th_times_lambda = th->multiply(inout_lambda[0]);
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

    for (typename std::vector<const Faust::ConstraintGeneric*>::iterator it = cons.begin() ; it != cons.end(); ++it)
        delete *it;

    return core;

}
