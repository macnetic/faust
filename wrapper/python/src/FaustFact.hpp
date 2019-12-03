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
            return true;
        default:
            return false;
    }
}

template<typename FPP>
void prox_mat(unsigned int cons_type, FPP* cons_param, FPP* mat_in, unsigned long num_rows, unsigned long num_cols, FPP* mat_out)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
    switch(static_cast<faust_constraint_name>(cons_type))
	{
		case CONSTRAINT_NAME_CONST: /**< Matrix equal to A ; MAT */
			// nothing to do, same mat returned
			break;
		case CONSTRAINT_NAME_BLKDIAG:
			//not impl. yet in cpp core
			break;
		case CONSTRAINT_NAME_SUPP: /**< Matrix which support is equal to A ; MAT ; (frobenius norm 1)*/
			Faust::prox_supp(fmat, Faust::MatDense<FPP,Cpu>(cons_param, num_rows, num_cols));
			break;
		default:
			throw invalid_argument("PyxConstraintMat::project() inconsistent constraint name");
	}
    memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
}

template<typename FPP>
void prox_int(unsigned int cons_type, unsigned long cons_param, FPP* mat_in, unsigned long num_rows,
        unsigned long num_cols, FPP* mat_out, const bool normalized /* default to true */)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
    switch(static_cast<faust_constraint_name>(cons_type))
    {
        case CONSTRAINT_NAME_SPCOL: /*!< fixed number of non zero elements per column INT (frobenius norm 1) */
            if(normalized)
                Faust::prox_spcol(fmat, (faust_unsigned_int) cons_param);
            else
                Faust::prox_spcol_normfree(fmat, (faust_unsigned_int) cons_param);
            break;
        case CONSTRAINT_NAME_SPLIN: /*!< fixed number of non zero elements per line INT (frobenius norm 1) */
            if(normalized)
                Faust::prox_splin(fmat, (faust_unsigned_int) cons_param);
            else
                Faust::prox_splin_normfree(fmat, (faust_unsigned_int) cons_param);
            break;
        case CONSTRAINT_NAME_SPLINCOL:
            if(normalized)
                Faust::prox_splincol(fmat, (faust_unsigned_int) cons_param);
            else
                Faust::prox_splincol_normfree(fmat, (faust_unsigned_int) cons_param);
            break;
        case CONSTRAINT_NAME_SP_POS:/**< fixed number of non zeros coefficients: INT (frobenius norm 1) */
            if(normalized)
                Faust::prox_sp_pos(fmat, (faust_unsigned_int) cons_param);
            else
                Faust::prox_sp_pos_normfree(fmat, (faust_unsigned_int) cons_param);
            break;
        case CONSTRAINT_NAME_SP:
            if(normalized)
                Faust::prox_sp(fmat, (faust_unsigned_int) cons_param);
            else
                Faust::prox_sp_normfree(fmat, (faust_unsigned_int) cons_param);
            break;
        default:
            throw invalid_argument("PyxConstraintInt::project() inconsistent constraint name");
    }
    memcpy(mat_out, fmat.getData(), sizeof(FPP)*num_rows*num_cols);
}

template<typename FPP, typename FPP2>
void prox_real(unsigned int cons_type, FPP2 cons_param, FPP* mat_in, unsigned long num_rows, unsigned long num_cols, FPP* mat_out)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
    switch(static_cast<faust_constraint_name>(cons_type))
	{
		case CONSTRAINT_NAME_NORMLIN:/**< 2nd norm of the lines of matrix A ; REAL  */
			Faust::prox_normlin(fmat, cons_param);
			break;
		case CONSTRAINT_NAME_NORMCOL:/*!< 2nd norm of the columns of A REAL */
			Faust::prox_normcol(fmat, cons_param);
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
            tmp_cons = new Faust::ConstraintMat<FPP,Cpu>(static_cast<faust_constraint_name>(p->constraints[i]->name), Faust::MatDense<FPP, Cpu>(cons_mat->parameter, p->constraints[i]->num_rows, p->constraints[i]->num_cols), p->constraints[i]->num_rows, p->constraints[i]->num_cols);
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
        params = new ParamsPalmFGFT<FPP,Cpu,FPP2>(inMat, p->num_facts, cons, initFacts, p_fft->init_D, crit, p->is_verbose, p->is_update_way_R2L, p->init_lambda, p->step_size);

        palm = new Palm4MSAFGFT<FPP,Cpu,FPP2>(*static_cast<ParamsPalmFGFT<FPP,Cpu,FPP2>*>(params),blasHandle,true);
    }
    else {
        params = new ParamsPalm<FPP,Cpu,FPP2>(inMat, p->num_facts, cons, initFacts, crit, p->is_verbose, p->is_update_way_R2L, p->init_lambda, p->constant_step_size, p->step_size);
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

    FPP lambda = palm->get_lambda();

    std::vector<Faust::MatDense<FPP,Cpu> > facts;
    std::vector<Faust::MatGeneric<FPP, Cpu>*> sp_facts;
    facts=palm->get_facts();

    for(typename std::vector<Faust::MatDense<FPP, Cpu>>::iterator it = facts.begin(); it != facts.end(); it++)
    {
        Faust::MatSparse<FPP, Cpu> * M = new Faust::MatSparse<FPP, Cpu>(*it);
        sp_facts.push_back(M);
    }

    Faust::TransformHelper<FPP, Cpu> *th = new Faust::TransformHelper<FPP,Cpu>(sp_facts, lambda, false);

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
        *out_buf = lambda;
    else
    {
        // retrieve D matrix from Palm4MSAFFT
        // out buffer must have been allocated from outside
        dynamic_cast<Palm4MSAFGFT<FPP,Cpu,FPP2>*>(palm)->get_D(out_buf+1);
        // add lambda at the first position
        out_buf[0] = lambda;
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
        params = new Faust::ParamsFGFT<FPP,Cpu,FPP2>(p->num_rows, p->num_cols, p->num_facts, cons_, initFacts_deft, p_fft->init_D, crit0, crit1, p->is_verbose, p->is_update_way_R2L, p->is_fact_side_left, p->init_lambda, p->constant_step_size, p->step_size);
        hierFact = new HierarchicalFactFGFT<FPP,Cpu,FPP2>(inMat, *inMat2, *(static_cast<ParamsFGFT<FPP,Cpu,FPP2>*>(params)), blasHandle, spblasHandle);
    }
    else
    {
        params = new Params<FPP,Cpu,FPP2>(p->num_rows, p->num_cols, p->num_facts, cons_, initFacts_deft, crit0, crit1, p->is_verbose, p->is_update_way_R2L, p->is_fact_side_left, p->init_lambda, p->constant_step_size, p->step_size);
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
    FPP lambda = hierFact->get_lambda();
    *out_buf =  lambda;
    (facts[0]) *= lambda;
    // transform the sparse matrix into matrix pointers
    std::vector<Faust::MatGeneric<FPP,Cpu> *> list_fact_generic;
    list_fact_generic.resize(facts.size());
    for (int i=0;i<list_fact_generic.size();i++)
        list_fact_generic[i]=facts[i].Clone();

    //don't delete list_fact_generic because they are directly used in
    //TransformHelper (without copy)

    Faust::TransformHelper<FPP, Cpu>* th = new Faust::TransformHelper<FPP,Cpu>(list_fact_generic, 1.0, false, false);

    if(p->is_verbose) th->display();
    core = new FaustCoreCpp<FPP>(th);

    blasHandle.Destroy();
    spblasHandle.Destroy();

    for (typename std::vector<const Faust::ConstraintGeneric*>::iterator it = cons.begin() ; it != cons.end(); ++it)
        delete *it;


    if(p_fft == nullptr)
        // palm4MSA basis case, just get lambda
        *out_buf = lambda;
    else
    {
        // retrieve D matrix from HierarchicalFactFFT
        // out buffer must have been allocated from outside
        dynamic_cast<HierarchicalFactFGFT<FPP,Cpu,FPP2>*>(hierFact)->get_D(out_buf+1);
        // add lambda at the first position
        out_buf[0] = lambda;
        delete inMat2;
    }

    delete params;
    delete hierFact;
    return core;

}
