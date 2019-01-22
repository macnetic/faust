#include <iostream>
#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include "faust_Params.h"
#include "faust_ParamsPalm.h"
#include "faust_StoppingCriterion.h"
#include "faust_Palm4MSA.h"
#include "faust_HierarchicalFact.h"
#include "faust_BlasHandle.h"
#include "faust_ConstraintGeneric.h"

using namespace std;

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
        unsigned long num_cols, FPP* mat_out)
{
    Faust::MatDense<FPP, Cpu> fmat(mat_in, num_rows, num_cols);
    switch(static_cast<faust_constraint_name>(cons_type))
    {
        case CONSTRAINT_NAME_SPCOL: /*!< fixed number of non zero elements per column INT (frobenius norm 1) */
            Faust::prox_spcol(fmat, (faust_unsigned_int) cons_param);
            break;
        case CONSTRAINT_NAME_SPLIN: /*!< fixed number of non zero elements per line INT (frobenius norm 1) */
            Faust::prox_splin(fmat, (faust_unsigned_int) cons_param);
            break;
        case CONSTRAINT_NAME_SPLINCOL:
            Faust::prox_splincol(fmat, (faust_unsigned_int) cons_param);
            break;
        case CONSTRAINT_NAME_SP_POS:/**< fixed number of non zeros coefficients: INT (frobenius norm 1) */
            Faust::prox_sp_pos(fmat, (faust_unsigned_int) cons_param);
            break;
        case CONSTRAINT_NAME_SP:
            Faust::prox_sp(fmat, (faust_unsigned_int) cons_param);            break;
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
FaustCoreCpp<FPP>* fact_palm4MSA(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsFactPalm4MSA<FPP,FPP2>* p, FPP* out_lambda)
{
    FaustCoreCpp<FPP>* core;
    Faust::MatDense<FPP,Cpu> inMat(mat, num_rows, num_cols);
    vector<const Faust::ConstraintGeneric*> cons;
    vector<Faust::MatDense<FPP,Cpu> > initFacts;
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

    Faust::ParamsPalm<FPP,Cpu,FPP2> params(inMat, p->num_facts, cons, initFacts, crit, p->is_verbose, p->is_update_way_R2L, p->init_lambda, p->constant_step_size, p->step_size);

    if(p->is_verbose) params.Display();

    Faust::BlasHandle<Cpu> blasHandle;

    Faust::Palm4MSA<FPP,Cpu, FPP2> palm(params,blasHandle,true);

    palm.compute_facts();

    FPP lambda = palm.get_lambda();

    std::vector<Faust::MatDense<FPP,Cpu> > facts;
    std::vector<Faust::MatGeneric<FPP, Cpu>*> sp_facts;
    facts=palm.get_facts();

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

    *out_lambda = lambda;

    return core;

}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_hierarchical(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsHierarchicalFact<FPP, FPP2>* p, FPP* out_lambda)
{
    FaustCoreCpp<FPP>* core;
    Faust::MatDense<FPP,Cpu> inMat(mat, num_rows, num_cols);
    vector<const Faust::ConstraintGeneric*> cons;
    vector<std::vector<const Faust::ConstraintGeneric*>> cons_;
    vector<const Faust::ConstraintGeneric*> fact_cons;
    vector<const Faust::ConstraintGeneric*> residuum_cons;
    vector<Faust::MatDense<FPP,Cpu> > initFacts_deft;
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

    Faust::Params<FPP,Cpu,FPP2> params(p->num_rows, p->num_cols, p->num_facts, cons_, initFacts_deft, crit0, crit1, p->is_verbose, p->is_update_way_R2L, p->is_fact_side_left, p->init_lambda, p->constant_step_size, p->step_size);
    
    if(p->is_verbose) params.Display();


    Faust::BlasHandle<Cpu> blasHandle;
    Faust::SpBlasHandle<Cpu> spblasHandle;

    Faust::HierarchicalFact<FPP,Cpu,FPP2> hierFact(inMat, params,blasHandle,spblasHandle);

    hierFact.compute_facts();

    vector<Faust::MatSparse<FPP,Cpu> > facts;
    hierFact.get_facts(facts);
    FPP lambda = hierFact.get_lambda();
    *out_lambda =  lambda;
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

    return core;

}
