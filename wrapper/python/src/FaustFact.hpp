#include <iostream>
#include "faust_MatDense.h"
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
void prepare_fact(const FPP* mat, const unsigned int num_rows, const unsigned int num_cols, const PyxParamsFact<FPP>* p,
        /* out args */vector<const Faust::ConstraintGeneric<FPP,Cpu>*>& cons, vector<Faust::MatDense<FPP,Cpu> > &initFacts)
{
    Faust::ConstraintGeneric<FPP, Cpu>* tmp_cons;
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
    PyxConstraintScalar<FPP>* cons_real;
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
            cons_real = static_cast<PyxConstraintScalar<FPP>*>(p->constraints[i]);
            if(p->is_verbose)
                cout << "constraint[" << i << "]->parameter: " << cons_real->parameter << endl;
            tmp_cons = new Faust::ConstraintFPP<FPP,Cpu>(static_cast<faust_constraint_name>(p->constraints[i]->name), cons_real->parameter, p->constraints[i]->num_rows, p->constraints[i]->num_cols);
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
    for(int i=0; i < p->num_facts;i++) {
        if(p->is_verbose)
        {
            cout << "init_facts[" << i << "] ele[0][0]: " << p->init_facts[i][0] << " size: " << p->init_fact_sizes[i*2] << "x"<< p->init_fact_sizes[i*2+1] << endl; 
        }
        initFacts.push_back(Faust::MatDense<FPP, Cpu>(p->init_facts[i], p->init_fact_sizes[i*2], p->init_fact_sizes[i*2+1]));
    }


}

template<typename FPP>
FaustCoreCpp<FPP>* fact_palm4MSA(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsFactPalm4MSA<FPP>* p)
{
    FaustCoreCpp<FPP>* core;
    Faust::MatDense<FPP,Cpu> inMat(mat, num_rows, num_cols);
    vector<const Faust::ConstraintGeneric<FPP,Cpu>*> cons;
    vector<Faust::MatDense<FPP,Cpu> > initFacts;
    if(p->is_verbose) {
        cout << "stop_crit.is_criterion_error: " << p->stop_crit.is_criterion_error << endl;
        cout << "stop_crit.error_treshold: " << p->stop_crit.error_treshold << endl;
        cout << "stop_crit.num_its: " << p->stop_crit.num_its << endl;
        cout << "stop_crit.max_num_its: " << p->stop_crit.max_num_its << endl;
    }
    prepare_fact(mat, num_rows, num_cols, p, cons, initFacts);

    // set all constructor arguments because they could be at non-default
    // values
    Faust::StoppingCriterion<FPP> crit(p->stop_crit.num_its, p->stop_crit.is_criterion_error, p->stop_crit.error_treshold, p->stop_crit.max_num_its);

    Faust::ParamsPalm<FPP,Cpu> params(inMat, p->num_facts, cons, initFacts, crit, p->is_verbose, p->is_update_way_R2L, p->init_lambda, p->constant_step_size, p->step_size);

    if(p->is_verbose) params.Display();

    Faust::BlasHandle<Cpu> blasHandle;

    Faust::Palm4MSA<FPP,Cpu> palm(params,blasHandle,true);

    palm.compute_facts();


//    palm.get_facts(faust);
//
    const std::vector<Faust::MatDense<FPP,Cpu> >& full_facts = palm.get_facts();
    Faust::Transform<FPP, Cpu> faust(full_facts, true);
    if(p->is_verbose) faust.Display();

    Faust::TransformHelper<FPP, Cpu> th(faust);
    core = new FaustCoreCpp<FPP>(th);

    blasHandle.Destroy();

    for (typename std::vector<const Faust::ConstraintGeneric<FPP,Cpu>*>::iterator it = cons.begin() ; it != cons.end(); ++it)
        delete *it;


    return core;

}

template<typename FPP>
FaustCoreCpp<FPP>* fact_hierarchical(FPP* mat, unsigned int num_rows, unsigned int num_cols, PyxParamsHierarchicalFact<FPP>* p)
{
    FaustCoreCpp<FPP>* core;
    Faust::MatDense<FPP,Cpu> inMat(mat, num_rows, num_cols);
    vector<const Faust::ConstraintGeneric<FPP,Cpu>*> cons;
    vector<std::vector<const Faust::ConstraintGeneric<FPP,Cpu>*>> cons_;
    vector<const Faust::ConstraintGeneric<FPP,Cpu>*> fact_cons;
    vector<const Faust::ConstraintGeneric<FPP,Cpu>*> residuum_cons;
    vector<Faust::MatDense<FPP,Cpu> > initFacts;
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
    prepare_fact(mat, num_rows, num_cols, p, cons, initFacts);

    // set all constructor arguments because they could be at non-default
    // values
    Faust::StoppingCriterion<FPP> crit0(p->stop_crits[0].num_its, p->stop_crits[0].is_criterion_error, p->stop_crits[0].error_treshold, p->stop_crits[0].max_num_its); //2 facts
    Faust::StoppingCriterion<FPP> crit1(p->stop_crits[1].num_its, p->stop_crits[1].is_criterion_error, p->stop_crits[1].error_treshold, p->stop_crits[1].max_num_its); //global

    for(int i=0;i<p->num_facts-1;i++)
        fact_cons.push_back(cons[i]);

    for(int i=0;i<p->num_facts-1;i++)
        residuum_cons.push_back(cons[i+p->num_facts-1]);

    cons_.push_back(fact_cons);
    cons_.push_back(residuum_cons);

    Faust::Params<FPP,Cpu> params(p->num_rows, p->num_cols, p->num_facts, cons_, initFacts, crit0, crit1, p->is_verbose, p->is_update_way_R2L, p->is_fact_side_left, p->init_lambda, p->constant_step_size, p->step_size);
    
    if(p->is_verbose) params.Display();


    Faust::BlasHandle<Cpu> blasHandle;
    Faust::SpBlasHandle<Cpu> spblasHandle;

    Faust::HierarchicalFact<FPP,Cpu> hierFact(inMat, params,blasHandle,spblasHandle);

    hierFact.compute_facts();

    vector<Faust::MatSparse<FPP,Cpu> > facts;
    hierFact.get_facts(facts);
    FPP lambda = hierFact.get_lambda();
    (facts[0]) *= lambda;
    // transform the sparse matrix into generic one
    std::vector<Faust::MatGeneric<FPP,Cpu> *> list_fact_generic;
    list_fact_generic.resize(facts.size());
    for (int i=0;i<list_fact_generic.size();i++)
        list_fact_generic[i]=facts[i].Clone();

    Faust::Transform<FPP,Cpu> faust(list_fact_generic);

    for (int i=0;i<list_fact_generic.size();i++)
        delete list_fact_generic[i];

    if(p->is_verbose) faust.Display();

    Faust::TransformHelper<FPP, Cpu> th(faust);
    core = new FaustCoreCpp<FPP>(th);

    blasHandle.Destroy();
    spblasHandle.Destroy();

    for (typename std::vector<const Faust::ConstraintGeneric<FPP,Cpu>*>::iterator it = cons.begin() ; it != cons.end(); ++it)
        delete *it;

    return core;

}
