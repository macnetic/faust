template<typename FPP, Device DEVICE>
Faust::TransformHelper<FPP,DEVICE>* Faust::hierarchical(const Faust::MatDense<FPP,DEVICE>& A,
        const int nites,
        std::vector<const Faust::ConstraintGeneric*> & fac_constraints,
        std::vector<const Faust::ConstraintGeneric*> & res_constraints,
        Real<FPP>& lambda,
        const bool is_update_way_R2L, const bool is_fact_side_left,
        const bool use_csr,
        const bool compute_2norm_on_array,
        const Real<FPP> norm2_threshold,
        const unsigned int norm2_max_iter)
{
    auto S = new Faust::TransformHelper<FPP,DEVICE>(); // A is copied
    S->push_back(&A);
    Faust::MatGeneric<FPP,DEVICE> *Si;
    const Faust::ConstraintGeneric *fac_cons, *res_cons;
    Faust::MatGeneric<FPP,DEVICE> *zero_mat, *id_mat;
    Faust::MatDense<FPP,DEVICE> * tmp_dense;
    Faust::MatSparse<FPP,DEVICE> * tmp_sparse;
    std::vector<Faust::MatGeneric<FPP,DEVICE>*> Si_vec;
    std::vector<Faust::ConstraintGeneric*> Si_cons;
    Real<FPP> lambda_ = 1;
    Real<FPP> glo_lambda = 1;
    for(int i=0;i < fac_constraints.size();i++)
    {
        cout << "Faust::hierarchical: " << i+1 << endl;
        if(is_fact_side_left)
        {
            Si = S->get_gen_fact_nonconst(0);

        }
        else
            Si = S->get_gen_fact_nonconst(i);
        fac_cons = fac_constraints[i];
        res_cons = res_constraints[i];
        // init factors for the local optimization (factorization of Si)
        //TODO: refactor into a separate function init_zero_id
        if(use_csr)
        {
            tmp_sparse = new Faust::MatSparse<FPP,DEVICE>(fac_constraints[i]->get_rows(), fac_constraints[i]->get_cols());
            tmp_sparse->setZeros();
            zero_mat = tmp_sparse;
            tmp_sparse = new Faust::MatSparse<FPP,DEVICE>(res_constraints[i]->get_rows(), res_constraints[i]->get_cols());
            tmp_sparse->setEyes();
            id_mat = tmp_sparse;
        }
        else
        {

            tmp_dense = new Faust::MatDense<FPP,DEVICE>(fac_constraints[i]->get_rows(), fac_constraints[i]->get_cols());
            tmp_dense->setZeros();
            zero_mat = tmp_dense;
            tmp_dense = new Faust::MatDense<FPP,DEVICE>(res_constraints[i]->get_rows(), res_constraints[i]->get_cols());
            tmp_dense->setEyes();
            id_mat = tmp_dense;
        }
        if(is_update_way_R2L)
            Si_vec = {id_mat, zero_mat};
        else
            Si_vec = {zero_mat, id_mat};
        Si_cons = { const_cast<Faust::ConstraintGeneric*>(fac_cons), const_cast<Faust::ConstraintGeneric*>(res_cons) };
        Faust::TransformHelper<FPP,DEVICE> Si_th(Si_vec, 1.0, false, false, true);
        lambda_ = 1;
        tmp_dense = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(Si);
        if(tmp_dense == nullptr)
        {
            tmp_sparse = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(Si);
            tmp_dense = new MatDense<FPP,Cpu>(*tmp_sparse);
        }
        else tmp_sparse = nullptr;
        Faust::palm4msa2(*tmp_dense, Si_cons, Si_th, lambda_, nites, is_update_way_R2L , use_csr, compute_2norm_on_array, norm2_threshold, norm2_max_iter);
        if(tmp_sparse != nullptr)
            // the Si factor has been converted into a MatDense in the memory
            // storage
            // delete it // TODO: palm4msa2 should handle this on its own
            delete tmp_dense;
//        cout << "Local opt. result:" << endl;
//        Si_th.display();
        //global optimization
        glo_lambda *= lambda_;
        if(is_fact_side_left)
        {
            S->pop_front();
            S->push_first(*(Si_th.begin()+1), false, false);
            S->push_first(*(Si_th.begin()), false, false);
        }
        else
        {
            S->pop_back();
            S->push_back(*(Si_th.begin()), false, false);
            S->push_back(*(Si_th.begin()+1), false, false);
        }
//        cout << "update global faust:" << endl;
//        S->display();
        //TODO: verify if the constraints order doesn't depend on
        //is_fact_side_left
        std::vector<ConstraintGeneric*> glo_cons;
        for(auto ite_cons=fac_constraints.begin(); ite_cons != fac_constraints.begin()+i+1;ite_cons++)
            glo_cons.push_back(const_cast<Faust::ConstraintGeneric*>(*ite_cons));
        glo_cons.push_back(const_cast<Faust::ConstraintGeneric*>(res_constraints[i]));
        // TODO: arguments with default values for norm threshold, norm num
        // iters
        // global optimization
//        cout << "S before global opt.:" << endl;
//        S->display();
        Faust::palm4msa2(A, glo_cons, *S, glo_lambda, nites, is_update_way_R2L, use_csr, compute_2norm_on_array, norm2_threshold, norm2_max_iter);
//        cout << "S after global opt.:" << endl;
//        S->display();
    }
    lambda = glo_lambda;
    return S;
}
