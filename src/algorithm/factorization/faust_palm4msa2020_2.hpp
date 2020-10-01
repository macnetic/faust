template <typename FPP, FDevice DEVICE>
void Faust::palm4msa2(const Faust::MatDense<FPP,DEVICE>& A,
		std::vector<Faust::ConstraintGeneric*> & constraints,
		Faust::TransformHelper<FPP,DEVICE>& S,
		Real<FPP>& lambda, //TODO: FPP lambda ? is useful to have a complex lamdba ?
		//const unsigned int nites,
		const StoppingCriterion<Real<FPP>>& sc,
		const bool is_update_way_R2L,
		const bool use_csr,
		const bool packing_RL,
		const bool compute_2norm_on_array,
		const Real<FPP> norm2_threshold,
		const unsigned int norm2_max_iter,
		const bool constant_step_size, const Real<FPP> step_size,
		const bool on_gpu /*=false*/)
{
	double norm1, norm2;
//	std::cout << "palm4msa2 "<< std::endl;
//	std::cout << "on_gpu: " << on_gpu << std::endl;
	if(constraints.size() == 0)
		throw out_of_range("No constraint passed to palm4msa.");
	const Real<FPP> lipschitz_multiplicator = 1.001;
	Faust::MatGeneric<FPP,DEVICE>* cur_fac;
	Faust::MatSparse<FPP,DEVICE>* scur_fac;
	Faust::MatDense<FPP,DEVICE>* dcur_fac;
	unsigned int nfacts = constraints.size();
	Real<FPP> error = -1; //negative error is ignored
	std::vector<std::pair<faust_unsigned_int,faust_unsigned_int>> dims;
	int norm2_flag; // return val
	for(auto c: constraints)
		dims.push_back(make_pair(c->get_rows(), c->get_cols()));
	//TODO: make it possible to receive a MatSparse A
	Faust::MatDense<FPP,DEVICE> A_H = A;
	A_H.adjoint();
	if(S.size() != nfacts)
		fill_of_eyes(S, nfacts, use_csr, dims, on_gpu);
	else if(on_gpu)
		S.enable_gpu_meth_for_mul();
	int i = 0, f_id;
	std::function<void()> init_ite, next_fid;
	std::function<bool()> updating_facs;
	std::function<bool()> is_last_fac_updated;
	// packed Fausts corresponding to each factor
	std::vector<TransformHelper<FPP,DEVICE>*> pL, pR;
	pL.resize(nfacts);// pL[i] is the Faust for all factors to the left of the factor *(S.begin()+i)
	pR.resize(nfacts);// pR[i] the same for the right of S_i
	for(int i=0;i<nfacts;i++)
	{
		pL[i] = new TransformHelper<FPP,DEVICE>();
		pR[i] = new TransformHelper<FPP,DEVICE>();
	}
	if(is_update_way_R2L)
	{
		init_ite = [&f_id, &nfacts, &pL, &S, &packing_RL, &on_gpu]()
		{
			//pre-compute left products for each Si
			if(pL[0] != nullptr) delete pL[0];
			pL[0] = new TransformHelper<FPP,DEVICE>(); // empty faust // no factors to the left of *(S.begin())
			for(int i=1;i < nfacts; i++)
			{
				auto vec_Si_minus_1 = { *(S.begin()+i-1) };
				if(pL[i] != nullptr) delete pL[i]; //TODO: maybe to replace by a TransformHelper stored in the stack to avoid deleting each time
				pL[i] = new TransformHelper<FPP,DEVICE>(*pL[i-1], vec_Si_minus_1);
				// if the ctor args are GPU-enabled so is pL[i]
//				if(on_gpu) assert(10 == pL[i]->get_mul_order_opt_mode());
				if(packing_RL) pL[i]->pack_factors();
			}
			// all pL[i] Fausts are composed at most of one factor matrix
			f_id = nfacts-1;
		};
		next_fid = [&f_id, &pR, &S, &packing_RL, &on_gpu]()
		{
			if(f_id > 0)
			{
				if(pR[f_id-1] != nullptr)
					delete pR[f_id-1];
				auto vec_Sj = { *(S.begin()+f_id) };
				pR[f_id-1] = new Faust::TransformHelper<FPP,DEVICE>(vec_Sj, *pR[f_id]);
//				if(on_gpu) assert(10 == pR[f_id-1]->get_mul_order_opt_mode());
				if(packing_RL) pR[f_id-1]->pack_factors();
			}
			f_id--;
		};
		is_last_fac_updated = [&f_id]() {return f_id == 0;};
		updating_facs = [&f_id]() {return f_id >= 0;};
	}
	else
	{
		init_ite = [&f_id, &pR, &S, &packing_RL, &on_gpu]()
		{
			if(pR[S.size()-1] != nullptr) delete pR[S.size()-1];
			pR[S.size()-1] = new TransformHelper<FPP,DEVICE>(); // empty faust // no factors to the right of *(S.begin()+S.size()]-1)
			for(int i=S.size()-2;i >= 0; i--)
			{
				auto vec_Si_plus_1 = { *(S.begin()+i+1) };
				if(pR[i] != nullptr) delete pR[i];
				pR[i] = new TransformHelper<FPP,DEVICE>(vec_Si_plus_1, *pR[i+1]);
//				if(on_gpu) assert(10 == pR[i]->get_mul_order_opt_mode());
				if(packing_RL) pR[i]->pack_factors();
			}
			f_id = 0;
		};
		next_fid = [&f_id, &S, &pL, &nfacts, packing_RL, &on_gpu]()
		{
			if(f_id < nfacts-1)
			{
				if(pL[f_id+1] != nullptr)
					delete pL[f_id+1];
				auto vec_Sj = { *(S.begin()+f_id) };
				pL[f_id+1] = new Faust::TransformHelper<FPP,DEVICE>(*pL[f_id], vec_Sj);
//				if(on_gpu) assert(10 == pL[f_id+1]->get_mul_order_opt_mode());
				if(packing_RL) pL[f_id+1]->pack_factors();
			}
			f_id++;
		};
		updating_facs = [&f_id, &nfacts]() {return f_id < nfacts;};
		is_last_fac_updated = [&f_id, &nfacts]() {return f_id == nfacts-1;};
	}
	Faust::MatDense<FPP,DEVICE> D, tmp;
	Faust::MatSparse<FPP,DEVICE> spD;
	Faust::MatDense<FPP,DEVICE> * LorR;
	Faust::MatDense<FPP,DEVICE> _LorR;
	Real<FPP> c = 1/step_size;
	while(sc.do_continue(i, error))
	{
		//		std::cout << "nfacts:" << nfacts << std::endl;

		init_ite();
		while(updating_facs())
		{
//			std::cout << "f_id: " << f_id << std::endl;
			cur_fac = S.get_gen_fact_nonconst(f_id);
			Real<FPP> nR=1,nL=1;
			if(! constant_step_size)
			{
				if(pR[f_id]->size() > 0)
					nR = pR[f_id]->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
				if(pL[f_id]->size() > 0)
					nL = pL[f_id]->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
				c = lipschitz_multiplicator*lambda*lambda*nR*nR*nL*nL;
			}
			auto S_j_vec = {*(S.begin()+f_id)};
			Faust::TransformHelper<FPP, DEVICE> _LSR(*pL[f_id], S_j_vec, *pR[f_id]);
//			tmp = _LSR.get_product();
			_LSR.get_product(tmp);
			tmp *= FPP(lambda);
			tmp -= A;
			if(sc.isCriterionErr())
				error = tmp.norm();
			FPP alpha_R = 1, alpha_L = 1, beta_R = 0, beta_L = 0; //decl in parent scope
			if(S.is_fact_sparse(f_id))
			{
				scur_fac = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(cur_fac);
				D = *scur_fac;
				dcur_fac = nullptr;
			}
			else
			{
				dcur_fac = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(cur_fac); // TOFIX: possible it is not Dense... but MatDiag (sanity check on function start)
				D = *dcur_fac;
				scur_fac = nullptr;
			}
			if(pR[f_id]->size() > 0)
			{
				if(packing_RL)
					LorR = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(pR[f_id]->get_gen_fact_nonconst(0)); //normally pR[f_id] is packed (hence reduced to a single MatDense)
				else
				{
					_LorR = pR[f_id]->get_product();
					LorR = &_LorR;
				}
				if(pL[f_id]->size() == 0)
				{ //no L factor for factor f_id
					alpha_R = - lambda/c;
					beta_R = 1;
					gemm(tmp, *LorR, D, alpha_R, beta_R, 'N', 'H');
				}
				else
					gemm(tmp, *LorR, tmp, alpha_R, beta_R, 'N', 'H');
			}
			if(pL[f_id]->size() > 0)
			{
				if(packing_RL)
					LorR = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(pL[f_id]->get_gen_fact_nonconst(0));
				else
				{
					_LorR = pL[f_id]->get_product();
					LorR = &_LorR;
				}
				alpha_L = -lambda/c;
				beta_L = 1;
				gemm(*LorR, tmp, D, alpha_L, beta_L, 'H', 'N');
			}
			// really update now
			constraints[f_id]->project<FPP,DEVICE,Real<FPP>>(D);
			if(use_csr && dcur_fac != nullptr || !use_csr && scur_fac != nullptr)
				throw std::runtime_error("Current factor is inconsistent with use_csr.");
			if(use_csr)
			{
				spD = D;
				S.update(spD, f_id); // update is at higher level than a simple assignment
			}
			else
			{
				S.update(D, f_id);
			}
			next_fid(); //f_id updated to iteration factor index (pL or pR too)
		}
		//update lambda
		//TODO: variable decl in parent scope
		Faust::MatDense<FPP,DEVICE> A_H_S = S.multiply(A_H);
		//		auto last_Sfac_vec = { *(S.begin()+nfacts-1), dynamic_cast<Faust::MatGeneric<FPP,DEVICE>*>(&A_H)};
		//		Faust::TransformHelper<FPP,DEVICE> A_H_S_(*pL[nfacts-1], last_Sfac_vec);
		//		A_H_S_.disable_dtor();
		//		Faust::MatDense<FPP,DEVICE> A_H_S = A_H_S_.get_product();
		Real<FPP> trr = std::real(A_H_S.trace());
		Real<FPP> n = S.normFro();
		lambda = trr/(n*n);
		//		std::cout << "debug lambda: " << lambda << std::endl;
		i++;
	}
	S.update_total_nnz();
}

