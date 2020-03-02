template <typename FPP, Device DEVICE>
void Faust::palm4msa(const Faust::MatDense<FPP,DEVICE>& A,
		std::vector<Faust::ConstraintGeneric*> & constraints,
		Faust::TransformHelper<FPP,DEVICE>& S,
		FPP& lambda,
		const unsigned int nites,
		const bool is_update_way_R2L,
		const bool use_csr,
		const bool compute_2norm_on_array,
		const Real<FPP> norm2_threshold,
		const unsigned int norm2_max_iter)
{
	if(constraints.size() == 0)
		throw out_of_range("No constraint passed to palm4msa.");
	const Real<FPP> lipschitz_multiplicator = 1.001;
	Faust::MatGeneric<FPP,DEVICE>* cur_fac;
	Faust::MatSparse<FPP,DEVICE>* scur_fac;
	Faust::MatDense<FPP,DEVICE>* dcur_fac;
	const unsigned int nfacts = constraints.size();
	std::vector<std::pair<faust_unsigned_int,faust_unsigned_int>> dims;
	int norm2_flag; // return val
	for(auto c: constraints)
		dims.push_back(make_pair(c->get_rows(), c->get_cols()));
	Faust::MatDense<FPP,DEVICE> A_H = A;
	A_H.conjugate(false);
	A_H.transpose();
	if(S.size() != nfacts)
	{
//		S = Faust::TransformHelper<FPP,DEVICE>();
		//TODO: refactor the id factor gen. into TransformHelper
		for(auto fdims : dims)
		{
			// init all facts as identity matrices
			// with proper dimensions
			Faust::MatGeneric<FPP,DEVICE>* fact;
			if(use_csr)
			{
				auto sfact = new Faust::MatSparse<FPP,DEVICE>(fdims.first, fdims.second);
				sfact->setEyes();
				fact = sfact;
			}
			else
			{
				auto dfact = new Faust::MatDense<FPP,DEVICE>(fdims.first, fdims.second);
				dfact->setEyes();
				fact = dfact;
			}
			S.push_back(fact); //TODO: copying=false
		}
	}
	int i = 0, f_id;
	std::function<void()> init_fid, next_fid;
	std::function<bool()> updating_facs;
	if(is_update_way_R2L)
	{
		init_fid = [&f_id, &nfacts]() {f_id = nfacts-1;};
		next_fid = [&f_id]() {f_id--;};
		updating_facs = [&f_id]() {return f_id >= 0;};
	}
	else
	{
		init_fid = [&f_id]() {f_id = 0;};
		next_fid = [&f_id]() {f_id++;};
		updating_facs = [&f_id, &nfacts]() {return f_id < nfacts;};
	}

	Faust::MatDense<FPP,Cpu> D, tmp;
	Faust::TransformHelper<FPP, Cpu> LSR;
	// lambda exp to update fact when its id is 0
	auto update_1stfac = [&A, &D, &tmp, &LSR, &scur_fac, &dcur_fac, &f_id, &S, &lipschitz_multiplicator, &lambda, &norm2_threshold, &norm2_flag, &norm2_max_iter](Faust::MatGeneric<FPP, DEVICE> *cur_fac)
	{
		auto R = S.right(f_id+1);
		Real<FPP> nR = R->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
		Real<FPP> c = lipschitz_multiplicator*lambda*lambda*nR*nR;
//		std::cout << "c=" << c << std::endl;
		//TODO: check if c is nan
		Faust::TransformHelper<FPP, Cpu> L;
		scur_fac = nullptr, dcur_fac = nullptr;
		if(S.is_fact_sparse(f_id))
		{
			scur_fac = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(cur_fac);
			D = *scur_fac;
			Faust::TransformHelper<FPP, Cpu> _L({scur_fac}, 1.0, false, false);
			L = _L;
		}
		else
		{
			dcur_fac = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(cur_fac); // TOFIX: possible it is not Dense... but MatDiag (sanity check on function start)
			D = *dcur_fac;
			Faust::TransformHelper<FPP, Cpu> _L({dcur_fac}, 1.0, false, false);
			L = _L;
		}
		Faust::TransformHelper<FPP, Cpu> _LSR(&L, R); // L contains cur_fac
//		LSR = _LSR;
//		_LSR.multiply(lambda);
//		tmp = _LSR.get_product();
		_LSR.get_product(tmp);
		tmp *= lambda;
		tmp -= A;
		//TODO: do something to lighten the double transpose conjugate
		tmp.conjugate(false);
		tmp.transpose();
		tmp = R->multiply(tmp, /* H */ false, false);
		tmp.conjugate(false);
		tmp.transpose();
		tmp *= lambda/c;
		D -= tmp;
	};
	auto update_lastfac = [&A, &D, &tmp, &LSR, &scur_fac, &dcur_fac, &f_id, &S, &lipschitz_multiplicator, &lambda, &norm2_threshold, &norm2_flag, &norm2_max_iter](Faust::MatGeneric<FPP, DEVICE> *cur_fac)
	{
		auto L = S.left(f_id-1);
		// TODO: factorize with other lambda exp
		Real<FPP> nL = L->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
		Real<FPP> c = lipschitz_multiplicator*lambda*lambda*nL*nL;
//		std::cout << "c=" << c << std::endl;
		//TODO: check if c is nan
		Faust::TransformHelper<FPP, Cpu> R;
		scur_fac = nullptr, dcur_fac = nullptr;
		if(S.is_fact_sparse(f_id))
		{
			scur_fac = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(cur_fac);
			D = *scur_fac;
			Faust::TransformHelper<FPP, Cpu> _R({scur_fac}, 1.0, false, false);
			R = _R;
		}
		else
		{
			dcur_fac = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(cur_fac); // TOFIX: possible it is not Dense... but MatDiag (sanity check on function start)
			D = *dcur_fac;
			Faust::TransformHelper<FPP, Cpu> _R({dcur_fac}, 1.0, false, false);
			R = _R;
		}
		Faust::TransformHelper<FPP, Cpu> _LSR(L, &R); // L contains cur_fac
//		LSR = _LSR;
//		_LSR.multiply(lambda);
		tmp = _LSR.get_product();
		tmp *= lambda;
		tmp -= A;
		tmp = L->multiply(tmp, /* NO H */ true, true);
		tmp *= lambda/c;
		D -= tmp;
	};
	auto update_interfac = [&A, &D, &tmp, &LSR, &scur_fac, &dcur_fac, &f_id, &S, &lipschitz_multiplicator, &lambda, &norm2_threshold, &norm2_flag, &norm2_max_iter](Faust::MatGeneric<FPP, DEVICE> *cur_fac)
	{
		auto R = S.right(f_id+1);
		auto L = S.left(f_id-1);
		Real<FPP> nR = R->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
		Real<FPP> nL = L->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
		Real<FPP> c = lipschitz_multiplicator*lambda*lambda*nR*nR*nL*nL;
		//TODO: check if c is nan
		scur_fac = nullptr, dcur_fac = nullptr;
		if(S.is_fact_sparse(f_id))
		{
			scur_fac = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(cur_fac);
			D = *scur_fac;
		}
		else
		{
			dcur_fac = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(cur_fac); // TOFIX: possible it is not Dense... but MatDiag (sanity check on function start)
			D = *dcur_fac;
		}
		Faust::TransformHelper<FPP, Cpu> _LSR(L, R); // L contains cur_fac
//		LSR = _LSR;
//		_LSR.multiply(lambda);
		tmp = _LSR.get_product();
		tmp *= lambda;
		tmp -= A;
		//TODO: do something to lighten the double transpose conjugate
		tmp.conjugate(false);
		tmp.transpose();
		tmp = R->multiply(tmp, /* NO H */ false, false);
		tmp.conjugate(false);
		tmp.transpose();
		tmp = L->multiply(tmp, true, true);
		tmp *= lambda/c;
		D -= tmp;
	};
	while(i < nites)
	{
		//		std::cout << "nfacts:" << nfacts << std::endl;
		init_fid();
		while(updating_facs())
		{
			//			std::cout << "f_id: " << f_id << std::endl;
			cur_fac = S.get_gen_fact_nonconst(f_id);
			if(f_id == 0)
				update_1stfac(cur_fac);
			else if(f_id == nfacts-1)
				update_lastfac(cur_fac);
			else
				update_interfac(cur_fac);
			// really update now

			constraints[f_id]->project<FPP,DEVICE,Real<FPP>>(D);
			if(! use_csr && scur_fac == nullptr)
			{
				// not CSR and cur_fac DENSE
				*dcur_fac = D;
			}
			else if(use_csr && dcur_fac == nullptr)
			{
				// CSR and cur_fac SPARSE
				*scur_fac = D;
			}
			else
				throw std::runtime_error("Current factor is inconsistent with use_csr.");
			cur_fac->set_id(false);
			S.update_total_nnz();
			next_fid(); //f_id updated to iteration factor index
		}
		//update lambda
		//TODO: variable decl in parent scope
		Faust::MatDense<FPP,DEVICE> A_H_S = S.multiply(A_H);
		FPP tr = A_H_S.trace();
		Real<FPP> trr = std::real(tr);
		Real<FPP> n = S.normFro();
		lambda = trr/(n*n);
//		std::cout << "debug lambda: " << lambda << std::endl;
		i++;
	}
}

template <typename FPP, Device DEVICE>
void Faust::palm4msa2(const Faust::MatDense<FPP,DEVICE>& A,
		std::vector<Faust::ConstraintGeneric*> & constraints,
		Faust::TransformHelper<FPP,DEVICE>& S,
		Real<FPP>& lambda, //TODO: FPP lambda ? is useful to have a complex lamdba ?
		const unsigned int nites,
		const bool is_update_way_R2L,
		const bool use_csr,
		const bool packing_RL,
		const bool compute_2norm_on_array,
		const Real<FPP> norm2_threshold,
		const unsigned int norm2_max_iter)
{
	if(constraints.size() == 0)
		throw out_of_range("No constraint passed to palm4msa.");
	const Real<FPP> lipschitz_multiplicator = 1.001;
	Faust::MatGeneric<FPP,DEVICE>* cur_fac;
	Faust::MatSparse<FPP,DEVICE>* scur_fac;
	Faust::MatDense<FPP,DEVICE>* dcur_fac;
	const unsigned int nfacts = constraints.size();
	std::vector<std::pair<faust_unsigned_int,faust_unsigned_int>> dims;
	int norm2_flag; // return val
	for(auto c: constraints)
		dims.push_back(make_pair(c->get_rows(), c->get_cols()));
	//TODO: make it possible to have a MatSparse A
	Faust::MatDense<FPP,DEVICE> A_H = A;
	A_H.conjugate(false);
	A_H.transpose();
	if(S.size() != nfacts)
	{
//		S = Faust::TransformHelper<FPP,DEVICE>();
		//TODO: refactor the id factor gen. into TransformHelper
		for(auto fdims : dims)
		{
			// init all facts as identity matrices
			// with proper dimensions
			Faust::MatGeneric<FPP,DEVICE>* fact;
			if(use_csr)
			{
				auto sfact = new Faust::MatSparse<FPP,DEVICE>(fdims.first, fdims.second);
				sfact->setEyes();
				fact = sfact;
			}
			else
			{
				auto dfact = new Faust::MatDense<FPP,DEVICE>(fdims.first, fdims.second);
				dfact->setEyes();
				fact = dfact;
			}
			S.push_back(fact); //TODO: copying=false
		}
	}
	int i = 0, f_id;
	std::function<void()> init_ite, next_fid;
	std::function<bool()> updating_facs;
	std::function<bool()> is_last_fac_updated;
	// packed Fausts corresponding to each factor
	TransformHelper<FPP,Cpu>* pL[nfacts]; // pL[i] is the Faust for all factors to the left of the factor *(S.begin()+i)
	TransformHelper<FPP,Cpu>* pR[nfacts]; // pR[i] the same for the right of S_i
	for(int i=0;i<nfacts;i++)
	{
		pL[i] = new TransformHelper<FPP,Cpu>();
		pR[i] = new TransformHelper<FPP,Cpu>();
	}
	if(is_update_way_R2L)
	{
		init_ite = [&f_id, &nfacts, &pL, &S, &packing_RL]()
		{
			//pre-compute left products for each Si
			if(pL[0] != nullptr) delete pL[0];
			pL[0] = new TransformHelper<FPP,Cpu>(); // empty faust // no factors to the left of *(S.begin())
			for(int i=1;i < nfacts; i++)
			{
				auto vec_Si_minus_1 = { *(S.begin()+i-1) };
				if(pL[i] != nullptr) delete pL[i]; //TODO: maybe to replace by a TransformHelper stored in the stack to avoid deleting each time
				pL[i] = new TransformHelper<FPP,Cpu>(*pL[i-1], vec_Si_minus_1);
				if(packing_RL) pL[i]->pack_factors();
			}
			// all pL[i] Fausts are composed at most of one factor matrix
			f_id = nfacts-1;
		};
		next_fid = [&f_id, &pR, &S, &packing_RL]()
		{
			if(f_id > 0)
			{
				if(pR[f_id-1] != nullptr)
					delete pR[f_id-1];
				auto vec_Sj = { *(S.begin()+f_id) };
				pR[f_id-1] = new Faust::TransformHelper<FPP,Cpu>(vec_Sj, *pR[f_id]);
				if(packing_RL) pR[f_id-1]->pack_factors();
			}
			f_id--;
		};
		is_last_fac_updated = [&f_id]() {return f_id == 0;};
		updating_facs = [&f_id]() {return f_id >= 0;};
	}
	else
	{
		init_ite = [&f_id, &pR, &S, &packing_RL]()
		{
			if(pR[S.size()-1] != nullptr) delete pR[S.size()-1];
			pR[S.size()-1] = new TransformHelper<FPP,Cpu>(); // empty faust // no factors to the right of *(S.begin()+S.size()]-1)
			for(int i=S.size()-2;i >= 0; i--)
			{
				auto vec_Si_plus_1 = { *(S.begin()+i+1) };
				if(pR[i] != nullptr) delete pR[i];
				pR[i] = new TransformHelper<FPP,Cpu>(vec_Si_plus_1, *pR[i+1]);
				if(packing_RL) pR[i]->pack_factors();
			}
			f_id = 0;
		};
		next_fid = [&f_id, &S, &pL, &nfacts, packing_RL]()
		{
			if(f_id < nfacts-1)
			{
				if(pL[f_id+1] != nullptr)
					delete pL[f_id+1];
				auto vec_Sj = { *(S.begin()+f_id) };
				pL[f_id+1] = new Faust::TransformHelper<FPP,Cpu>(*pL[f_id], vec_Sj);
				if(packing_RL) pL[f_id+1]->pack_factors();
			}
			f_id++;
		};
		updating_facs = [&f_id, &nfacts]() {return f_id < nfacts;};
		is_last_fac_updated = [&f_id, &nfacts]() {return f_id == nfacts-1;};
	}
	Faust::MatDense<FPP,Cpu> D, tmp;
	Faust::MatDense<FPP,Cpu> * LorR;
	Faust::MatDense<FPP,Cpu> _LorR;
	while(i < nites)
	{
		//		std::cout << "nfacts:" << nfacts << std::endl;

		init_ite();
		while(updating_facs())
		{
			//			std::cout << "f_id: " << f_id << std::endl;
			cur_fac = S.get_gen_fact_nonconst(f_id);
			Real<FPP> nR=1,nL=1;
			if(pR[f_id]->size() > 0)
				nR = pR[f_id]->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
			if(pL[f_id]->size() > 0)
				nL = pL[f_id]->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
			Real<FPP> c = lipschitz_multiplicator*lambda*lambda*nR*nR*nL*nL;
			auto S_j_vec = {*(S.begin()+f_id)};
			Faust::TransformHelper<FPP, Cpu> _LSR(*pL[f_id], S_j_vec, *pR[f_id]);
//			tmp = _LSR.get_product();
			_LSR.get_product(tmp);
			tmp *= FPP(lambda);
			tmp -= A;
			FPP alpha_R = 1, alpha_L = 1, beta_R = 0, beta_L = 0; //decl in parent scope
			if(S.is_fact_sparse(f_id))
			{
				scur_fac = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(cur_fac);
				D = *scur_fac;
				dcur_fac = nullptr;
			}
			else
			{
				dcur_fac = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(cur_fac); // TOFIX: possible it is not Dense... but MatDiag (sanity check on function start)
				D = *dcur_fac;
				scur_fac = nullptr;
			}
			if(pR[f_id]->size() > 0)
			{
				if(packing_RL)
					LorR = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(pR[f_id]->get_gen_fact_nonconst(0)); //normally pR[f_id] is packed (hence reduced to a single MatDense)
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
					LorR = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(pL[f_id]->get_gen_fact_nonconst(0));
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
//			cout << "scur_fac=" << scur_fac << " dcur_fac=" << dcur_fac << endl;
			if(! use_csr && scur_fac == nullptr)
			{
				// not CSR and cur_fac DENSE
				*dcur_fac = D;
			}
			else if(use_csr && dcur_fac == nullptr)
			{
				// CSR and cur_fac SPARSE
				*scur_fac = D;
			}
			else
				throw std::runtime_error("Current factor is inconsistent with use_csr.");
			cur_fac->set_id(false);
			next_fid(); //f_id updated to iteration factor index (pL or pR too)
		}
		//update lambda
		//TODO: variable decl in parent scope
		Faust::MatDense<FPP,DEVICE> A_H_S = S.multiply(A_H);
		//		auto last_Sfac_vec = { *(S.begin()+nfacts-1), dynamic_cast<Faust::MatGeneric<FPP,Cpu>*>(&A_H)};
		//		Faust::TransformHelper<FPP,Cpu> A_H_S_(*pL[nfacts-1], last_Sfac_vec);
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
