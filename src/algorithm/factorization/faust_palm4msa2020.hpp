template <typename FPP, FDevice DEVICE>
void Faust::palm4msa(const Faust::MatDense<FPP,DEVICE>& A,
		std::vector<Faust::ConstraintGeneric*> & constraints,
		Faust::TransformHelper<FPP,DEVICE>& S,
		FPP& lambda,
		//const unsigned int nites,
		const StoppingCriterion<Real<FPP>>& sc,
		const bool is_update_way_R2L,
		const bool use_csr,
		const bool compute_2norm_on_array,
		const Real<FPP> norm2_threshold,
		const unsigned int norm2_max_iter,
		const bool constant_step_size, const Real<FPP> step_size,
		const bool on_gpu /*=false*/,
		const bool is_verbose/*=false*/)
{
	if(constraints.size() == 0)
		throw out_of_range("No constraint passed to palm4msa.");
	const Real<FPP> lipschitz_multiplicator = 1.001;
	Faust::MatGeneric<FPP,DEVICE>* cur_fac;
	Faust::MatSparse<FPP,DEVICE>* scur_fac;
	Faust::MatDense<FPP,DEVICE>* dcur_fac;
	Faust::MatSparse<FPP,DEVICE> spD;
	Real<FPP> error = -1; // negative error is ignored
	const unsigned int nfacts = constraints.size();
	std::vector<std::pair<faust_unsigned_int,faust_unsigned_int>> dims;
	int norm2_flag; // return val
	for(auto c: constraints)
		dims.push_back(make_pair(c->get_rows(), c->get_cols()));
	Faust::MatDense<FPP,DEVICE> A_H = A;
	A_H.adjoint();
	if(S.size() != nfacts)
		fill_of_eyes(S, nfacts, use_csr, dims, on_gpu);
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
	Real<FPP> c = 1/step_size;
	// lambda exp to update fact when its id is 0
	auto update_1stfac = [&sc, &error, &constant_step_size, &c, &A, &D, &tmp, &LSR, &scur_fac, &dcur_fac, &f_id, &S, &lipschitz_multiplicator, &lambda, &norm2_threshold, &norm2_flag, &norm2_max_iter, &on_gpu](Faust::MatGeneric<FPP, DEVICE> *cur_fac)
	{
		auto R = S.right(f_id+1);
		if(! constant_step_size)
		{
			Real<FPP> nR = R->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
			c = lipschitz_multiplicator*lambda*lambda*nR*nR;
		}
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
		Faust::TransformHelper<FPP, Cpu> _LSR(L, *R); // L contains cur_fac
//		LSR = _LSR;
//		_LSR.multiply(lambda);
//		tmp = _LSR.get_product();
		_LSR.get_product(tmp);
		tmp *= lambda;
		tmp -= A;
		if(sc.isCriterionErr())
			error = tmp.norm();
		//TODO: do something to lighten the double transpose conjugate
		tmp.adjoint();
		tmp = R->multiply(tmp, /* H */ false, false);
		tmp.adjoint();
		tmp *= lambda/c;
		D -= tmp;
	};
	auto update_lastfac = [&sc, &error, &constant_step_size, &c, &A, &D, &tmp, &LSR, &scur_fac, &dcur_fac, &f_id, &S, &lipschitz_multiplicator, &lambda, &norm2_threshold, &norm2_flag, &norm2_max_iter, &on_gpu](Faust::MatGeneric<FPP, DEVICE> *cur_fac)
	{
		auto L = S.left(f_id-1);
		// TODO: factorize with other lambda exp
		if(! constant_step_size)
		{
			Real<FPP> nL = L->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
			c = lipschitz_multiplicator*lambda*lambda*nL*nL;
		}
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
		Faust::TransformHelper<FPP, Cpu> _LSR(*L, R); // L contains cur_fac
//		LSR = _LSR;
//		_LSR.multiply(lambda);
		tmp = _LSR.get_product();
		tmp *= lambda;
		tmp -= A;
		if(sc.isCriterionErr())
			error = tmp.norm();
		tmp = L->multiply(tmp, /* NO H */ true, true);
		tmp *= lambda/c;
		D -= tmp;
	};
	auto update_interfac = [&sc, &error, &constant_step_size, &c, &A, &D, &tmp, &LSR, &scur_fac, &dcur_fac, &f_id, &S, &lipschitz_multiplicator, &lambda, &norm2_threshold, &norm2_flag, &norm2_max_iter, &on_gpu](Faust::MatGeneric<FPP, DEVICE> *cur_fac)
	{
		auto R = S.right(f_id+1);
		auto L = S.left(f_id-1);
		if(! constant_step_size)
		{
			Real<FPP> nR = R->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
			Real<FPP> nL = L->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
			c = lipschitz_multiplicator*lambda*lambda*nR*nR*nL*nL;
		}
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
		Faust::TransformHelper<FPP, Cpu> _LSR(*L, *R); // L contains cur_fac
//		LSR = _LSR;
//		_LSR.multiply(lambda);
		tmp = _LSR.get_product();
		tmp *= lambda;
		tmp -= A;
		if(sc.isCriterionErr())
			error = tmp.norm();
		//TODO: do something to lighten the double transpose conjugate
		tmp.adjoint();
		tmp = R->multiply(tmp, /* NO H */ false, false);
		tmp.adjoint();
		tmp = L->multiply(tmp, true, true);
		tmp *= lambda/c;
		D -= tmp;
	};
	while(sc.do_continue(i, error))
	{
		//		std::cout << "nfacts:" << nfacts << std::endl;
		init_fid();
		while(updating_facs())
		{
			//			std::cout << "f_id: " << f_id << std::endl;
			cur_fac = S.get_gen_fact_nonconst(f_id);
			if(f_id == 0)
			{
				update_1stfac(cur_fac);
			}
			else if(f_id == nfacts-1)
				update_lastfac(cur_fac);
			else
				update_interfac(cur_fac);
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
				S.update(D, f_id);
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

	template <typename FPP, FDevice DEVICE>
void Faust::fill_of_eyes(
		Faust::TransformHelper<FPP,DEVICE>& S,
		const unsigned int nfacts,
		const bool sparse,
		const std::vector<std::pair<faust_unsigned_int,faust_unsigned_int>> dims,
		const bool on_gpu)
{
	if(S.size() > 0)
		throw std::runtime_error("The Faust must be empty for intializing it to eyes factors.");
	for(auto fdims : dims)
	{
//		std::cout << "fill_of_eyes() fdims: " << fdims.first << " " << fdims.second << std::endl;
		// init all facts as identity matrices
		// with proper dimensions
		Faust::MatGeneric<FPP,DEVICE>* fact;
		if(sparse)
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
