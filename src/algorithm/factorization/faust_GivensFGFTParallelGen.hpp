
#include <cmath>
#ifdef OPT_UPDATE_L_OMP
#include <omp.h>
#endif

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
GivensFGFTParallelGen<FPP,DEVICE,FPP2,FPP4>::GivensFGFTParallelGen(int t, Faust::GivensFGFTGen<FPP, DEVICE, FPP2, FPP4> & alg) : alg(alg), t(t), fact_nrots(0)
{

	if(alg.verbosity > 1)
	{
		std::cout << "GivensFGFTGenParallelGen ctor:" << std::endl;
		std::cout << "J: " << alg.J << std::endl;
		std::cout << "tol: " << alg.stoppingError << std::endl;
		std::cout << "stopcrit is error: " << alg.stoppingCritIsError << std::endl;
		std::cout << "relErr: " << alg.errIsRel << std::endl;
		std::cout << "order: " << alg.D_order_dir << std::endl;
		std::cout << "enable_large_Faust: " << alg.enable_large_Faust << std::endl;
		auto dLap = dynamic_cast<MatDense<FPP, DEVICE>*>(&alg.Lap);
		if(dLap)
			std::cout << "matrix norm: " << dLap->norm() << std::endl;
		else
		{
			auto sLap = dynamic_cast<MatSparse<FPP, DEVICE>*>(&alg.Lap);
			if(sLap)
				std::cout << "matrix norm: " << sLap->norm() << std::endl;
		}
	}
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void GivensFGFTParallelGen<FPP,DEVICE,FPP2,FPP4>::max_L()
{
	// 	  matlab ref
	// 	  L_low = abs(tril(L,-1));
	//    ind_nnz = find(L_low);
	//    [~,ind_sorted] = sort(L_low(ind_nnz),'descend');
	//    [Rvec, Svec] = ind2sub([n,n],ind_nnz(ind_sorted));
	//    RSmat = [Rvec, Svec];
	bool is_L_sparse;
	MatSparse<FPP4,DEVICE> * sL;
	MatDense<FPP4,DEVICE> * dL = dynamic_cast<MatDense<FPP4,DEVICE>*>(alg.L);
	MatDense<FPP4,DEVICE> L_low;
	if(is_L_sparse = (dL == nullptr))
	{
		sL = dynamic_cast<MatSparse<FPP4,DEVICE>*>(alg.L);
		L_low = *sL;
	}
	else
		L_low = *dL;
	L_low.abs();
	L_low =	L_low.lower_tri(false);
	fact_nz_inds = L_low.nonzeros_indices();


#ifdef DEBUG_GIVENS
	for(auto &p : fact_nz_inds)
		cout << "GivensFGFTParallel::max_L() before sort (" << p.first+1 <<  "," << p.second+1 << ") :" << L_low(p.first, p.second) << endl;
#endif

	//ENOTE: can't use sort(it,it, lambda) as vector because list doesn't provide random access it.
	fact_nz_inds.sort([&](const pair<int,int> &a, const pair<int,int> &b)
			{
			return this->fact_nz_inds_sort_func(a,b, L_low);
			});
#ifdef DEBUG_GIVENS
	for(auto &p : fact_nz_inds)
		cout << "GivensFGFTParallel::max_L() after sort (" << p.first+1 <<  "," << p.second+1 << ") :" << L_low(p.first, p.second) << endl;
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void GivensFGFTParallelGen<FPP,DEVICE,FPP2, FPP4>::update_fact_nz_inds(int p, int q)
{
	// matlab ref
	//	other = RSmat==p | RSmat==q;
	//	tokill = any(other,2);
	//	%RSmat( tokill, : ) = [];
	//	RSmatnew = RSmat(~tokill,:);
	//	RSmat = RSmatnew;

	//remove all pairs containing p or q
	for(auto i = fact_nz_inds.begin(); i != fact_nz_inds.end();)
	{
		if((*i).first == p || (*i).second == q || (*i).first == q || (*i).second == p)
			i = fact_nz_inds.erase(i);
		else
			i++;
	}
#ifdef DEBUG_GIVENS
	std::cout << "GivensFGFTParallel::update_fact_nz_inds() after purge: ";
	for(auto &p_ : fact_nz_inds)
		std::cout << "(" << p_.first+1 <<  "," << p_.second+1 << ") :";
	std::cout <<  std::endl;
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void GivensFGFTParallelGen<FPP,DEVICE,FPP2, FPP4>::loop_update_fact()
{
	fact_nrots = 0;
	while(fact_nrots < t && 0 < fact_nz_inds.size())
	{
		choose_pivot();
		update_fact_nz_inds(alg.p, alg.q);
		alg.calc_theta();
		alg.update_fact();
		fact_nrots++;
	}
	finish_fact();
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void GivensFGFTParallelGen<FPP,DEVICE,FPP2, FPP4>::choose_pivot()
{
	// fact_nz_inds is assumed to be sorted in max_L()
	// the first one in the list maps the biggest abs val in L
	pair<int,int> max_elt = *(fact_nz_inds.begin());
	alg.p = max_elt.first;
	alg.q = max_elt.second;
	alg.coord_choices.push_back(std::pair<int,int>(alg.p, alg.q));
#ifdef DEBUG_GIVENS
	cout << "choose_pivot() alg.p: " << alg.p+1 << " alg.q:" << alg.q+1 << " " << "L(alg.p,alg.q): " << alg.L(alg.p,alg.q) << " nrots: " << fact_nrots << endl;
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void GivensFGFTParallelGen<FPP,DEVICE,FPP2,FPP4>::update_L(Faust::MatSparse<FPP,Cpu> & L)
{
	// L = S'*L*S
#ifdef OPT_UPDATE__SPARSE_L
	int choice_id;
	FPP2 c,s;
	int p, q, i;
	// L = S'*L
// HINT to enable OpenMP: add flag -fopenmp (to compilation and linker flags in CMakeCache)
// likewise in setup.py for python wrapper (into extra_compile_args and extra_link_args list -- setuptools)
// NOTE: OpenMP was an attempt not really useful because no speedup was noticed in two or four threads (running on 4-core CPU)
// The code is nevertheless kept here, just in case, to use it you need to set the compilation/preprocessor constant OPT_UPDATE_L_OMP
#ifdef OPT_UPDATE_L_OMP
#pragma omp parallel private(c, s, choice_id, p, q, i)
#endif
	{
		Eigen::SparseMatrix<FPP, Eigen::RowMajor> L_vec_p, L_vec_q;
#ifdef OPT_UPDATE_L_OMP
		int num_per_th = fact_nrots/omp_get_num_threads();
		int th_id = omp_get_thread_num();
		int th_offset = fact_nrots/omp_get_num_threads()*th_id;
		if(th_id == omp_get_num_threads() - 1) num_per_th += fact_nrots % omp_get_num_threads();
#else
		int num_per_th = fact_nrots; // a single thread is used (no openmp)
		int th_id = 0;
		int th_offset = 0;
#endif
//first attempt to openmp-ize the loops (abandoned for the single section started above)
//#ifdef OPT_UPDATE_L_OMP
//#pragma omp parallel for schedule(static) private(c, s, choice_id, p, q, i) default(shared) // num_threads(4)
//#endif OPT_UPDATE_L_OMP
		for(i=th_offset; i < th_offset+num_per_th; i++)
		{
			// applying first the last rotation submatrix
			choice_id = this->coord_choices.size()-1-i;
			c = *(this->fact_mod_values.end()-1-4*i); // cos(theta)
			s = *(this->fact_mod_values.end()-2*(2*i+1)); // sin(theta)
			p = this->coord_choices[choice_id].first;
			q = this->coord_choices[choice_id].second;
			//rely on parent for doing the job
			this->update_L_first(L_vec_p, L_vec_q, c, s, p, q, L);
		}
#ifdef OPT_UPDATE_L_OMP
#pragma omp barrier //S'*L must be computed before starting L*S computation
#endif
		// L = L*S
//#ifdef OPT_UPDATE_L_OMP
//#pragma omp parallel for schedule(static) private(c, s, choice_id, p, q, i) default(shared) // num_threads(4)
//#endif
//		for(i=th_offset; i < th_offset+num_per_th; i++)
//		{
//			// applying first the last rotation submatrix
//			choice_id = this->coord_choices.size()-1-i;
//			c = *(this->fact_mod_values.end()-1-4*i); // cos(theta)
//			s = *(this->fact_mod_values.end()-2*(2*i+1)); // sin(theta)
//			p = this->coord_choices[choice_id].first;
//			q = this->coord_choices[choice_id].second;
//			//rely on parent for doing the job
//			this->update_L_second(L_vec_p, L_vec_q, c, s, p, q, L);
//			L.update_dim();
//		}
		L.multiplyRight(this->facts[this->ite]);
//		L.update_dim();
	}
#else
	// L = S'*L*S
	this->facts[this->ite].multiply(L, 'T');
	L.multiplyRight(this->facts[this->ite]);
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void GivensFGFTParallelGen<FPP,DEVICE,FPP2,FPP4>::finish_fact()
{
	int n = alg.Lap.getNbRow();
	alg.facts[alg.ite] = Faust::MatSparse<FPP4,DEVICE>(alg.fact_mod_row_ids, alg.fact_mod_col_ids, alg.fact_mod_values, n, n);
	alg.facts[alg.ite].set_orthogonal(true);
#ifdef DEBUG_GIVENS
	cout << "GivensFGFTParallel::finish_fact() ite: " << alg.ite << " fact norm: " << alg.facts[alg.ite].norm() << endl;
	alg.facts[alg.ite].Display();
#endif
}
