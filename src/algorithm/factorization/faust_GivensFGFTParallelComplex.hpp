using namespace Faust;

#include <cmath>
#ifdef OPT_UPDATE_L_OMP
#include <omp.h>
#endif

template<typename FPP, Device DEVICE, typename FPP2>
GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::GivensFGFTParallelComplex(Faust::MatDense<FPP,DEVICE>& Lap, int J, int t, unsigned int verbosity, const double stoppingError, const bool errIsRel) : GivensFGFTComplex<FPP,DEVICE,FPP2>(Lap, J, verbosity, stoppingError, errIsRel), t(t), fact_nrots(0)
{
	this->facts.resize(round(J/(float)t));
	this->coord_choices.resize(0);
}

template<typename FPP, Device DEVICE, typename FPP2>
GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::GivensFGFTParallelComplex(Faust::MatSparse<FPP,DEVICE>& Lap, int J, int t, unsigned int verbosity, const double stoppingError, const bool errIsRel) : GivensFGFTComplex<FPP,DEVICE,FPP2>(Lap, J, verbosity, stoppingError, errIsRel), t(t), fact_nrots(0)
{
	this->facts.resize(round(J/(float)t));
	this->coord_choices.resize(0);
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::max_L()
{
	// 	  matlab ref
	// 	  L_low = abs(tril(L,-1));
	//    ind_nnz = find(L_low);
	//    [~,ind_sorted] = sort(L_low(ind_nnz),'descend');
	//    [Rvec, Svec] = ind2sub([n,n],ind_nnz(ind_sorted));
	//    RSmat = [Rvec, Svec];
	bool is_L_sparse;
	MatSparse<FPP,DEVICE> * sL;
	MatDense<FPP,DEVICE> * dL = dynamic_cast<MatDense<FPP,DEVICE>*>(this->L);
	MatDense<FPP,DEVICE> L_low;
	if(is_L_sparse = (dL == nullptr))
	{
		sL = dynamic_cast<MatSparse<FPP,DEVICE>*>(this->L);
		L_low = *sL;
	}
	else
		L_low = *dL;
	L_low.abs();
	L_low =	L_low.lower_tri(false);
	fact_nz_inds = L_low.nonzeros_indices();


#ifdef DEBUG_GIVENS
	for(auto &p : fact_nz_inds)
		cout << "GivensFGFTParallelComplex::max_L() before sort (" << p.first+1 <<  "," << p.second+1 << ") :" << L_low(p.first, p.second) << endl;
#endif

	//ENOTE: can't use sort(it,it, lambda) as vector because list doesn't provide random access it.
	fact_nz_inds.sort([&L_low](const pair<int,int> &a, const pair<int,int> &b)
			{
			FPP2 a_ = Faust::fabs(L_low(a.first, a.second));
			FPP2 b_ = Faust::fabs(L_low(b.first, b.second));
			if(isnan(b_)) return true; //select a
			else if(isnan(a_)) return false; //select b
			else return Faust::fabs(L_low(a.first, a.second)) > Faust::fabs(L_low(b.first, b.second));
			});
#ifdef DEBUG_GIVENS
	for(auto &p : fact_nz_inds)
		cout << "GivensFGFTParallelComplex::max_L() after sort (" << p.first+1 <<  "," << p.second+1 << ") :" << L_low(p.first, p.second) << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::update_fact_nz_inds()
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
		if((*i).first == this->p || (*i).second == this->q || (*i).first == this->q || (*i).second == this->p)
			i = fact_nz_inds.erase(i);
		else
			i++;
	}
#ifdef DEBUG_GIVENS
	for(auto &p : fact_nz_inds)
		cout << "GivensFGFTParallelComplex::update_fact_nz_inds() after purge (" << p.first+1 <<  "," << p.second+1 << ") :" << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::next_step()
{

	substep_fun substep[] = {
		&GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::max_L,
		&GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::loop_update_fact, //responsible to call choose_pivot(), calc_theta() and update_fact()
		&GivensFGFTComplex<FPP,DEVICE,FPP2>::update_L,
		&GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::update_D,
		&GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::update_err};

	for(int i=0;i<sizeof(substep)/sizeof(substep_fun);i++)
	{
#ifdef DEBUG_GIVENS
		cout << "GivensFGFTParallelComplex ite=" << this->ite << " substep i=" << i << endl;
#endif
		(this->*substep[i])();
	}
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::loop_update_fact()
{
	fact_nrots = 0;
	while(fact_nrots < t && 0 < fact_nz_inds.size())
	{
		choose_pivot();
		update_fact_nz_inds();
		this->calc_theta();
		this->update_fact();
		fact_nrots++;
	}
	finish_fact();
#ifdef DEBUG_GIVENS
	int n = this->L->getNbRow();
//	cout << "n=" << n << endl;
	MatSparse<FPP,DEVICE> test1(this->facts[this->ite]);
	MatSparse<FPP,DEVICE> test2(this->facts[this->ite]);
	test2.conjugate();
	test2.transpose();
	test1.multiply(test2, 'N');
	for(int j = 0; j < n; j++) {
		//		for(int k = 0; k < n; k++)
//		cout << "ite=" << this->ite << "S*S'(" << j << "," << j << ")=" << test2(j,j) << endl;
		assert(Faust::abs(test2(j,j)-FPP(1,0)) < 1e-3);
	}
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::choose_pivot()
{
	// fact_nz_inds is assumed to be sorted in max_L()
	// the first one in the list maps the biggest abs val in L
	pair<int,int> max_elt = *(fact_nz_inds.begin());
	this->p = max_elt.first;
	this->q = max_elt.second;
//	cout << "max_elt=" << (*(this->L))(this->p,this->q) << endl;
	this->coord_choices.push_back(pair<int,int>(this->p,this->q));
#ifdef DEBUG_GIVENS
	cout << "choose_pivot() p: " << this->p+1 << " q:" << this->q+1 << " " << "L(p,q): " << (*(this->L))(this->p,this->q) << " nrots: " << fact_nrots << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::update_fact()
{
	FPP c_pp, c_pq, c_qp, c_qq;
	FPP i = complex<typename FPP::value_type>(0,1);
	FPP sin_theta2, cos_theta2;
	sin_theta2 = sin(this->theta2);
	cos_theta2 = cos(this->theta2);
	c_pp = - i * exp(-i*this->theta1) * sin_theta2;
	c_pq = - exp(i*this->theta1) * cos_theta2;
	c_qp = exp(-i*this->theta1) * cos_theta2;
	c_qq = i*exp(i*this->theta1) * sin_theta2;

	FPP tmp = c_pq;
	c_pp = conj(c_pp);
	c_qq = conj(c_qq);
	c_pq = conj(c_qp);
	c_qp = conj(tmp);

	this->check_pivot_image(c_pp, c_pq, c_qp, c_qq);

	assert(!std::isnan(Faust::fabs(c_pp)));
	if(fact_nrots == 0)
	{
		int n = this->Lap.getNbRow();
		// forget previous rotation coeffs
		// and keep identity part (n first coeffs)
		this->fact_mod_row_ids.resize(n);
		this->fact_mod_col_ids.resize(n);
		this->fact_mod_values.resize(n);
	}
	// write new coeffs
	// 1st one
	this->fact_mod_row_ids.push_back(this->p);
	this->fact_mod_col_ids.push_back(this->p);
	this->fact_mod_values.push_back(c_pp);
	// 2nd
	this->fact_mod_row_ids.push_back(this->p);
	this->fact_mod_col_ids.push_back(this->q);
	this->fact_mod_values.push_back(c_pq);
	// 3rd
	this->fact_mod_row_ids.push_back(this->q);
	this->fact_mod_col_ids.push_back(this->p);
	this->fact_mod_values.push_back(c_qp);
	// 4th
	this->fact_mod_row_ids.push_back(this->q);
	this->fact_mod_col_ids.push_back(this->q);
	this->fact_mod_values.push_back(c_qq);
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::update_L(Faust::MatDense<FPP,Cpu> & L)
{
	// L = S'*L*S
//#undef OPT_UPDATE_L
//#ifdef OPT_UPDATE_L
//	int choice_id;
//	FPP  c_pp, c_pq, c_qp, c_qq;
//	int p, q, i;
//	// L = S'*L
//// HINT to enable OpenMP: add flag -fopenmp (to compilation and linker flags in CMakeCache)
//// likewise in setup.py for python wrapper (into extra_compile_args and extra_link_args list -- setuptools)
//// NOTE: OpenMP was an attempt not really useful because no speedup was noticed in two or four threads (running on 4-core CPU)
//// The code is nevertheless kept here, just in case, to use it you need to set the compilation/preprocessor constant OPT_UPDATE_L_OMP
//#ifdef OPT_UPDATE_L_OMP
//#pragma omp parallel private(c_pp, c_pq, c_qp, c_qq, choice_id, p, q, i)
//#endif
//	{
//		Vect<FPP,DEVICE> L_vec_p, L_vec_q;
//#ifdef OPT_UPDATE_L_OMP
//		int num_per_th = fact_nrots/omp_get_num_threads();
//		int th_id = omp_get_thread_num();
//		int th_offset = fact_nrots/omp_get_num_threads()*th_id;
//		if(th_id == omp_get_num_threads() - 1) num_per_th += fact_nrots % omp_get_num_threads();
//#else
//		int num_per_th = fact_nrots; // a single thread is used (no openmp)
//		int th_id = 0;
//		int th_offset = 0;
//#endif
////first attempt to openmp-ize the loops (abandoned for the single section started above)
////#ifdef OPT_UPDATE_L_OMP
////#pragma omp parallel for schedule(static) private(c, s, choice_id, p, q, i) default(shared) // num_threads(4)
////#endif OPT_UPDATE_L_OMP
//		for(i=th_offset; i < th_offset+num_per_th; i++)
//		{
//			// applying first the last rotation submatrix
//			choice_id = this->coord_choices.size()-1-i;
//			c_qq = *(this->fact_mod_values.end()-1-4*i);
//			c_qp = *(this->fact_mod_values.end()-2*(2*i+1));
//			c_pq = *(this->fact_mod_values.end()-4*i-3);
//			c_pp = *(this->fact_mod_values.end()-4*i-4);
//			p = this->coord_choices[choice_id].first;
//			q = this->coord_choices[choice_id].second;
//			//rely on parent for doing the job
//			this->update_L_first(L_vec_p, L_vec_q, c_pp, c_pq, c_qp, c_qq, p, q, L);
//		}
//#ifdef OPT_UPDATE_L_OMP
//#pragma omp barrier //S'*L must be computed before starting L*S computation
//#endif
//		// L = L*S
////#ifdef OPT_UPDATE_L_OMP
////#pragma omp parallel for schedule(static) private(c, s, choice_id, p, q, i) default(shared) // num_threads(4)
////#endif
//		for(i=th_offset; i < th_offset+num_per_th; i++)
//		{
//			// applying first the last rotation submatrix
//			choice_id = this->coord_choices.size()-1-i;
//			c_qq = *(this->fact_mod_values.end()-1-4*i);
//			c_qp = *(this->fact_mod_values.end()-2*(2*i+1));
//			c_pq = *(this->fact_mod_values.end()-4*i-3);
//			c_pp = *(this->fact_mod_values.end()-4*i-4);
//			//rely on parent for doing the job
//			this->update_L_second(L_vec_p, L_vec_q, c_pp, c_pq, c_qp, c_qq, p, q, L);
//		}
//	}
//#else
	this->facts[this->ite].conjugate();
	this->facts[this->ite].multiply(L, 'T');
	this->facts[this->ite].conjugate();
	L.multiplyRight(this->facts[this->ite]);
#ifdef DEBUG_GIVENS
	cout << "L(p,q) after update_L():" << L(this->p,this->q) << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::update_L(Faust::MatSparse<FPP,Cpu> & L)
{
	// L = S'*L*S
#ifdef OPT_UPDATE_SPARSE_L
	int choice_id;
	int p, q, i;
	FPP c_qq, c_qp, c_pq, c_pp;
	// L = S'*L
// HINT to enable OpenMP: add flag -fopenmp (to compilation and linker flags in CMakeCache)
// likewise in setup.py for python wrapper (into extra_compile_args and extra_link_args list -- setuptools)
// NOTE: OpenMP was an attempt not really useful because no speedup was noticed in two or four threads (running on 4-core CPU)
// The code is nevertheless kept here, just in case, to use it you need to set the compilation/preprocessor constant OPT_UPDATE_L_OMP
#ifdef OPT_UPDATE_L_OMP
#pragma omp parallel private(c_pp, c_pq, c_qp, c_qq, choice_id, p, q, i)
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
			c_qq = *(this->fact_mod_values.end()-1-4*i);
			c_qp = *(this->fact_mod_values.end()-2*(2*i+1));
			c_pq = *(this->fact_mod_values.end()-4*i-3);
			c_pp = *(this->fact_mod_values.end()-4*i-4);
			p = this->coord_choices[choice_id].first;
			q = this->coord_choices[choice_id].second;
			//rely on parent for doing the job
			this->update_L_first(L_vec_p, L_vec_q, c_pp, c_pq, c_qp, c_qq, p, q, L);
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
	this->facts[this->ite].multiply(L, 'H');
	L.multiplyRight(this->facts[this->ite]);
#endif
}


template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::finish_fact()
{
	int n = this->Lap.getNbRow();
	this->facts[this->ite] = MatSparse<FPP,DEVICE>(this->fact_mod_row_ids, this->fact_mod_col_ids, this->fact_mod_values, n, n);
	this->facts[this->ite].set_orthogonal(true);
#ifdef DEBUG_GIVENS
	cout << "GivensFGFTParallelComplex::finish_fact() ite: " << this->ite << " fact norm: " << this->facts[this->ite].norm() << endl;
	this->facts[this->ite].Display();
#endif
}
