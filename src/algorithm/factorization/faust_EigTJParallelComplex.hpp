namespace Faust
{
	template<typename FPP, FDevice DEVICE, typename FPP2>
		EigTJParallelComplex<FPP,DEVICE,FPP2>::EigTJParallelComplex(MatDense<FPP,DEVICE>& Lap, int J, int t, unsigned int verbosity, const double stoppingError, const bool errIsRel,  const bool enable_large_Faust, const int err_period/*=100*/) : EigTJComplex<FPP,DEVICE,FPP2>(Lap, J, verbosity, stoppingError, errIsRel, enable_large_Faust, err_period), /*t(t), fact_nrots(0)*/ EigTJParallelGen<typename FPP::value_type, DEVICE, FPP2, FPP>(t, *this)
	{
		if(J > 0) this->facts.resize(round(J/(float)t));
		this->coord_choices.resize(0);
		init_fact_nz_inds_sort_func();
	}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		EigTJParallelComplex<FPP,DEVICE,FPP2>::EigTJParallelComplex(MatSparse<FPP,DEVICE>& Lap, int J, int t, unsigned int verbosity, const double stoppingError, const bool errIsRel, const bool enable_large_Faust, const int err_period/*=100*/) : EigTJComplex<FPP,DEVICE,FPP2>(Lap, J, verbosity, stoppingError, errIsRel, enable_large_Faust, err_period), /* t(t), fact_nrots(0) */ EigTJParallelGen<typename FPP::value_type, DEVICE, FPP2, FPP>(t, *this)
	{
		if(J > 0) this->facts.resize(round(J/(float)t));
		this->coord_choices.resize(0);
		init_fact_nz_inds_sort_func();
	}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void EigTJParallelComplex<FPP,DEVICE,FPP2>::init_fact_nz_inds_sort_func()
		{
			this->fact_nz_inds_sort_func= [](const pair<int,int> &a, const pair<int,int> &b, MatDense<FPP,DEVICE> & L_low)
			{
				FPP2 a_ = Faust::fabs(L_low(a.first, a.second));
				FPP2 b_ = Faust::fabs(L_low(b.first, b.second));
				if(isnan(b_)) return true; //select a
				else if(isnan(a_)) return false; //select b
				else return Faust::fabs(L_low(a.first, a.second)) > Faust::fabs(L_low(b.first, b.second));
			};
		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void EigTJParallelComplex<FPP,DEVICE,FPP2>::next_step()
		{

			substep_fun substep[] = {
				&EigTJParallelGen<typename FPP::value_type,DEVICE,FPP2,FPP>::max_L,
				&EigTJParallelComplex<FPP,DEVICE,FPP2>::loop_update_fact, //responsible to call choose_pivot(), calc_theta() and update_fact()
				&EigTJComplex<FPP,DEVICE,FPP2>::update_L,
				&EigTJParallelComplex<FPP,DEVICE,FPP2>::update_D,
				&EigTJParallelComplex<FPP,DEVICE,FPP2>::update_err};

			for(int i=0;i<sizeof(substep)/sizeof(substep_fun);i++)
			{
#ifdef DEBUG_GIVENS
				cout << "EigTJParallelComplex ite=" << this->ite << " substep i=" << i << endl;
#endif
				(this->*substep[i])();
			}
			this->ite++;
		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void EigTJParallelComplex<FPP,DEVICE,FPP2>::update_fact()
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
			if(this->fact_nrots == 0)
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
			if(this->J <= this->ite+1) this->facts.resize(this->ite+1);
		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void EigTJParallelComplex<FPP,DEVICE,FPP2>::update_L(MatDense<FPP,Cpu> & L)
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
}
