using namespace Faust;

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::next_step()
{

	substep_fun substep[] = {
		&EigTJComplex<FPP,DEVICE,FPP2>::max_L,
		&EigTJComplex<FPP,DEVICE,FPP2>::choose_pivot,
		&EigTJComplex<FPP,DEVICE,FPP2>::calc_theta,
		&EigTJComplex<FPP,DEVICE,FPP2>::update_fact,
		&EigTJComplex<FPP,DEVICE,FPP2>::update_L,
		&EigTJComplex<FPP,DEVICE,FPP2>::update_D,
		&EigTJComplex<FPP,DEVICE,FPP2>::update_err};

	for(int i=0;i<sizeof(substep)/sizeof(substep_fun);i++)
	{
#ifdef DEBUG_GIVENS
		cout << "EigTJComplex ite=" << this->ite << " substep i=" << i << endl;
#endif
		(this->*substep[i])();
	}
	this->ite++;
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::choose_pivot()
{
//Matlab ref. code:
//        [~,p] = min(C_min_row);
//        q = q_candidates(p);
//        coord_choices(1,j) = p;
//        coord_choices(2,j) = q;
	C_min_row.max_coeff(&(this->p));
	this->q = this->q_candidates[this->p];
	this->coord_choices.push_back(pair<int,int>(this->p,this->q));
#ifdef DEBUG_GIVENS
	cout << "EigTJComplex::choose_pivot() ite: " << this->ite << " p: " << this->p << " q: " << this->q << endl;
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::max_L()
{
	// Matlab ref. code:
	//**************  at initialization
	//    for r=1:n
	//        for s=r+1:n
	//            C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
	//        end
	//    end
	//    [C_min_row, q_candidates] = min(C,[],2);
	//
	int n = this->dim_size, r, s;
	if(!this->ite)
	{
		//first iteration
		for(r=0;r<n;r++)
			for(s=r+1;s<n;s++)
			{
				//TODO: really slow (element-by-element copy) see if we can use a max and directly copy L per block
				// MatDense's data is in column-major order
				C.getData()[s*n+r] = Faust::fabs((*this->L)(r,s));
			}
		C_min_row = C.rowwise_max(this->q_candidates);
	}

	/************ at end of ite.
	 *       for r=[p,q]
	 *           for s=r+1:n
	 *               C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
	 *           end
	 *           [C_min_row(r), q_candidates(r)] = min(C(r,:));
	 *       end
	 *
	 * 		for s=[p,q]
	 *            for r=1:s-1
	 *                C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
	 *                if C(r,s)<C_min_row(r)
	 *                    C_min_row(r) = C(r,s);
	 *                    q_candidates(r) = s;
	 *                elseif q_candidates(r) == s
	 *                    [C_min_row(r), q_candidates(r)] = min(C(r,:));
	 *                end
	 *            end
	 *        end
	 */
	else
	{ // 2nd to last iteration
		int pq[2] = { this->p , this->q };
		int rid;
		for(int i=0;i<2;i++)
		{
			r = pq[i];
			for(int s=r+1;s<n;s++)
				C.getData()[s*n+r] = Faust::fabs((*this->L)(r,s));
			C_min_row.getData()[r] = C.get_row(r).max_coeff(&rid);
			this->q_candidates[r] = rid;
		}
		for(int i=0;i<2;i++)
		{
			s = pq[i];
			for(r=0;r<s-1;r++)
			{
				C.getData()[s*n+r] =  Faust::fabs((*this->L)(r,s));
				if(Faust::fabs(C(r,s)) > Faust::fabs(C_min_row[r]))
				{
					C_min_row[r] = C(r,s);
					this->q_candidates[r] = s;
				}
				else if(this->q_candidates[r] == s)
				{
					C_min_row.getData()[r] = C.get_row(r).max_coeff(&rid);
					this->q_candidates[r] = rid;
				}
			}
		}
	}
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::calc_theta()
{
	FPP phi1, phi2;

	phi1 = atan((*this->L)(this->p,this->q).imag() / (*this->L)(this->p,this->q).real());
	phi2 = atan(2*Faust::fabs((*this->L)(this->p,this->q)) / ((*this->L)(this->p,this->p) - (*this->L)(this->q,this->q)));


	theta1 = (FPP(2)*phi1 - FPP(M_PI))/FPP(4);
	theta2 = phi2/FPP(2);
#ifdef DEBUG_GIVENS
	cout << "theta1=" << theta1 << "theta2=" << theta2 << endl;
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::check_pivot_image(FPP& c_pp, FPP& c_pq, FPP& c_qp, FPP& c_qq)
{
	FPP im_pivot_pq = (conj(c_pp)*(*this->L)(this->p,this->p)+conj(c_qp)*(*this->L)(this->q,this->p))*c_pq +(conj(c_pp)*(*this->L)(this->p,this->q)+conj(c_qp)*(*this->L)(this->q,this->q))*c_qq;

#ifdef DEBUG_GIVENS
	cout << "First value of theta2 gives pivot image im_L(p,q)=" << im_pivot_pq << endl;
#endif
	if(Faust::fabs(im_pivot_pq) > 1e-3)
	{
		c_pp = - c_pp;
		c_qq = - c_qq;
	}
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::update_fact()
{
	// Matlab ref. code:
	//        S = eye(n);
	//        S(p,p) = cos(theta); S(p,q) = -sin(theta);
	//        S(q,p) = sin(theta); S(q,q) = cos(theta);
	//        S = sparse(S);
	//        facts{j} = S;
	//
	int n = this->dim_size;
	FPP tmp, sin_theta2, cos_theta2;
	FPP c_pp, c_pq, c_qp, c_qq;
	FPP i = complex<typename FPP::value_type>(0,1); //TODO: replace by 1i if C++14
	//ref: https://en.cppreference.com/w/cpp/numeric/complex/operator%22%22i
	sin_theta2 = sin(theta2);
	cos_theta2 = cos(theta2);
	c_pp = - i * exp(-i*theta1) * sin_theta2;
	c_pq = - exp(i*theta1) * cos_theta2;
	c_qp = exp(-i*theta1) * cos_theta2;
	c_qq = i*exp(i*theta1) * sin_theta2;

	tmp = c_pq;
//	c_pp = complex<typename FPP::value_type>(c_pp.real(), - c_pp.imag());
	c_pp = conj(c_pp);
//	c_pq = complex<typename FPP::value_type>(c_qp.real(), - c_qp.imag());
	c_pq = conj(c_qp);
//	c_qp = complex<typename FPP::value_type>(tmp.real(), - tmp.imag());
	c_qp = conj(tmp);
//	c_qq = complex<typename FPP::value_type>(c_qq.real(), - c_qq.imag());
	c_qq = conj(c_qq);
#ifdef DEBUG_GIVENS
	cout << "c_pp=" << c_pp << "c_pq=" << c_pq << endl;
	cout << "c_qp=" << c_qp << "c_qq=" << c_qq << endl;
#endif
	check_pivot_image(c_pp, c_pq, c_qp, c_qq);

	// forget previous rotation coeffs
	// and keep identity part (n first coeffs)
	this->fact_mod_row_ids.resize(n);
	this->fact_mod_col_ids.resize(n);
	this->fact_mod_values.resize(n);
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
	this->facts[this->ite] = MatSparse<FPP,DEVICE>(this->fact_mod_row_ids, this->fact_mod_col_ids, this->fact_mod_values, n, n);
	this->facts[this->ite].set_orthogonal(true);
#ifdef DEBUG_GIVENS
	MatSparse<FPP,DEVICE> test1(this->facts[this->ite]);
	MatSparse<FPP,DEVICE> test2(this->facts[this->ite]);
	test2.conjugate();
	test2.transpose();
	test1.multiply(test2, 'N');
	for(int j = 0; j < n; j++) {
//		for(int k = 0; k < n; k++)
//			cout << "ite=" << ite << "S*S'(" << j << "," << j << ")=" << test2(j,j) << endl;
			assert(Faust::abs(test2(j,j)-FPP(1,0)) < .01);
	}
	cout << "EigTJComplex::update_fact() ite: " << this->ite << " fact norm: " << this->facts[this->ite].norm() << endl;
	this->facts[this->ite].Display();
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::update_L(Faust::MatDense<FPP,Cpu> & L)
{
	// L = S'*L*S
#ifdef DEBUG_GIVENS
	//Faust::MatDense<FPP,Cpu> L_copy = L;
	cout << "L(p,q) before update_L():" << L(this->p,this->q) << endl;
#endif
#ifdef NO_OPT_UPDATE_L_CPLX
	this->facts[this->ite].multiply(L, 'H');
	L.multiplyRight(this->facts[this->ite]);
//	facts[ite].multiply(L, 'N');
//	facts[ite].transpose();
//	facts[ite].conjugate();
//	L.multiplyRight(facts[ite]);
//	facts[ite].transpose();
//	facts[ite].conjugate();
#else
	Faust::Vect<FPP,DEVICE> L_vec_p, L_vec_q;
	FPP c_qq = *(this->fact_mod_values.end()-1);
	FPP c_qp = *(this->fact_mod_values.end()-2);
	FPP c_pq = *(this->fact_mod_values.end()-3);
	FPP c_pp = *(this->fact_mod_values.end()-4);
	update_L_first(L_vec_p, L_vec_q, c_pp, c_pq, c_qp, c_qq, this->p, this->q, L);
	update_L_second(L_vec_p, L_vec_q, c_pp, c_pq, c_qp, c_qq, this->p, this->q, L);
#endif
#ifdef DEBUG_GIVENS
	cout << "L(p,q) after update_L():" << L(p,q) << endl;
	cout << "L(p,p): "<< L(p,p) << endl;
	cout << "L(q,q): "<< L(q,q) << endl;
//	for(int i=0;i<L.getNbRow();i++)
//		for(int j=i+1;j<L.getNbCol();j++)
//		{
//			cout << Faust::fabs(L(i,j)-conj(L(j,i))) << endl;
//			assert(Faust::fabs(L(i,j)-conj(L(j,i))) < 1e-3);
//		}
//	cout << endl;
//#endif
//	if(Faust::fabs(L(p,q)) > 1e-3)
//	{
//		cout << "launch recomputing after not null L(p,q)" << endl;
//		L = L_copy;
//		this->theta2 = - this->theta2;
//		this->update_fact();
//		this->update_L(L);
//	}
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::update_L(Faust::MatSparse<FPP,Cpu> & L)
{
//#define OPT_UPDATE_SPARSE_L
#undef OPT_UPDATE_SPARSE_L
#ifndef OPT_UPDATE_SPARSE_L
	// L = S'*L*S
	this->facts[this->ite].multiply(L, 'H');
	L.multiplyRight(this->facts[this->ite]);
#else
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> L_vec_p, L_vec_q;
	FPP c_qq = *(this->fact_mod_values.end()-1);
	FPP c_qp = *(this->fact_mod_values.end()-2);
	FPP c_pq = *(this->fact_mod_values.end()-3);
	FPP c_pp = *(this->fact_mod_values.end()-4);

	update_L_first(L_vec_p, L_vec_q, c_pp, c_pq, c_qp, c_qq, this->p, this->q, L);
	L.update_dim();
//	update_L_second(L_vec_p, L_vec_q, c, s, this->p, this->q, L);
	L.multiplyRight(this->facts[this->ite]);
#endif
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::update_L()
{
	MatSparse<FPP, DEVICE>* sL = dynamic_cast<MatSparse<FPP, DEVICE>*>(this->L);
	if(sL) update_L(*sL);
	else update_L(*dynamic_cast<MatDense<FPP, DEVICE>*>(this->L));
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::update_L_first(Faust::Vect<FPP,DEVICE>& L_vec_p, Faust::Vect<FPP,DEVICE>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, Faust::MatDense<FPP,DEVICE> & L)
{
#define copy_vec2Lrow(vec,rowi) \
	for(int i=0;i<L.getNbCol();i++) L.getData()[L.getNbRow()*i+rowi] = tmp[i]

	Faust::Vect<FPP,DEVICE> tmp, tmp2;
	/*========== L = S'*L */
	L_vec_p = L.get_row(p), L_vec_q = L.get_row(q);
	// L(p,:) = c*L(p,:) + s*L(q,:)
	tmp = L_vec_p;
	tmp *= conj(c_pp);
	tmp2 = L_vec_q;
	tmp2 *= conj(c_qp);
	tmp += tmp2;
	copy_vec2Lrow(tmp,p);

	// L(q,:) = -s*L(p,:) + c*L(q,:)
	tmp = L_vec_p;
	tmp *= conj(c_pq);
	tmp2 = L_vec_q;
	tmp2 *= conj(c_qq);
	tmp += tmp2;
	copy_vec2Lrow(tmp, q);
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::update_L_second(Faust::Vect<FPP,DEVICE>& L_vec_p, Faust::Vect<FPP,DEVICE>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, Faust::MatDense<FPP,DEVICE> & L)
{
	Faust::Vect<FPP,DEVICE> tmp, tmp2;
	/*========== L *= S */
	L_vec_p = L.get_col(p), L_vec_q = L.get_col(q);
	// L(:,p) = c*L(:,p) + s*L(:,q)
	tmp = L_vec_p;
	tmp *= c_pp;
	tmp2 = L_vec_q;
	tmp2 *= c_qp;
	tmp += tmp2;

	memcpy(L.getData()+L.getNbRow()*p, tmp.getData(), sizeof(FPP)*L.getNbRow());
	// L(:,q) = -s*L(:,p) + c*L(:,q)
	tmp = L_vec_p;
	tmp *= c_pq;
	tmp2 = L_vec_q;
	tmp2 *= c_qq;
	tmp += tmp2;
	memcpy(L.getData()+L.getNbRow()*q, tmp.getData(), sizeof(FPP)*L.getNbRow());
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::update_L_second(Eigen::SparseMatrix<FPP,Eigen::RowMajor > & L_vec_p, Eigen::SparseMatrix<FPP, Eigen::RowMajor>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, Faust::MatSparse<FPP,DEVICE> & L)
{
	Eigen::SparseMatrix<FPP, Eigen::RowMajor> tmp, tmp2;
	/*========== L *= S */
	//L_vec_p = L.get_col(p), L_vec_q = L.get_col(q);
	L_vec_p = L.mat.block(0, p, L.getNbRow(), 1);
	L_vec_q = L.mat.block(0, q, L.getNbRow(), 1);
	// L(:,p) = c*L(:,p) + s*L(:,q)
	tmp = L_vec_p;
	tmp *= c_pp;
	tmp2 = L_vec_q;
	tmp2 *= c_qp;
	tmp += tmp2;
	L.mat.col(p) = tmp;

	// L(:,q) = -s*L(:,p) + c*L(:,q)
	tmp = L_vec_p;
	tmp *= c_pq;
	tmp2 = L_vec_q;
	tmp2 *= c_qq;
	tmp += tmp2;
	L.mat.col(q) = tmp;
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void EigTJComplex<FPP,DEVICE,FPP2>::update_L_first(Eigen::SparseMatrix<FPP, Eigen::RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, Eigen::RowMajor>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, Faust::MatSparse<FPP,DEVICE> & L)
{
	Eigen::SparseMatrix<FPP, Eigen::RowMajor> tmp, tmp2, tmp3;
	/*========== L = S'*L */
	L_vec_p = L.mat.innerVector(p), L_vec_q = L.mat.innerVector(q);
//	L_vec_p = L.mat.block(p, 0, 1, L.getNbCol());
//	L_vec_q = L.mat.block(q, 0, 1, L.getNbCol());
	// L(p,:) = c*L(p,:) + s*L(q,:)
	tmp = L_vec_p;
	tmp *= c_pp;
	tmp2 = L_vec_q;
	tmp2 *= c_qp;
	tmp += tmp2;
	L.mat.innerVector(p) = tmp;
	//
	//
	// L(q,:) = -s*L(p,:) + c*L(q,:)
	tmp = L_vec_p;
	tmp *= c_pq;
	tmp2 = L_vec_q;
	tmp2 *= c_qq;
	tmp += tmp2;
	L.mat.innerVector(q) = tmp;
}

template<typename FPP, FDevice DEVICE, typename FPP2>
EigTJComplex<FPP,DEVICE,FPP2>::EigTJComplex(Faust::MatSparse<FPP,DEVICE>& Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError /* default to 0.0 */, const bool errIsRel, const bool enable_large_Faust, const int err_period/*=100*/) : Faust::EigTJGen<typename FPP::value_type, DEVICE, FPP2, FPP>(Lap, J, verbosity, stoppingError, errIsRel, enable_large_Faust, err_period), C(Lap.getNbRow(), Lap.getNbCol())
{
	//see parent ctor
	this->C.setZeros();
}

template<typename FPP, FDevice DEVICE, typename FPP2>
EigTJComplex<FPP,DEVICE,FPP2>::EigTJComplex(Faust::MatDense<FPP,DEVICE>& Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust, const int err_period/*=100*/) : Faust::EigTJGen<typename FPP::value_type, DEVICE, FPP2, FPP>(Lap, J, verbosity, stoppingError, errIsRel, enable_large_Faust, err_period), C(Lap.getNbRow(), Lap.getNbCol())
{
	// see parent ctor
	this->C.setZeros();
}

template<typename FPP, FDevice DEVICE, typename FPP2>
const Faust::MatSparse<FPP,DEVICE> EigTJComplex<FPP,DEVICE,FPP2>::get_Dspm(const bool ord /* default to false*/)
{
	MatSparse <FPP,DEVICE> spD;
	EigTJGen<typename FPP::value_type,DEVICE,FPP2,FPP>::get_Dspm(spD, ord);
	return spD;
}


