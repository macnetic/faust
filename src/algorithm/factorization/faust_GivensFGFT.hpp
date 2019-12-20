using namespace Faust;

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::next_step()
{

	substep_fun substep[] = {
		&GivensFGFT<FPP,DEVICE,FPP2>::max_L,
		&GivensFGFT<FPP,DEVICE,FPP2>::choose_pivot,
		&GivensFGFT<FPP,DEVICE,FPP2>::calc_theta,
		&GivensFGFT<FPP,DEVICE,FPP2>::update_fact,
		&GivensFGFT<FPP,DEVICE,FPP2>::update_L,
		&GivensFGFT<FPP,DEVICE,FPP2>::update_D,
		&GivensFGFT<FPP,DEVICE,FPP2>::update_err};

	for(int i=0;i<sizeof(substep)/sizeof(substep_fun);i++)
	{
#ifdef DEBUG_GIVENS
		cout << "GivensFGFT ite=" << ite << " substep i=" << i << endl;
#endif
		(this->*substep[i])();
	}
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::choose_pivot()
{
//Matlab ref. code:
//        [~,p] = min(C_min_row);
//        q = q_candidates(p);
//        coord_choices(1,j) = p;
//        coord_choices(2,j) = q;
	C_min_row.min_coeff(&this->p);
	this->q = this->q_candidates[this->p];
	this->coord_choices.push_back(pair<int,int>(this->p,this->q));
#ifdef DEBUG_GIVENS
	cout << "GivensFGFT::choose_pivot() ite: " << ite << " p: " << this->p << " q: " << this->q << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::max_L()
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
	int n = this->Lap.getNbRow(), r, s;
	if(!this->ite)
	{
		//first iteration
		for(r=0;r<n;r++)
			for(s=r+1;s<n;s++)
			{
				//TODO: really slow (element-by-element copy) see if we can use a max and directly copy L per block
				// MatDense's data is in column-major order
				C.getData()[s*n+r] = -Faust::fabs((*this->L)(r,s));
			}
		C_min_row = C.rowwise_min(this->q_candidates);
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
				C.getData()[s*n+r] = - Faust::fabs((*this->L)(r,s));
			C_min_row.getData()[r] = C.get_row(r).min_coeff(&rid);
			this->q_candidates[r] = rid;
		}
		for(int i=0;i<2;i++)
		{
			s = pq[i];
			for(r=0;r<s-1;r++)
			{
				C.getData()[s*n+r] = - Faust::fabs((*this->L)(r,s));
				if(C(r,s) < C_min_row[r])
				{
					C_min_row[r] = C(r,s);
					this->q_candidates[r] = s;
				}
				else if(this->q_candidates[r] == s)
				{
					C_min_row.getData()[r] = C.get_row(r).min_coeff(&rid);
					this->q_candidates[r] = rid;
				}
			}
		}
	}
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::calc_theta()
{
// Matlab ref. code:
//		theta1 = atan2(L(q,q) - L(p,p),2*L(p,q))/2 ;
//        err_theta1 = (L(p,q)*cos(2*theta1) + 0.5*(L(q,q) - L(p,p))*sin(2*theta1))^2;
//        theta2 = atan2(L(q,q) - L(p,p),2*L(p,q))/2 + pi/4;
//        err_theta2 = (L(p,q)*cos(2*theta2) + 0.5*(L(q,q) - L(p,p))*sin(2*theta2))^2;
//        if err_theta1<err_theta2
//            theta=theta1;
//        else
//            theta=theta2;
//        end
//

#define calc_err(theta) (*this->L)(this->p,this->q)*cos(2*theta) + 0.5*((*this->L)(this->q,this->q) - (*this->L)(this->p,this->p))*sin(2*theta)

	FPP2 theta1, theta2, err_theta1, err_theta2;

	theta1 = atan2((*this->L)(this->q,this->q) - (*this->L)(this->p,this->p),(2*(*this->L)(this->p,this->q)))/2;
	theta2 = theta1 + M_PI_4; // from cmath
	err_theta1 = calc_err(theta1);
	err_theta2 = calc_err(theta2);
	if(err_theta1 < err_theta2 && !always_theta2)
		theta = theta1;
	else
		theta = theta2;
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_fact()
{
	// Matlab ref. code:
	//        S = eye(n);
	//        S(p,p) = cos(theta); S(p,q) = -sin(theta);
	//        S(q,p) = sin(theta); S(q,q) = cos(theta);
	//        S = sparse(S);
	//        facts{j} = S;
	//
	int n = this->dim_size;
	FPP2 c = cos(theta);
	FPP2 s = sin(theta);
	// forget previous rotation coeffs
	// and keep identity part (n first coeffs)
	this->fact_mod_row_ids.resize(n);
	this->fact_mod_col_ids.resize(n);
	this->fact_mod_values.resize(n);
	// write new coeffs
	// 1st one
	this->fact_mod_row_ids.push_back(this->p);
	this->fact_mod_col_ids.push_back(this->p);
	this->fact_mod_values.push_back(c);
	// 2nd
	this->fact_mod_row_ids.push_back(this->p);
	this->fact_mod_col_ids.push_back(this->q);
	this->fact_mod_values.push_back(-s);
	// 3rd
	this->fact_mod_row_ids.push_back(this->q);
	this->fact_mod_col_ids.push_back(this->p);
	this->fact_mod_values.push_back(s);
	// 4th
	this->fact_mod_row_ids.push_back(this->q);
	this->fact_mod_col_ids.push_back(this->q);
	this->fact_mod_values.push_back(c);
	if(this->J == 0) this->facts.resize(this->ite+1);
	this->facts[this->ite] = MatSparse<FPP,DEVICE>(this->fact_mod_row_ids, this->fact_mod_col_ids, this->fact_mod_values, n, n);
	this->facts[this->ite].set_orthogonal(true);
#ifdef DEBUG_GIVENS
	cout << "GivensFGFT::update_fact() ite: " << ite << " fact norm: " << facts[this->ite].norm() << endl;
	facts[this->ite].Display();
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L(Faust::MatDense<FPP,Cpu> & L)
{
	// L = S'*L*S
#ifdef DEBUG_GIVENS
	cout << "L(p,q) before update_L():" << *this->L(p,q) << endl;
#endif
#define OPT_UPDATE_L
#ifndef OPT_UPDATE_L
	this->facts[this->ite].multiply(L, 'T');
	L.multiplyRight(this->facts[this->ite]);
#else
	Faust::Vect<FPP,DEVICE> L_vec_p, L_vec_q;
	FPP2 c = *(this->fact_mod_values.end()-1); // cos(theta)
	FPP2 s = *(this->fact_mod_values.end()-2); // sin(theta)
	update_L_first(L_vec_p, L_vec_q, c, s, this->p, this->q, L);
	update_L_second(L_vec_p, L_vec_q, c, s, this->p, this->q, L);
#endif
#ifdef DEBUG_GIVENS
	cout << "this->L(p,q) after update_L():" << L(p,q) << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L(Faust::MatSparse<FPP,Cpu> & L)
{
//#define OPT_UPDATE_SPARSE_L
#undef OPT_UPDATE_SPARSE_L
#ifndef OPT_UPDATE_SPARSE_L
	// L = S'*L*S
	this->facts[this->ite].multiply(L, 'T');
	L.multiplyRight(this->facts[this->ite]);
#else
	Eigen::SparseMatrix<FPP,RowMajor> L_vec_p, L_vec_q;
	FPP2 c = *(this->fact_mod_values.end()-1); // cos(theta)
	FPP2 s = *(this->fact_mod_values.end()-2); // sin(theta)
	update_L_first(L_vec_p, L_vec_q, c, s, this->p, this->q, L);
	L.update_dim();
//	update_L_second(L_vec_p, L_vec_q, c, s, this->p, this->q, L);
	L.multiplyRight(this->facts[ite]);
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L()
{
	MatSparse<FPP, DEVICE>* sL = dynamic_cast<MatSparse<FPP, DEVICE>*>(this->L);
	if(sL) update_L(*sL);
	else update_L(*dynamic_cast<MatDense<FPP, DEVICE>*>(this->L));
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L_first(Faust::Vect<FPP,DEVICE>& L_vec_p, Faust::Vect<FPP,DEVICE>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, Faust::MatDense<FPP,DEVICE> & L)
{
#define copy_vec2Lrow(vec,rowi) \
	for(int i=0;i<L.getNbCol();i++) L.getData()[L.getNbRow()*i+rowi] = tmp[i]

	Faust::Vect<FPP,DEVICE> tmp, tmp2;
	/*========== L = S'*L */
	L_vec_p = L.get_row(p), L_vec_q = L.get_row(q);
	// L(p,:) = c*L(p,:) + s*L(q,:)
	tmp = L_vec_p;
	tmp *= c;
	tmp2 = L_vec_q;
	tmp2 *= s;
	tmp += tmp2;
	copy_vec2Lrow(tmp,p);

	// L(q,:) = -s*L(p,:) + c*L(q,:)
	tmp = L_vec_p;
	tmp *= -s;
	tmp2 = L_vec_q;
	tmp2 *= c;
	tmp += tmp2;
	copy_vec2Lrow(tmp, q);
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L_second(Faust::Vect<FPP,DEVICE>& L_vec_p, Faust::Vect<FPP,DEVICE>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, Faust::MatDense<FPP,DEVICE> & L)
{
	Faust::Vect<FPP,DEVICE> tmp, tmp2;
	/*========== L *= S */
	L_vec_p = L.get_col(p), L_vec_q = L.get_col(q);
	// L(:,p) = c*L(:,p) + s*L(:,q)
	tmp = L_vec_p;
	tmp *= c;
	tmp2 = L_vec_q;
	tmp2 *= s;
	tmp += tmp2;

	memcpy(L.getData()+L.getNbRow()*p, tmp.getData(), sizeof(FPP)*L.getNbRow());
	// L(:,q) = -s*L(:,p) + c*L(:,q)
	tmp = L_vec_p;
	tmp *= -s;
	tmp2 = L_vec_q;
	tmp2 *= c;
	tmp += tmp2;
	memcpy(L.getData()+L.getNbRow()*q, tmp.getData(), sizeof(FPP)*L.getNbRow());
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L_second(Eigen::SparseMatrix<FPP,RowMajor > & L_vec_p, Eigen::SparseMatrix<FPP, RowMajor>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, Faust::MatSparse<FPP,DEVICE> & L)
{
	Eigen::SparseMatrix<FPP, RowMajor> tmp, tmp2;
	/*========== L *= S */
	//L_vec_p = L.get_col(p), L_vec_q = L.get_col(q);
	L_vec_p = L.mat.block(0, p, L.getNbRow(), 1);
	L_vec_q = L.mat.block(0, q, L.getNbRow(), 1);
	// L(:,p) = c*L(:,p) + s*L(:,q)
	tmp = L_vec_p;
	tmp *= c;
	tmp2 = L_vec_q;
	tmp2 *= s;
	tmp += tmp2;
	L.mat.col(p) = tmp;

	// L(:,q) = -s*L(:,p) + c*L(:,q)
	tmp = L_vec_p;
	tmp *= -s;
	tmp2 = L_vec_q;
	tmp2 *= c;
	tmp += tmp2;
	L.mat.col(q) = tmp;
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L_first(Eigen::SparseMatrix<FPP, RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, RowMajor>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, Faust::MatSparse<FPP,DEVICE> & L)
{
	Eigen::SparseMatrix<FPP, RowMajor> tmp, tmp2, tmp3;
	/*========== L = S'*L */
	L_vec_p = L.mat.innerVector(p), L_vec_q = L.mat.innerVector(q);
//	L_vec_p = L.mat.block(p, 0, 1, L.getNbCol());
//	L_vec_q = L.mat.block(q, 0, 1, L.getNbCol());
	// L(p,:) = c*L(p,:) + s*L(q,:)
	tmp = L_vec_p;
	tmp *= c;
	tmp2 = L_vec_q;
	tmp2 *= s;
	tmp += tmp2;
	L.mat.innerVector(p) = tmp;
	//
	//
	// L(q,:) = -s*L(p,:) + c*L(q,:)
	tmp = L_vec_p;
	tmp *= -s;
	tmp2 = L_vec_q;
	tmp2 *= c;
	tmp += tmp2;
	L.mat.innerVector(q) = tmp;
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_err()
{
	// Matlab ref. code:
	//        if mod(j,100)==0
	//            %err(j) = norm(D-L,'fro')^2/norm(L,'fro')^2;
	//            err(j) = norm(D-L,'fro')^2/norm(Lap,'fro')^2;
	//            %err(j) = norm(D-L)/norm(Lap);
	//            disp(['Iter ' num2str(j) ', error = ' num2str(err(j))])
	//            %    disp(['Number of edges: ' num2str(N_edges)])
	//        end
	//
	if(!((this->ite+1)%GivensFGFT<FPP,DEVICE,FPP2>::ERROR_CALC_PERIOD) || this->stoppingCritIsError || this->verbosity > 1)
	{
		FPP2 err = 0, err_d;
		for(int i=0;i<this->D.size();i++)
			err += this->D(i)*this->D(i);
		if(this->Lap_squared_fro_norm == FPP(0))
		{
			err_d = Faust::fabs(this->Lap.norm());
			err_d = err_d*err_d;
			this->Lap_squared_fro_norm = err_d;
		}
		else
			err_d = Faust::fabs(this->Lap_squared_fro_norm);
		err = Faust::fabs(err_d - err);
		if(this->errIsRel) err /= err_d;
		if(this->verbosity)
		{
			cout << "factor : "<< this->ite <<  ", " << ((this->errIsRel)?"relative ":"absolute ") << "err.: " << err;
			if(this->stoppingCritIsError) cout << " stoppingError: " << this->stoppingError << ")";
			cout << endl;
		}
		this->errs.push_back(err);
	}
}

template<typename FPP, Device DEVICE, typename FPP2>
GivensFGFT<FPP,DEVICE,FPP2>::GivensFGFT(Faust::MatSparse<FPP,DEVICE>& Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError /* default to 0.0 */, const bool errIsRel, const bool enable_large_Faust/* deft to false */) :  GivensFGFTGen<FPP,DEVICE,FPP2>(Lap, J, verbosity, stoppingError, errIsRel, enable_large_Faust), C(Lap.getNbRow(), Lap.getNbCol()), always_theta2(false)
{
	//see parent ctor
	C.setOnes();
	C.scalarMultiply(15); // purely abitrary
}

template<typename FPP, Device DEVICE, typename FPP2>
GivensFGFT<FPP,DEVICE,FPP2>::GivensFGFT(Faust::MatDense<FPP,DEVICE>& Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust/* deft to false */) : GivensFGFTGen<FPP,DEVICE,FPP2>(Lap, J, verbosity, stoppingError, errIsRel, enable_large_Faust), C(Lap.getNbRow(), Lap.getNbCol()), always_theta2(false)
{
	// see parent ctor
	C.setOnes();
	C.scalarMultiply(15); // purely abitrary
}

template<typename FPP, Device DEVICE, typename FPP2>
const Faust::MatSparse<FPP,DEVICE> GivensFGFT<FPP,DEVICE,FPP2>::get_Dspm(const bool ord /* default to false */)
{
	MatSparse <FPP,DEVICE> spD;
	GivensFGFTGen<FPP,DEVICE,FPP2>::get_Dspm(spD, ord);
	return spD;
}

