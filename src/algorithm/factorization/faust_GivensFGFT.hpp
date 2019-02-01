using namespace Faust;

#include <cmath>

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::next_step()
{

	substep_fun substep[] = {
		&GivensFGFT<FPP,DEVICE,FPP2>::max_L_into_C,
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
	C_min_row.min(&p);
	q = q_candidates[p];
	coord_choices[ite] = pair<int,int>(p,q);
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::max_L_into_C()
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
	int n = Lap.getNbRow(), r, s;
	if(!ite)
	{
		//first iteration
		for(r=0;r<n;r++)
			for(s=r+1;s<n;s++)
			{
				//TODO: really slow (element-by-element copy) see if we can use a max and directly copy L per block
				// MatDense's data is in column-major order
				C.getData()[s*n+r] = -Faust::fabs(L(r,s));
			}
		C_min_row = C.rowwise_min(q_candidates);
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
		int pq[2] = { p , q };
		int rid;
		for(int i=0;i<2;i++)
		{
			r = pq[i];
			for(int s=r+1;s<n;s++)
				C.getData()[s*n+r] = - Faust::fabs(L(r,s));
			C_min_row.getData()[r] = C.get_row(r).min(&rid);
			q_candidates[r] = rid;
		}
		for(int i=0;i<2;i++)
		{
			s = pq[i];
			for(r=0;r<s-1;r++)
			{
				C.getData()[s*n+r] = - Faust::fabs(L(r,s));
				if(C(r,s) < C_min_row[r])
				{
					C_min_row[r] = C(r,s);
					q_candidates[r] = s;
				}
				else if(q_candidates[r] == s)
				{
					C_min_row.getData()[r] = C.get_row(r).min(&rid);
					q_candidates[r] = rid;
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

#define calc_err(theta) L(p,q)*cos(2*theta) + 0.5*(L(q,q) - L(p,p))*sin(2*theta)

	FPP2 theta1, theta2, err_theta1, err_theta2;

	theta1 = atan2(L(q,q) - L(p,p),(2*L(p,q)))/2;
	theta2 = theta1 + M_PI_4; // from cmath
	err_theta1 = calc_err(theta1);
	err_theta2 = calc_err(theta2);
	if(err_theta1 < err_theta2)
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
	int n = Lap.getNbRow();
	// forget previous rotation coeffs
	// and keep identity part (n first coeffs)
	fact_mod_row_ids.resize(n);
	fact_mod_col_ids.resize(n);
	fact_mod_values.resize(n);
	// write new coeffs
	// 1st one
	fact_mod_row_ids.push_back(p);
	fact_mod_col_ids.push_back(p);
	fact_mod_values.push_back(cos(theta));
	// 2nd
	fact_mod_row_ids.push_back(p);
	fact_mod_col_ids.push_back(q);
	fact_mod_values.push_back(-sin(theta));
	// 3rd
	fact_mod_row_ids.push_back(q);
	fact_mod_col_ids.push_back(p);
	fact_mod_values.push_back(sin(theta));
	// 4th
	fact_mod_row_ids.push_back(q);
	fact_mod_col_ids.push_back(q);
	fact_mod_values.push_back(cos(theta));
	facts[ite] = MatSparse<FPP,DEVICE>(fact_mod_row_ids, fact_mod_col_ids, fact_mod_values, n, n);
	cout << "norm facts ite: " << ite << facts[ite].norm() << endl;
	facts[ite].Display();
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L()
{
	// L = S'*L*S
	facts[ite].multiply(L, 'T');
	L.multiplyRight(MatDense<FPP,Cpu>(facts[ite]));
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_D()
{
	// D = spdiag(diag(L))
	for(int i=0;i<L.getNbRow();i++)
		D.setCoeff(i,i, L(i,i));
#ifdef DEBUG_GIVENS
	D.Display();
	cout << "D fro. norm: " << D.norm() << endl;
#endif
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
	MatDense<FPP,Cpu> tmp = D;
	tmp -= Lap;
	FPP2 err = tmp.norm(), err_d;
//	err *= err;
	err_d = Lap.norm();
//	err_d *= err_d;
	err /= err_d;
	cout << "ite. i: "<< ite << " err.: " << err << endl;
	//TODO: boolean to display error or verbose mode
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::compute_facts()
{
	ite = 0;
	while(ite < facts.size())
	{
		next_step();
		ite++;
	}
}

template<typename FPP, Device DEVICE, typename FPP2>
GivensFGFT<FPP,DEVICE,FPP2>::GivensFGFT(Faust::MatDense<FPP,DEVICE>& Lap, int J) : Lap(Lap), facts(J), D(Lap.getNbRow(), Lap.getNbCol()), C(Lap.getNbRow(), Lap.getNbCol()), errs(J), coord_choices(J), L(Lap), q_candidates(new int[Lap.getNbCol()])
{
	/** Matlab ref. code:
	 *     facts = cell(1,J);
	 *     n=size(Lap,1);
	 *     L=Lap;
	 *     C = 15*ones(n);
	 *     err=zeros(1,J);
	 *     coord_choices = zeros(2,J);
	 *
	 */
	C.setOnes();
	C.scalarMultiply(15); // purely abitrary
	if(Lap.getNbCol() != Lap.getNbRow())
		handleError("Faust::GivensFGFT", "Laplacian must be a square matrix.");

	// init the identity part of the factor buffer model
	// allocate the mem. space for the 4 additional rotation part coeffs
	for(int i=0;i<Lap.getNbRow();i++)
	{
		fact_mod_values.push_back(FPP(1));
		fact_mod_col_ids.push_back(i);
		fact_mod_row_ids.push_back(i);
	}

	// init. D
	D.setZeros();
}

template<typename FPP, Device DEVICE, typename FPP2>
FPP2 GivensFGFT<FPP,DEVICE,FPP2>::get_err(int j) const
{
	if(j > 0 && j < errs.size())
		return errs[j];
	else
		throw out_of_range("GivensFGFT::get_err(j): j is out of range.");
}

template<typename FPP, Device DEVICE, typename FPP2>
const vector<FPP2>& GivensFGFT<FPP,DEVICE,FPP2>::get_errs() const
{
	return errs;
}

template<typename FPP, Device DEVICE, typename FPP2>
const MatDense<FPP,DEVICE> GivensFGFT<FPP,DEVICE,FPP2>::get_D() const
{
	return D;
}

template<typename FPP, Device DEVICE, typename FPP2>
const MatDense<FPP,DEVICE>& GivensFGFT<FPP,DEVICE,FPP2>::get_L() const
{
	return L;
}

template<typename FPP, Device DEVICE, typename FPP2>
const vector<pair<int,int>>& GivensFGFT<FPP,DEVICE,FPP2>::get_coord_choices() const
{
	return coord_choices;
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::get_coord_choice(int j, int& p, int& q) const
{
	if(j > 0 && j < coord_choices.size())
	{
		p = coord_choices[j].first;
		q = coord_choices[j].second;
	}
	else
		throw out_of_range("GivensFGFT::get_coord_choice(j,p,q): j is out of range.");
}

template<typename FPP, Device DEVICE, typename FPP2>
const MatDense<FPP,DEVICE>& GivensFGFT<FPP,DEVICE,FPP2>::get_Lap() const
{
	return Lap;
}

template<typename FPP, Device DEVICE, typename FPP2>
const vector<Faust::MatSparse<FPP,DEVICE>>& GivensFGFT<FPP,DEVICE,FPP2>::get_facts() const
{
	return facts;
}


