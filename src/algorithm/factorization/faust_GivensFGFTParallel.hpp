using namespace Faust;

#include <cmath>

template<typename FPP, Device DEVICE, typename FPP2>
GivensFGFTParallel<FPP,DEVICE,FPP2>::GivensFGFTParallel(Faust::MatDense<FPP,DEVICE>& Lap, int J, int t) : GivensFGFT<FPP,DEVICE,FPP2>(Lap,J), t(t), fact_nrots(0)
{
	this->facts.resize(round(J/(float)t));
	this->always_theta2 = true;
	this->coord_choices.resize(0);
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallel<FPP,DEVICE,FPP2>::max_L()
{
	// 	  matlab ref
	// 	  L_low = abs(tril(L,-1));
	//    ind_nnz = find(L_low);
	//    [~,ind_sorted] = sort(L_low(ind_nnz),'descend');
	//    [Rvec, Svec] = ind2sub([n,n],ind_nnz(ind_sorted));
	//    RSmat = [Rvec, Svec];
	MatDense<FPP,DEVICE> L_low = this->L;
	L_low.abs();
	L_low =	L_low.lower_tri(false);
	fact_nz_inds = L_low.nonzeros_indices();

#ifdef DEBUG_GIVENS
	for(auto &p : fact_nz_inds)
		cout << "GivensFGFTParallel::max_L() before sort (" << p.first+1 <<  "," << p.second+1 << ") :" << L_low(p.first, p.second) << endl;
#endif

	//ENOTE: can't use sort(it,it, lambda) as vector because list doesn't provide random access it.
	fact_nz_inds.sort([&L_low](const pair<int,int> &a, const pair<int,int> &b)
			{
			return L_low(a.first, a.second) > L_low(b.first, b.second);
			});
#ifdef DEBUG_GIVENS
	for(auto &p : fact_nz_inds)
		cout << "GivensFGFTParallel::max_L() after sort (" << p.first+1 <<  "," << p.second+1 << ") :" << L_low(p.first, p.second) << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallel<FPP,DEVICE,FPP2>::update_fact_nz_inds()
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
		cout << "GivensFGFTParallel::update_fact_nz_inds() after purge (" << p.first+1 <<  "," << p.second+1 << ") :" << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallel<FPP,DEVICE,FPP2>::next_step()
{

	substep_fun substep[] = {
		&GivensFGFTParallel<FPP,DEVICE,FPP2>::max_L,
		&GivensFGFTParallel<FPP,DEVICE,FPP2>::loop_update_fact, //responsible to call choose_pivot() and calc_theta(), update_fact()
		&GivensFGFTParallel<FPP,DEVICE,FPP2>::update_L,
		&GivensFGFTParallel<FPP,DEVICE,FPP2>::update_D,
		&GivensFGFTParallel<FPP,DEVICE,FPP2>::update_err};

	for(int i=0;i<sizeof(substep)/sizeof(substep_fun);i++)
	{
#ifdef DEBUG_GIVENS
		cout << "GivensFGFTParallel ite=" << this->ite << " substep i=" << i << endl;
#endif
		(this->*substep[i])();
	}
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallel<FPP,DEVICE,FPP2>::loop_update_fact()
{
	fact_nrots = 0;
	int k = 0;
	while(k < t && k < fact_nz_inds.size())
	{
		choose_pivot();
		update_fact_nz_inds();
		this->calc_theta();
		this->update_fact();
		fact_nrots++;
		k++;
	}
	finish_fact();
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallel<FPP,DEVICE,FPP2>::choose_pivot()
{
	// fact_nz_inds is assumed to be sorted in max_L()
	// the first one in the list maps the biggest abs val in L
	pair<int,int> max_elt = *(fact_nz_inds.begin());
	this->p = max_elt.first;
	this->q = max_elt.second;
	this->coord_choices.push_back(pair<int,int>(this->p,this->q));
#ifdef DEBUG_GIVENS
	cout << "choose_pivot() p: " << this->p+1 << " q:" << this->q+1 << " " << "L(p,q): " << this->L(this->p,this->q) << " nrots: " << fact_nrots << endl;
#endif
}


template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallel<FPP,DEVICE,FPP2>::update_fact()
{
	if(fact_nrots == 0){
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
	this->fact_mod_values.push_back(cos(this->theta));
	// 2nd
	this->fact_mod_row_ids.push_back(this->p);
	this->fact_mod_col_ids.push_back(this->q);
	this->fact_mod_values.push_back(-sin(this->theta));
	// 3rd
	this->fact_mod_row_ids.push_back(this->q);
	this->fact_mod_col_ids.push_back(this->p);
	this->fact_mod_values.push_back(sin(this->theta));
	// 4th
	this->fact_mod_row_ids.push_back(this->q);
	this->fact_mod_col_ids.push_back(this->q);
	this->fact_mod_values.push_back(cos(this->theta));
}


template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFTParallel<FPP,DEVICE,FPP2>::finish_fact()
{
	int n = this->Lap.getNbRow();
	this->facts[this->ite] = MatSparse<FPP,DEVICE>(this->fact_mod_row_ids, this->fact_mod_col_ids, this->fact_mod_values, n, n);
#ifdef DEBUG_GIVENS
	cout << "GivensFGFTParallel::finish_fact() ite: " << this->ite << " fact norm: " << this->facts[this->ite].norm() << endl;
	this->facts[this->ite].Display();
#endif
}
