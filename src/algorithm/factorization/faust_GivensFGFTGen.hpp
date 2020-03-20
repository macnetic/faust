using namespace Faust;

// handle Visual Studio's whim https://msdn.microsoft.com/en-us/library/4hwaceh6.aspx
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI_4
#define M_PI_4 (M_PI/4.0)
#endif


template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
void GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::update_D()
{
	// D = spdiag(diag(L))
	for(int i=0;i<D.size();i++)
		D.getData()[i] = complex<FPP>((*L)(i,i)).real();
	//ENOTE: the cast to complex is necessary when FPP3 is double/float
#ifdef DEBUG_GIVENS
	D.Display();
	cout << "D fro. norm: " << D.norm() << endl;
#endif
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
void GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::order_D()
{
	order_D(1);
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
void GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::order_D(const int order /* -1 for descending order, 1 for ascending order */)
{
	ordered_D = Faust::Vect<FPP,DEVICE>(D.size());
	ord_indices.resize(0);
	for(int i=0;i<D.size();i++)
		ord_indices.push_back(i);
	sort(ord_indices.begin(), ord_indices.end(), [this, &order](int i, int j) {
			return order>0?D.getData()[i] < D.getData()[j]:(order <0?D.getData()[i] > D.getData()[j]:0);
			});
	for(int i=0;i<ord_indices.size();i++)
		ordered_D.getData()[i] = D.getData()[ord_indices[i]];
	// compute inverse permutation to keep easily possible retrieving undefined order
	inv_ord_indices.resize(ord_indices.size());
	int j = 0, i = 0;
	while(j<ord_indices.size())
		if(ord_indices[i] == j)
		{
			inv_ord_indices[j++] = i;
			i = 0;
		}
		else
			i++;
	is_D_ordered = true;
	D_order_dir = order;
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const vector<int>& GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_ord_indices()
{
	if(! is_D_ordered)
		order_D();
	return ord_indices;
}


template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
void GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::compute_facts()
{
	is_D_ordered = false; // facts (re)computed then D must be reordered
	ite = 0;
	bool stopping = false;
	if(stopping = !enable_large_Faust && ! stoppingCritIsError && dim_size*dim_size <= J*4)
	{
		cerr << "WARNING: the eigtj algorithm stopped because the transform to be computed doesn't worth it according to its complexity (in space and time) relatively to the size of the matrix to decompose. Still, if you want to force the computation, please use the enable_large_Faust flag." << endl;
		facts.resize(0);
	}
	while(! stopping && (J == 0 || ite < facts.size())) // when J == 0 the stopping criterion is the error against Lap
	{
		next_step();
		ite++;
		if(stopping = (ite > 1 && stoppingCritIsError && errs.size() > 2 && errs[ite-1]-errs[ite-2] > FLT_EPSILON))
			/*if(verbosity>0) */cerr << "WARNING: the eigtj algorithm stopped because the last error is greater than the previous one." << endl;
		if(stopping || stoppingCritIsError && errs.size() > 0 && (*(errs.end()-1) - stoppingError ) < FLT_EPSILON)
		{
			facts.resize(ite);
			break;
		}
	}
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::GivensFGFTGen(MatGeneric<FPP4,DEVICE>* Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust /* deft to false */) :
Lap(*Lap),  D(Lap->getNbRow()), errs(0), coord_choices(0), q_candidates(new int[Lap->getNbRow()]), is_D_ordered(false), verbosity(verbosity), stoppingCritIsError(stoppingError != 0.0), stoppingError(stoppingError), errIsRel(errIsRel), Lap_squared_fro_norm(0), facts(J>0?(J*4<Lap->getNbRow()*Lap->getNbRow()||enable_large_Faust?J:0):1 /* don't allocate if the complexity doesn't worth it and enable_large_Faust is false*/), last_fact_permuted(false), J(J), dim_size(Lap->getNbRow()), enable_large_Faust(enable_large_Faust)
{
	if(Lap->getNbCol() != Lap->getNbRow())
		handleError("Faust::GivensFGFTComplex", "Laplacian must be a square matrix.");

	if(this->J == 0 && ! this->stoppingCritIsError) handleError("GivensFGFT", "Either J or stoppingError must be > 0");
	// init the identity part of the factor buffer model
	// allocate the mem. space for the 4 additional rotation part coeffs
	for(int i=0;i<dim_size;i++)
	{
		fact_mod_values.push_back(FPP(1));
		fact_mod_col_ids.push_back(i);
		fact_mod_row_ids.push_back(i);
	}
	//init. D
	memset(this->D.getData(), 0, sizeof(FPP)*dim_size);
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::GivensFGFTGen(Faust::MatSparse<FPP4, DEVICE> & Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust/* deft to false */) :  GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>(&Lap, J, verbosity, stoppingError, errIsRel, enable_large_Faust)
{
	L = new MatSparse<FPP4,DEVICE>(Lap);
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::GivensFGFTGen(Faust::MatDense<FPP4, DEVICE> & Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust/* deft to false */) : GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>(&Lap, J, verbosity, stoppingError, errIsRel,  enable_large_Faust)
{
	L = new MatDense<FPP4,DEVICE>(Lap);
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
FPP2 GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_err(int j) const
{
	if(j > 0 && j < errs.size())
		return errs[j];
	else
		throw out_of_range("GivensFGFTGen::get_err(j): j is out of range.");
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const vector<FPP2>& GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_errs() const
{
	return errs;
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const Faust::Vect<FPP,DEVICE>& GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_D(const bool ord /* default to false */)
{
	if(ord)
	{
		if(!is_D_ordered)
			order_D();
		return ordered_D;
	}
	return D;
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const Faust::Vect<FPP,DEVICE>& GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_D(const int ord /* default to false */)
{
	if(ord != 0)
	{
		if(!is_D_ordered || ord != D_order_dir)
			order_D(ord);
		return ordered_D;
	}
	return D;
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
template<typename FPP3>
void GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_Dspm(Faust::MatSparse<FPP3,DEVICE> & spD, const bool ord /* default to false */)
{
	const Faust::Vect<FPP,DEVICE>& D_ = this->get_D(ord);
	vector<int> nat_ord_indices;
	vector<FPP3> diag_D;
	for(int i=0;i<D_.size();i++)
	{
		nat_ord_indices.push_back(i);
		diag_D.push_back(D_.getData()[i]);
	}
	spD = MatSparse<FPP3,DEVICE>(nat_ord_indices, nat_ord_indices, diag_D, nat_ord_indices.size(), nat_ord_indices.size());
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
void GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_D(FPP* diag_data, const bool ord /* default to false */)
{
	get_D(diag_data, ord?1:0);
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
void GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_D(FPP* diag_data, const int ord /* default to false */)
{
	const Faust::Vect<FPP,DEVICE>& D_ = get_D(ord);
	const FPP* src_data_ptr = D_.getData();
	memcpy(diag_data, src_data_ptr, sizeof(FPP)*D_.size());
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const Faust::MatDense<FPP4,DEVICE> GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::compute_fourier(const bool ord /* default to false */)
{
	Faust::MatDense<FPP4,Cpu> fourier(L->getNbRow(), L->getNbCol());
	Faust::MatDense<FPP4,Cpu>* ord_fourier;
	fourier.setEyes();
	for(int i=facts.size()-1;i>=0;i--)
		facts[i].multiply(fourier, 'N');
	if(ord)
	{
		if(!is_D_ordered)
			order_D();
		ord_fourier = fourier.get_cols(ord_indices);
		fourier = *ord_fourier;
		delete ord_fourier;
	}
	return fourier;
}


template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const Faust::MatGeneric<FPP,DEVICE>& GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_L() const
{
	return *L;
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const vector<pair<int,int>>& GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_coord_choices() const
{
	return coord_choices;
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
void GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_coord_choice(int j, int& p, int& q) const
{
	if(j > 0 && j < coord_choices.size())
	{
		p = coord_choices[j].first;
		q = coord_choices[j].second;
	}
	else
		throw out_of_range("GivensFGFTGen::get_coord_choice(j,p,q): j is out of range.");
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const Faust::MatDense<FPP4,DEVICE>& GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_Lap() const
{
	return Lap;
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
const vector<Faust::MatSparse<FPP4,DEVICE>>& GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_facts() const
{
	return facts;
}

template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
Faust::Transform<FPP4,DEVICE> GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_transform(bool ord)
{
	return get_transform(ord?1:0);
}


template<typename FPP, Device DEVICE, typename FPP2, typename FPP4>
Faust::Transform<FPP4,DEVICE> GivensFGFTGen<FPP,DEVICE,FPP2,FPP4>::get_transform(int ord)
{
	//TODO: an optimization is possible by changing type of facts to vector<MatGeneric*> it would avoid copying facts into Transform and rather use them directly. It will need a destructor that deletes them eventually if they weren't transfered to a Transform object before.
	if(facts.size() == 0)
		throw out_of_range("The transform is empty. The algorithm stopped before computing any factor.");
	MatSparse<FPP4,DEVICE> & last_fact = *(facts.end()-1);
	MatSparse<FPP4,DEVICE> P(last_fact.getNbCol(), last_fact.getNbCol()); //last_fact permutation matrix
	// (to reorder eigenvector transform according to ordered D)

	if(last_fact_permuted)
	{
		// get back to undefined order
		for(int i=0;i<inv_ord_indices.size();i++)
			P.setCoeff(inv_ord_indices[i],i, FPP4(1.0));
		last_fact.multiplyRight(P);
		last_fact_permuted = false;
	}

	if(ord)
	{
		if(!is_D_ordered || ord != D_order_dir)
			order_D(ord);
		for(int i=0;i<ord_indices.size();i++)
			P.setCoeff(ord_indices[i],i, FPP4(1.0));
		//		P.set_orthogonal(true);
		//		facts.push_back(P); // we prefer to directly multiply the last factor by P
		last_fact_permuted = true;
		last_fact.multiplyRight(P);
	}
	Faust::Transform<FPP4,DEVICE> t = Faust::Transform<FPP4,DEVICE>(facts);
	//	// remove the permutation factor if added temporarily for reordering
	//	return ord?facts.erase(facts.end()-1),t:t;
	return t;
}
