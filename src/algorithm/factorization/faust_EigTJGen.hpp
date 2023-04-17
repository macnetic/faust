using namespace Faust;

// handle Visual Studio's whim https://msdn.microsoft.com/en-us/library/4hwaceh6.aspx
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI_4
#define M_PI_4 (M_PI/4.0)
#endif


template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void EigTJGen<FPP,DEVICE,FPP2,FPP4>::update_D()
{
	// D = spdiag(diag(L))
	for(int i=0;i<D.size();i++)
		D.getData()[i] = complex<FPP>((*L)(i,i)).real();
	//ENOTE: the cast to complex is necessary when FPP3 is double/float
#ifdef DEBUG_GIVENS
	D.Display();
	cout << "D fro. norm: " << D.norm() << endl;
#endif
	is_D_ordered = false;
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void EigTJGen<FPP,DEVICE,FPP2,FPP4>::order_D()
{
	order_D(1);
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void EigTJGen<FPP,DEVICE,FPP2,FPP4>::order_D(const int order /* -1 for descending order, 1 for ascending order, 0 no order*/)
{
	ordered_D = Faust::Vect<FPP,DEVICE>(D.size());
	ord_indices.resize(0);
	for(int i=0;i<D.size();i++)
		ord_indices.push_back(i);
	if(order)
		sort(ord_indices.begin(), ord_indices.end(), [this, &order](int i, int j) {
				auto di = std::abs(D.getData()[i]);
				auto dj = std::abs(D.getData()[j]);
				return order>0? di < dj:(order <0?di > dj:0);
				});
	// order == 0, no order
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

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const vector<int>& EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_ord_indices(const int order/*=1*/)
{
	if(! is_D_ordered || D_order_dir != order)
		order_D(order);
	return ord_indices;
}


template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void EigTJGen<FPP,DEVICE,FPP2,FPP4>::compute_facts()
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
		if(stopping = (ite > 1 && stoppingCritIsError && errs.size() > 2 && errs[errs.size()-1]-errs[errs.size()-2] > FLT_EPSILON))
			/*if(verbosity>0) */cerr << "WARNING: the eigtj algorithm stopped because the last error is greater than the previous one." << endl;
		if(stopping || stoppingCritIsError && errs.size() > 0 && (*(errs.end()-1) - stoppingError ) < FLT_EPSILON)
		{
			facts.resize(ite);
			break;
		}
	}
	if(verbosity > 1)
	{
		std::cout << "EigTJGen::compute_facts() end" << std::endl;
		std::cout << "J: " << J << std::endl;
		std::cout << "tol: " << stoppingError << std::endl;
		std::cout << "stopcrit is error: " << stoppingCritIsError << std::endl;
		std::cout << "relErr: " << errIsRel << std::endl;
		std::cout << "err_period: " << err_period << std::endl;
		std::cout << "order: " << is_D_ordered << std::endl;
		std::cout << "enable_large_Faust: " << enable_large_Faust << std::endl;
	}
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
EigTJGen<FPP,DEVICE,FPP2,FPP4>::EigTJGen(MatGeneric<FPP4,DEVICE>* Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust /* deft to false */, const int err_period/*=100*/) :
facts(J>0?(J*4<Lap->getNbRow()*Lap->getNbRow()||enable_large_Faust?J:0):0 /* don't allocate if the complexity doesn't worth it and enable_large_Faust is false*/), q_candidates(new int[Lap->getNbRow()]), J(J), D(Lap->getNbRow()), errs(0), coord_choices(0), Lap(*Lap), dim_size(Lap->getNbRow()), Lap_squared_fro_norm(0), is_D_ordered(false), D_order_dir(0), verbosity(verbosity), stoppingCritIsError(stoppingError != 0.0), stoppingError(stoppingError), errIsRel(errIsRel), enable_large_Faust(enable_large_Faust), ite(0), err_period(err_period), tag("")
{
	if(Lap->getNbCol() != Lap->getNbRow())
		handleError("Faust::EigTJComplex", "Laplacian must be a square matrix.");

	if(this->J == 0 && ! this->stoppingCritIsError) handleError("EigTJ", "Either J or stoppingError must be > 0");
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

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
EigTJGen<FPP,DEVICE,FPP2,FPP4>::EigTJGen(Faust::MatSparse<FPP4, DEVICE> & Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust/* deft to false */, const int err_period/*=100*/) :  EigTJGen<FPP,DEVICE,FPP2,FPP4>(&Lap, J, verbosity, stoppingError, errIsRel, enable_large_Faust, err_period)
{
	L = new MatSparse<FPP4,DEVICE>(Lap);
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
EigTJGen<FPP,DEVICE,FPP2,FPP4>::EigTJGen(Faust::MatDense<FPP4, DEVICE> & Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust/* deft to false */, const int err_period/*=100*/) : EigTJGen<FPP,DEVICE,FPP2,FPP4>(&Lap, J, verbosity, stoppingError, errIsRel,  enable_large_Faust, err_period)
{
	L = new MatDense<FPP4,DEVICE>(Lap);
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
FPP2 EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_err(int j) const
{
	if(j > 0 && j < errs.size())
		return errs[j];
	else
		throw out_of_range("EigTJGen::get_err(j): j is out of range.");
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const vector<FPP2>& EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_errs() const
{
	return errs;
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const Faust::Vect<FPP,DEVICE>& EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_D(const bool ord /* default to false */)
{
	if(ord)
	{
		if(!is_D_ordered)
			order_D();
		return ordered_D;
	}
	return D;
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const Faust::Vect<FPP,DEVICE>& EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_D(const int ord /* default to false */)
{
	if(ord != 0)
	{
		if(!is_D_ordered || ord != D_order_dir)
			order_D(ord);
		return ordered_D;
	}
	return D;
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
template<typename FPP3>
void EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_Dspm(Faust::MatSparse<FPP3,DEVICE> & spD, const bool ord /* default to false */)
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

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_D(FPP* diag_data, const bool ord /* default to false */)
{
	get_D(diag_data, ord?1:0);
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_D(FPP* diag_data, const int ord /* default to false */)
{
	const Faust::Vect<FPP,DEVICE>& D_ = get_D(ord);
	const FPP* src_data_ptr = D_.getData();
	memcpy(diag_data, src_data_ptr, sizeof(FPP)*D_.size());
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const Faust::MatDense<FPP4,DEVICE> EigTJGen<FPP,DEVICE,FPP2,FPP4>::compute_fourier(const bool ord /* default to false */)
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


template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const Faust::MatGeneric<FPP,DEVICE>& EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_L() const
{
	return *L;
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const vector<pair<int,int>>& EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_coord_choices() const
{
	return coord_choices;
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_coord_choice(int j, int& p, int& q) const
{
	if(j > 0 && j < coord_choices.size())
	{
		p = coord_choices[j].first;
		q = coord_choices[j].second;
	}
	else
		throw out_of_range("EigTJGen::get_coord_choice(j,p,q): j is out of range.");
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const Faust::MatDense<FPP4,DEVICE>& EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_Lap() const
{
	return Lap;
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
const vector<Faust::MatSparse<FPP4,DEVICE>>& EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_facts() const
{
	return facts;
}

template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
Faust::Transform<FPP4,DEVICE> EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_transform(const bool ord)
{
	return get_transform(ord?1:0);
}


template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
Faust::Transform<FPP4,DEVICE> EigTJGen<FPP,DEVICE,FPP2,FPP4>::get_transform(const int ord, const bool copy/*=true*/, const int first_nfacts/*=-1*/)
{
	//TODO: facts should be a Transform or a TransformHelper to avoid the risk of memory leak in caller code when ord == true and copy == false
	if(facts.size() == 0)
		throw out_of_range("The transform is empty. The algorithm stopped before computing any factor.");

	std::vector<MatGeneric<FPP4, DEVICE>*> facts;
	int end_id;

	size_t nfacts;
	if(first_nfacts < 0)
		// we want only effectively computed facts
		nfacts = this->ite;
	else
		// we want the number of facts asked by the caller
		nfacts = first_nfacts;

	if (ord)
		// we process the last factor separately because the columns must be reordered
		end_id = nfacts - 1;
	else
		end_id = nfacts;

	for(int i=0; i < end_id; i++)
		facts.push_back(&this->facts[i]);

	auto t = Faust::Transform<FPP4,DEVICE>(facts, /* lambda */ FPP4(1.0), /* optimizedCopy */false, /* cloning_fact */ copy);

	if(! copy)
		t.disable_dtor(); // don't want to free facts when t is out of scope because we might still need them in this object

	if(ord)
	{
		// (reorder eigenvector transform according to ordered D)
		MatSparse<FPP4,DEVICE> & last_fact = *(this->facts.begin() + end_id);
		auto P = new MatSparse<FPP4, DEVICE>(last_fact.getNbCol(), last_fact.getNbCol()); //last_fact permutation matrix
		if(!is_D_ordered || ord != D_order_dir)
			order_D(ord);
		for(int i=0;i<ord_indices.size();i++)
			P->setCoeff(ord_indices[i],i, FPP4(1.0));
		//		P.set_orthogonal(true);
		//		facts.push_back(P); // we prefer to directly multiply the last factor by P
//		last_fact.multiplyRight(P);
		last_fact.multiply(*P, 'N');
		t.push_back(P, false, false, false, copy); //default to optimizedCopy=false, transpose=false, conjugate=false, copying=false, verify_dims_agree=true
		if(copy) delete P;
		// Warning: if copy == false the caller is responsible to free the P
	}

	return t;

}


template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
FPP2 EigTJGen<FPP,DEVICE,FPP2,FPP4>::calc_err()
{
	// please read: https://gitlab.inria.fr/faustgrp/faust/-/issues/316
	// for proofs that the function computes the squared absolute error
	// i.e. \| Lap -  get_transform() * D  get_transform()' \|_F^2
	// (resp. the squared relative error if this->errIsRel is true)
	FPP2 err = 0, err_d;
	for(int i=0;i<this->D.size();i++)
		err += /*Faust::fabs(*/this->D(i)*this->D(i)/*)*/;
	if(this->Lap_squared_fro_norm == FPP2(0))
	{
		err_d = Faust::fabs(this->Lap.norm());
		err_d = err_d*err_d;
		this->Lap_squared_fro_norm = err_d;
	}
	else
		err_d = Faust::fabs(this->Lap_squared_fro_norm);
	err = Faust::fabs(err_d - err);
	if(this->errIsRel) err /= err_d;
	return err;
}


template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
void EigTJGen<FPP, DEVICE, FPP2, FPP4>::update_err()
{
	if(!((this->ite+1) % this->err_period) && this->stoppingCritIsError || this->verbosity > 1)
	{
		auto err = this->calc_err();
		if(this->verbosity)
		{
			if(tag.size())
				std::cout << tag << " ";
			std::cout << "factor : "<< this->ite <<  ", " << ((this->errIsRel)?"relative ":"absolute ") << "err.: " << err;
			std::cout << std::endl;
		}
		this->errs.push_back(err);
	}
}
