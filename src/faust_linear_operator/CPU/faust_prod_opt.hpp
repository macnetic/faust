#include "faust_openmp.h"
//TODO: allow nflags == 1 and all factors using with same flag
/**
 *	\brief Multiply all the matrices together (as a product of n factors) with an optimization based on the associativity of the matrix product ; following an order that minimizes the number of scalar operations.
 *	\note: the std::vector facts is altered after the call! Don't reuse it.
 */
template<typename FPP, FDevice DEVICE>
void Faust::multiply_order_opt_all_ends(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	Faust::MatDense<FPP, DEVICE> tmpr, tmpl;
	int nfacts = facts.size();
	int ri = 0, li = nfacts-1;
	Faust::MatDense<FPP,DEVICE> *R1, *R2, *L1, *L2;
	faust_unsigned_int R1nr, R1nc, L1nr, L1nc;
	while(li-ri > 1)
	{
		R1 = facts[ri];
		R2 = facts[ri+1];
		L1 = facts[li-1];
		L2 = facts[li];
		R1nr = R1->getNbRow();
		R1nc = R1->getNbCol();
		L1nr = L1->getNbRow();
		L1nc = L1->getNbCol();
		if(R1nr * R1nc * R2->getNbCol() < L1nr * L1nc * L2->getNbCol())
		{
			gemm(*R1, *R2, tmpr, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>ri?ri:0], transconj_flags[transconj_flags.size()>ri+1?ri+1:0]);
			ri++;
			facts[ri] = &tmpr;
			if(transconj_flags.size() > ri)
				transconj_flags[ri] = 'N';
		}
		else
		{
			gemm(*L1, *L2, tmpl, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>li-1?li-1:0], transconj_flags[transconj_flags.size()>li?li:0]);
			li--;
			facts[li] = &tmpl;
			if(transconj_flags.size() > li)
				transconj_flags[li] = 'N';
		}
	}
	// last mul
	gemm(*facts[ri], *facts[li], out, alpha, beta_out, ri==0?transconj_flags[0]:'N', li==nfacts-1&&transconj_flags.size()>li?transconj_flags[li]:'N');
	facts.erase(facts.begin(), facts.end());
}

template<typename FPP, FDevice DEVICE>
void Faust::multiply_order_opt_all_ends(std::vector<Faust::MatGeneric<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	Faust::MatDense<FPP, DEVICE> tmpr, tmpl;
	int nfacts = facts.size();
	int ri = 0, li = nfacts-1;
	Faust::MatGeneric<FPP,DEVICE> *R1, *R2, *L1, *L2;
	faust_unsigned_int R1nr, R1nc, L1nr, L1nc;
	while(li-ri > 1)
	{
		R1 = facts[ri];
		R2 = facts[ri+1];
		L1 = facts[li-1];
		L2 = facts[li];
		R1nr = R1->getNbRow();
		R1nc = R1->getNbCol();
		L1nr = L1->getNbRow();
		L1nc = L1->getNbCol();
		if(R1nr * R1nc * R2->getNbCol() < L1nr * L1nc * L2->getNbCol())
		{
			gemm_gen(*R1, *R2, tmpr, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>ri?ri:0], transconj_flags[transconj_flags.size()>ri+1?ri+1:0]);
			ri++;
			facts[ri] = &tmpr;
			if(transconj_flags.size() > ri)
				transconj_flags[ri] = 'N';
		}
		else
		{
			gemm_gen(*L1, *L2, tmpl, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>li-1?li-1:0], transconj_flags[transconj_flags.size()>li?li:0]);
			li--;
			facts[li] = &tmpl;
			if(transconj_flags.size() > li)
				transconj_flags[li] = 'N';
		}
	}
	// last mul
	gemm_gen(*facts[ri], *facts[li], out, alpha, beta_out, ri==0?transconj_flags[0]:'N', li==nfacts-1&&transconj_flags.size()>li?transconj_flags[li]:'N');
	facts.erase(facts.begin(), facts.end());
}

template<typename FPP, FDevice DEVICE>
void Faust::multiply_order_opt_all_best(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	std::vector<Faust::MatDense<FPP,DEVICE>*> tmp_facts; //temporary product results
	Faust::MatDense<FPP, DEVICE>* tmp;
	int nfacts = facts.size();
	Faust::MatDense<FPP,DEVICE> *Si, *Sj;
	std::vector<int> complexity(nfacts-1);
	for(int i = 0; i <nfacts-1; i++)
	{
		Si = facts[i];
		Sj = facts[i+1];
		complexity[i] = Si->getNbRow() * Si->getNbCol() * Sj->getNbCol();
	}
	int idx; // marks the factor to update with a product of contiguous factors
	bool multiplying_tmp_factor = false; // allows to avoid to allocate uselessly a tmp factor if Si or Sj are already tmp factors
	while(facts.size() > 2)
	{
		// find the least complex product facts[idx]*facts[idx+1]
		idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));
		Si = facts[idx];
		Sj = facts[idx+1];
		for(auto Tit = tmp_facts.begin(); Tit != tmp_facts.end(); Tit++)
		{
			// TODO: both Sj and Si can be tmp_factor
			if(Sj == *Tit)
			{// Sj is original fact
				multiplying_tmp_factor = true;
				tmp = Sj;
				break;
			}
			else if(Si == *Tit)
			{
				multiplying_tmp_factor = true;
				tmp = Si;
				break;
			}
		}
		if(! multiplying_tmp_factor)
		{
			tmp = new Faust::MatDense<FPP, DEVICE>();
			tmp_facts.push_back(tmp);
		}
		//else no need to instantiate a new tmp, erasing Sj which is a tmp
		gemm(*Si, *Sj, *tmp, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>idx?idx:0], transconj_flags[transconj_flags.size()>idx+1?idx+1:0]);
		facts.erase(facts.begin()+idx+1);
		complexity.erase(complexity.begin()+idx); //complexity size == facts size - 1
		facts[idx] = tmp;
		if(transconj_flags.size() > idx)
			transconj_flags[idx] = 'N';
		// update complexity around the new factor

		if(facts.size() > 2)
		{
			if(idx > 0)
				complexity[idx-1] = facts[idx-1]->getNbRow() * facts[idx-1]->getNbCol() * facts[idx]->getNbCol(); //whether facts[idx-1] is transposed or not, it's the same
			if(idx < facts.size()-1)
				complexity[idx] = facts[idx]->getNbRow() * facts[idx]->getNbCol();
				if(transconj_flags.size() > idx && transconj_flags[idx+1] != 'N')
					complexity[idx] *= facts[idx+1]->getNbRow();
				else
					complexity[idx] *= facts[idx+1]->getNbCol();
		}
		multiplying_tmp_factor = false;
	}
	// last mul
	gemm(*facts[0], *facts[1], out, alpha, beta_out, transconj_flags[0], transconj_flags.size()>1?transconj_flags[1]:'N');
	facts.erase(facts.begin(), facts.end());
	// delete all tmp facts
	for(auto Tit = tmp_facts.begin(); Tit != tmp_facts.end(); Tit++)
	{
		delete *Tit;
	}
}

template<typename FPP, FDevice DEVICE>
void Faust::multiply_order_opt_all_best(std::vector<Faust::MatGeneric<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	Faust::MatSparse<FPP,Cpu> * tmp_sp1,* tmp_sp2;
	Faust::MatDense<FPP,Cpu> * tmp_ds1,* tmp_ds2;
	std::vector<Faust::MatGeneric<FPP,DEVICE>*> tmp_facts; //temporary product results
	Faust::MatDense<FPP, DEVICE>* tmp;
	int nfacts = facts.size();
	Faust::MatGeneric<FPP,DEVICE> *Si, *Sj;
	std::vector<int> complexity(nfacts-1);
	auto calc_cost = [](Faust::MatGeneric<FPP,DEVICE> *Si, Faust::MatGeneric<FPP,DEVICE> *Sj)
	{
		Faust::MatSparse<FPP,Cpu> * tmp_sp;
		int complexity;
		if((tmp_sp = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(Si)))
			complexity = Si->getNonZeros() * Sj->getNbCol();
		else // Si dense
			if((tmp_sp = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(Sj)))
				complexity = Si->getNbRow() * Sj->getNonZeros();
		else // Si and Sj dense
			complexity = Si->getNbRow() * Si->getNbCol() * Sj->getNbCol();
		return complexity;
	};
	std::string type_err = "Sj shouldn't be anything else than a MatSparse or MatDense.";
	auto reverse_transp_flag = [](char & flag)
	{
		switch(flag)
		{
			case 'T':
			case 'H':
				return 'N';
			case 'N':
				return 'H';
		}
	};
	for(int i = 0; i <nfacts-1; i++)
	{
		Si = facts[i];
		Sj = facts[i+1];
		complexity[i] = calc_cost(Si, Sj);
	}
	int idx; // marks the factor to update with a product of contiguous factors
	bool multiplying_tmp_factor = false; // allows to avoid to allocate uselessly a tmp factor if Si or Sj are already tmp factors
	while(facts.size() > 2)
	{
		// find the least complex product facts[idx]*facts[idx+1]
		idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));
		Si = facts[idx];
		Sj = facts[idx+1];
		tmp = new Faust::MatDense<FPP, DEVICE>();
		tmp_facts.push_back(tmp);
		gemm_gen(*Si, *Sj, *tmp, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>idx?idx:0], transconj_flags[transconj_flags.size()>idx+1?idx+1:0]);
		facts.erase(facts.begin()+idx+1);
		complexity.erase(complexity.begin()+idx); //complexity size == facts size - 1
		facts[idx] = tmp;
		if(transconj_flags.size() > idx)
		{
			if(transconj_flags.size() > idx+1)
				transconj_flags.erase(transconj_flags.begin()+idx+1);
			transconj_flags[idx] = 'N';
		}
		// update complexity around the new factor
		if(facts.size() > 2)
		{
			if(idx > 0)
				complexity[idx-1] = calc_cost(facts[idx-1], facts[idx]);
			if(idx < facts.size()-1)
				complexity[idx] = calc_cost(facts[idx], facts[idx+1]);
		}
		multiplying_tmp_factor = false;
	}
	// last mul
	gemm_gen(*facts[0], *facts[1], out, alpha, beta_out, transconj_flags[0], transconj_flags.size()>1?transconj_flags[1]:'N');
	facts.erase(facts.begin(), facts.end());
	// delete all tmp facts
	for(auto Tit = tmp_facts.begin(); Tit != tmp_facts.end(); Tit++)
	{
		delete *Tit;
	}
}

template<typename FPP, FDevice DEVICE>
void Faust::multiply_order_opt_first_best(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	int nfacts = facts.size();
	Faust::MatDense<FPP,DEVICE> *Si, *Sj, tmp;
	if(nfacts == 1)
	{
		tmp = *facts[0];
		if(transconj_flags[0] != 'N')
		{
			tmp.conjugate(false);
			tmp.transpose();
		}
		tmp *= alpha;
		if(beta_out != FPP(0))
		{
			out *= beta_out;
			out += tmp;
		}
		else
			out = tmp;

	}
	if(nfacts <= 2)
	{
		gemm(*facts[0], *facts[1], out, alpha, beta_out, transconj_flags[0], transconj_flags[transconj_flags.size()>1?1:0]);
		return;
	}
	std::vector<int> complexity(nfacts-1);
//	int min_cplx = std::numeric_limits<int>::max(); //doesn't work with VS14
	int min_cplx = 1 << 30;
	int i, idx = 0, lasti; // idx marks the factor to update with a product of contiguous factors
	for(int i = 0; i <nfacts-1; i++)
	{
		Si = facts[i];
		Sj = facts[i+1];
		complexity[i] = Si->getNbRow() * Si->getNbCol() * Sj->getNbCol();
		if(complexity[i] < min_cplx)
		{
			min_cplx = complexity[i];
			idx = i;
		}
	}
	// find the least complex product facts[idx]*facts[idx+1]
//	idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));
//	std::cout << "idx: " << idx << std::endl;
	Si = facts[idx];
	Sj = facts[idx+1];
	gemm(*Si, *Sj, tmp, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>idx?idx:0], transconj_flags[transconj_flags.size()>idx+1?idx+1:0]);
	lasti = idx+2<facts.size()?0:1;
	i = idx-1;
	while(i >= lasti)
	{
		// multiply on the left
		gemm(*facts[i], tmp, tmp, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>i?i:0], 'N');
		i--;
	}
	if(lasti == 0)
	{
		// now multiply on the right
		i = idx+2;
		while(i < facts.size()-1)
		{
			gemm(tmp, *facts[i], tmp, (FPP)1.0, (FPP)0, 'N', transconj_flags[transconj_flags.size()>i?i:0]);
			i++;
		}
		//last mul
		gemm(tmp, *facts[i], out, alpha, beta_out, 'N', transconj_flags[transconj_flags.size()>i?i:0]);
	}
	else
		gemm(*facts[i], tmp, out, alpha, beta_out, transconj_flags[transconj_flags.size()>i?i:0], 'N');
}

template<typename FPP, FDevice DEVICE>
void Faust::multiply_order_opt(const int mode, std::vector<Faust::MatGeneric<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	if(mode == GREEDY_ALL_BEST_GENMAT)
	{
		// no need to copy/convert for this method
		Faust::multiply_order_opt_all_best(facts, out, alpha, beta_out, transconj_flags);
		return;
	}
	std::vector<Faust::MatDense<FPP,DEVICE>*> dfacts;
	std::vector<Faust::MatDense<FPP,DEVICE>*> dfacts_to_del;
	Faust::MatDense<FPP,DEVICE> * df = nullptr;
	for(auto f: facts)
	{
		if(!(df = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(f)))
		{
			Faust::MatSparse<FPP,DEVICE> *sf = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(f);
			df = new Faust::MatDense<FPP,DEVICE>(*sf);
			dfacts_to_del.push_back(df);
		}
		dfacts.push_back(df);
	}
	switch(mode)
	{
		case GREEDY_ALL_ENDS:
			Faust::multiply_order_opt_all_ends(dfacts, out, alpha, beta_out, transconj_flags);
			break;
		case GREEDY_1ST_BEST:
			Faust::multiply_order_opt_first_best(dfacts, out, alpha, beta_out, transconj_flags);
			break;
		case GREEDY_ALL_BEST_CONVDENSE:
			Faust::multiply_order_opt_all_best(dfacts, out, alpha, beta_out, transconj_flags);
			break;
		default:
			throw std::out_of_range("unknown optimization method asked");
	}
	// free dense copies
	for(auto df: dfacts_to_del)
		delete df;
}

template<typename FPP, FDevice DEVICE>
Faust::MatDense<FPP,DEVICE> Faust::multiply_omp(const std::vector<Faust::MatGeneric<FPP,DEVICE>*> data, const Faust::MatDense<FPP,DEVICE> A, const char opThis)
{
	Faust::MatDense<FPP, DEVICE> *M = nullptr;
#ifdef OMP_ENABLED
				// until this this method is disabled at compilation unless we manually define the constant in CFLAGS (for example).
	int nth, start_nth, thid, num_per_th, data_size;
	Faust::MatDense<FPP,DEVICE>* mats[8];
	std::vector<Faust::MatGeneric<FPP, DEVICE> *> _data(data.size()+1);
	Faust::MatSparse<FPP, DEVICE> * sM;
	Faust::MatDense<FPP,DEVICE>* tmp; // (_data[end_id-1]);
//	std::mutex m;
	int i = 0;
	for(auto ptr: data)
	{
//		std::cout << "i=" << i << " _data[i]=" << ptr << endl;
		_data[i++] = ptr;
	}
	_data[i] = const_cast<Faust::MatDense<FPP,DEVICE>*>(&A);
#pragma omp parallel private(nth, thid, tmp, num_per_th, data_size, start_nth)
	{
		data_size = data.size()+1; // _data + 1
		start_nth = omp_get_num_threads()<=data_size>>1?omp_get_num_threads():data_size>>1;
		// assert start_nth is even ?
		nth	= start_nth;
		num_per_th = data_size;
		num_per_th = num_per_th%2?num_per_th-1:num_per_th;
		num_per_th /= nth;
		num_per_th = num_per_th<2?2:num_per_th;
//		std::cout << "omp th: " << nth << " num_per_th: " << num_per_th << std::endl;
		thid = omp_get_thread_num();
//		std::cout << "num_per_th: " << num_per_th << std::endl;
		while(num_per_th > 1)
		{
			int first_id, end_id, id;
			first_id = num_per_th*thid;
			if(nth == start_nth && thid == nth - 1)
				end_id = data_size;
			else
				end_id = num_per_th*(thid+1);
			if(first_id < data_size)
			{
//				{
//					std::lock_guard<std::mutex> guard(m);
//					std::cout << "thread id: " << thid << " first_id: " << first_id << " end_id: " << end_id << std::endl;
//				}
				id = 'N' == opThis? end_id-1: first_id;
				if(sM = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(_data[id]))
					tmp = new Faust::MatDense<FPP,DEVICE>(*sM);
				else
					tmp = new Faust::MatDense<FPP,DEVICE>(*(Faust::MatDense<FPP,DEVICE>*)(_data[id]));
				mats[thid] = tmp;
				//				std::cout << "thid=" << thid << "mats[thid]=" << mats[thid] << "tmp=" << tmp << endl;
				if(opThis == 'N')
					for(int i=end_id-2;i>=first_id; i--)
					{
						//						std::cout << "mul:" << _data[i] << endl;
						if(sM = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(_data[i]))
							sM->multiply(*mats[thid], opThis);
						else
							dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(_data[i])->multiply(*mats[thid], opThis);
					}
				else
					for(int i=first_id+1;i < end_id; i++)
						if(sM = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(_data[i]))
							sM->multiply(*mats[thid],opThis);
						else
							dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(_data[i])->multiply(*mats[thid], opThis);
				_data[thid] = mats[thid];
				//				mats[thid]->Display();
				//				std::cout << mats[thid]->norm() << endl;
				//				std::cout << "thid=" << thid << "mats[thid]=" << mats[thid] << "_data[thid]=" << data[thid] << endl;
				data_size = nth;
			}
			if(nth > 1)
				num_per_th = nth/(nth>>1);
			else
				num_per_th = 1;
			nth >>= 1;

#pragma omp barrier
		}
		//		if(thid == 0) M = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(_data[0]);
	}
	M = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(mats[0]);
//	std::cout << M << endl;
	MatDense<FPP,DEVICE> M_;
	if(! M)
	{
		M_ = Faust::MatDense<FPP,DEVICE>(*dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(mats[0]));
		M = &M_;
	}
	//TODO: delete other thread mats
	//TODO: return a ptr instead of a copy
#else
	throw std::runtime_error("It's not possible to call Faust::multiply_omp because the library hasn't been compiled with this function enabled.");
#endif
	return *M;
}

template<typename FPP, FDevice DEVICE>
Faust::MatDense<FPP,DEVICE> Faust::multiply_par(const std::vector<Faust::MatGeneric<FPP,DEVICE>*>& data, const Faust::MatDense<FPP,DEVICE> A, const char opThis)
{

#ifdef OMP_ENABLED // this function doesn't use OpenMP but C++ threads are implemented using POSIX threads on Linux gcc as OpenMP, so for FAµST we assume it's the same
	int nth = std::thread::hardware_concurrency(); // https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency
	int barrier_count = nth;
	Faust::MatDense<FPP, DEVICE> *M;
	std::vector<std::thread*> threads;
	std::vector<Faust::MatDense<FPP,DEVICE>*> mats(nth);
	std::vector<Faust::MatGeneric<FPP, DEVICE>*> _data(data.size()+1);
	std::mutex barrier_mutex;
	std::condition_variable barrier_cv;
	int data_size = data.size()+1; // counting A in addition to Transform factors
	int num_per_th = data_size;
	num_per_th /= nth;
	num_per_th = num_per_th<2?2:num_per_th;
	// recompute effective needed nth
	nth = data_size / num_per_th;
	nth = nth%2?nth-1:nth; // nth must be even
	nth = nth == 0?1:nth;
//	std::cout << "nth=" << nth << endl;
	barrier_count = nth;
	int i = 0;
	if(opThis != 'N')
	{
		_data[0] = const_cast<Faust::MatDense<FPP,DEVICE>*>(&A);
		i = data_size-1;
	}
	else
		_data[data_size-1] = const_cast<Faust::MatDense<FPP,DEVICE>*>(&A);
	for(auto ptr: data)
	{
//		std::cout << "fac i=" << i << ": " << ptr << endl;
		if(opThis == 'N')
			_data[i++] = ptr;
		else
			_data[i--] = ptr;
	}
//	std::cout << "A ptr: " << &A << endl;
//	std::cout << A.norm() << endl;
	std::thread *th;
	for(i=0; i < nth; i++)
	{
//		th = new std::thread(&Faust::Transform<FPP,DEVICE>::multiply_par_run, nth, i, num_per_th, data_size, opThis, _data, mats, threads);
		th = new std::thread(Faust::multiply_par_run<FPP, Cpu>, nth, i, num_per_th, data_size, opThis, std::ref(_data), std::ref(mats), std::ref(threads), std::ref(barrier_mutex), std::ref(barrier_cv), std::ref(barrier_count));
		threads.push_back(th);
	}
	for(auto t: threads)
		if(t->get_id() != std::this_thread::get_id() )
			t->join();

	M = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(mats[0]);
//	std::cout << "M:" << M << endl;
	MatDense<FPP,DEVICE> M_;
	if(! M)
	{
		M_ = Faust::MatDense<FPP,DEVICE>(*dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(mats[0]));
		M = &M_;
	}
	//TODO: delete other thread mats
	//TODO: return a ptr instead of a copy
	//TODO: delete threads
	return *M;
#else
	throw std::runtime_error("It's not possible to call Faust::multiply_par because the library hasn't been compiled with this function enabled.");
#endif
}

template<typename FPP, FDevice DEVICE>
void Faust::multiply_par_run(int nth, int thid, int num_per_th, int data_size, char opThis, std::vector<Faust::MatGeneric<FPP, DEVICE>*>& data, std::vector<Faust::MatDense<FPP,DEVICE>*>& mats, std::vector<std::thread*>& threads, std::mutex & barrier_mutex, std::condition_variable & barrier_cv, int & barrier_count)
{
	Faust::MatSparse<FPP, DEVICE> * sM;
	Faust::MatDense<FPP, DEVICE> *M;
	Faust::MatDense<FPP,DEVICE>* tmp = nullptr, *tmp2 = nullptr; // (data[end_id-1]);
//	std::cout << "num_per_th: " << num_per_th << std::endl;
	int start_nth = nth;
	while(num_per_th > 1)
	{
		if(thid < nth)
		{
			int first_id, end_id, id;
			first_id = num_per_th*thid;
			end_id = num_per_th*(thid+1);
			if(thid == start_nth-1 && end_id < data_size) 
				end_id = data_size;
			else
				end_id = end_id>data_size?data_size:end_id;
			id = 'N' == opThis? end_id-1: first_id;
//			std::cout << "thid=" << thid << " end orig adr.:" << data[id] << endl;
			if(tmp != nullptr) tmp2 = tmp; //keep track of previous tmp to delete it afterward
			if((sM = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(data[id])))
				tmp = new Faust::MatDense<FPP,DEVICE>(*sM);
			else
				tmp = new Faust::MatDense<FPP,DEVICE>(*(Faust::MatDense<FPP,DEVICE>*)(data[id]));
			mats[thid] = tmp;
//			std::cout << "thid=" << thid << "mats[thid]=" << mats[thid] << "tmp=" << tmp << endl;
//			std::cout << "tmp.norm()" << tmp->norm() << endl;
			if(opThis == 'N')
				for(int i=end_id-2;i>=first_id; i--)
				{
//					std::cout << "mul:" << data[i] << " thid: " << thid << endl;
					if((sM = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(data[i])))
						sM->multiply(*mats[thid], opThis);
					else
						dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(data[i])->multiply(*mats[thid], opThis);
				}
			else
				for(int i=first_id+1;i < end_id; i++)
					if((sM = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(data[i])))
						sM->multiply(*mats[thid],opThis);
					else
						dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(data[i])->multiply(*mats[thid], opThis);
			//				mats[thid]->Display();
			//				std::cout << mats[thid]->norm() << endl;
//			std::cout << "thid=" << thid << "mats[thid]=" << mats[thid] << "data[thid]=" << data[thid] << endl;
			if(tmp2 != nullptr) delete tmp2;
		}
#define barrier() \
		{ \
			std::unique_lock<std::mutex> l1(barrier_mutex);  /* mutex block */ \
			barrier_count--;\
/*			string ite_str = i>=0?string(" (ite = ") + to_string(i) + ") ": "";*/ \
			if(barrier_count > 0)\
			{\
				/*	 the threads wait here for the last one to finish its iteration i*/ \
/*				std::cout << "thread " << thid << ite_str  << " is waiting on barrier (" << barrier_count << ")." << std::endl;*/ \
				barrier_cv.wait(l1/*, [&barrier_count, &num_threads]{return barrier_count == num_threads;}*/);\
/*				std::cout << "thread " << thid << ite_str <<" continues." << std::endl;*/\
			}\
			else {\
				/* the last thread to finish the iteration i wakes the other to continue the execution */ \
				/*				std::cout << "thread " << thid << ite_str << " reinits barrier." << std::endl;*/\
				barrier_count = start_nth;\
				barrier_cv.notify_all();\
			}\
		}
		barrier();
		if(thid < nth)
		{
			data[thid] = mats[thid];
		}
		barrier();
		if(nth > 1)
		{
			num_per_th = 2;
		}
		else
			num_per_th = 1;

		nth >>= 1;

		//		for(auto t: threads)
		//			if(t->get_id() != std::this_thread::get_id() )
		//				t->join();
//		std::cout <<"num_per_th: " <<  num_per_th << " thid:" << thid << endl;


	}
}

namespace Faust
{
	template<typename FPP>
		MatGeneric<FPP, Cpu>* dynprog_multiply_rec(const std::vector<MatGeneric<FPP, Cpu>*>& factors, int** inds, int i, int j, const char op='N', const char last_op='N')
		{
			int j_minus_i = j-i;
			int p_nrows, p_ncols; // prod numbers rows, cols
			char op_left, op_right;
			if(j_minus_i == 0)
				return factors[i];
			else if(j_minus_i == 1)
			{
				op_left = op; // i < factors.size()-1 because j_minus_i == 1
				if(j < factors.size()-1)
					op_right = op;
				else // j == last factor
					op_right = last_op;
				if(op_left == 'N')
					p_nrows = factors[i]->getNbRow();

				else
					p_nrows = factors[i]->getNbCol();
				if(op_right == 'N')
					p_ncols = factors[j]->getNbCol();
				else
					p_ncols = factors[j]->getNbRow();
				auto prod = new MatDense<FPP, Cpu>(p_nrows, p_ncols);
//				std::cout << factors[i]->getNbRow() << "x" << factors[i]->getNbCol() << "*" << factors[j]->getNbRow() << "x" << factors[j]->getNbCol() << std::endl;
				gemm_gen(*factors[i], *factors[j], *prod, (FPP)1.0, (FPP)0.0, op_left, op_right);
				return prod;
			}
			else
			{
				auto k = inds[i][j];


				auto prod1 = dynprog_multiply_rec(factors, inds, i, k, op, last_op);
				auto prod2 = dynprog_multiply_rec(factors, inds, k+1, j, op, last_op);

				if(i == k)
					op_left = op;
				else // prod1 was computed here, it's not one of factors[i], so no need to take care of op
					op_left = 'N';
				if(k+1 == j)
					if(j == factors.size()-1)
						op_right = last_op;
					else
						op_right = op;
				else // prod2 was computed here, it's not one of factors[i], so no need to take care of op
					op_right = 'N';

				if(op_left == 'N')
					p_nrows = prod1->getNbRow();
				else
					p_nrows = prod1->getNbCol();

				if(op_right == 'N')
					p_ncols = prod2->getNbCol();
				else
					p_ncols = prod2->getNbRow();

				auto prod12 = new MatDense<FPP, Cpu>(p_nrows, p_ncols);
//				std::cout << prod1->getNbRow() << "x" << prod1->getNbCol() << "*" << prod2->getNbRow() << "x" << prod2->getNbCol() << std::endl;
				gemm_gen(*prod1, *prod2, *prod12, FPP(1.0), FPP(0.0), op_left, op_right);
				//TODO/ verify if a memory leak exists
				// delete matrices allocated in this function's recursive calls
				if(k-i > 1)
					delete prod1;
				if(j-k-1 > 1)
					delete prod2;
				return prod12;
			}
		}


	template<typename FPP>
		MatDense<FPP, Cpu> dynprog_multiply(std::vector<MatGeneric<FPP, Cpu>*>& factors, const char op/*='N'*/, const MatGeneric<FPP, Cpu>* A/*=nullptr*/)
		{
			//TODO: manage useless cases when factors.size is too small
			char last_op = op;
			if(A != nullptr)
			{
				factors.push_back(const_cast<MatGeneric<FPP, Cpu>*>(A)); // A won't be touched
				last_op = 'N';
			}
			const int n = factors.size();
			int** c = new int*[n]; // TODO: reduce the memory used because only the upper triangle of the array is used
			int** inds = new int*[n]; // TODO: idem
			int j, k, cost;
			for(int i=0;i<n;i++)
			{
				c[i] = new int[n];
				inds[i] = new int[n];
				c[i][i] = 0;
			}
			for(int p=1;p<n;p++)
				for(int i=0;i<n-p;i++)
				{
					j = p + i;
					k = i;
					c[i][j] = std::numeric_limits<int>::max();
					while(k < j)
					{
						cost = c[i][k] + c[k+1][j];
						if(j < factors.size()-1 || A == nullptr)
						{
							if(op == 'N')
								cost += factors[i]->getNbRow() * factors[k]->getNbCol() * factors[j]->getNbCol();
							else
								cost += factors[i]->getNbCol() * factors[k]->getNbRow() * factors[j]->getNbRow();
						}
						else
						{ // j == factors.size()-1 && A == factors[n-1]
							// last_op == 'N'
							if(op == 'N')
								cost += factors[i]->getNbRow() * factors[k]->getNbCol() * factors[j]->getNbCol();
							else
								cost += factors[i]->getNbCol() * factors[k]->getNbRow() * factors[j]->getNbCol();
						}
						if(cost < c[i][j])
						{
							c[i][j] = cost;
							inds[i][j] = k;
						}
						k++;
					}
				}
//			std::cout << "inds:" << std::endl;
//			for(int i=0;i<n;i++)
//			{
//				for(int j=i+1;j<n;j++)
//				{
//					std::cout << inds[i][j] << " ";
//				}
//				std::cout << std::endl;
//			}
			auto prod = dynamic_cast<MatDense<FPP, Cpu>*>(dynprog_multiply_rec(factors, inds, 0, n-1, op, last_op));
			for(int i=0;i<n;i++)
			{
				delete[] c[i];
				delete[] inds[i];
			}
			delete[] c;
			delete[] inds;
			auto M = std::move(*prod);
			delete prod;
			if(A != nullptr)
				factors.erase(factors.end()-1);
			return std::move(M);
		}
}
