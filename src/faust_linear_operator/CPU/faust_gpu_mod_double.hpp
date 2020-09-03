#include "faust_TransformHelper.h"
namespace Faust
{
	template <>
	void FaustGPU<double>::load_gm_functions()
	{
		gm_DenseMatFunc_double* dsm_funcs;
		gm_MatArrayFunc_double* marr_funcs;
		gm_GenPurposeFunc_double* gp_funcs;

		if(FaustGPU<double>::marr_funcs == nullptr)
		{
			marr_funcs = new gm_MatArrayFunc_double(); // on the heap because because it cannot be shared among FaustGPU instances if on the stack
			dsm_funcs = new gm_DenseMatFunc_double();
			gp_funcs = new gm_GenPurposeFunc_double();
			load_marr_funcs_double(gm_handle, marr_funcs);
			load_dsm_funcs_double(gm_handle, dsm_funcs);
			load_gp_funcs_double(gm_handle, gp_funcs);
			FaustGPU<double>::marr_funcs = marr_funcs;
			FaustGPU<double>::dsm_funcs = dsm_funcs;
			FaustGPU<double>::gp_funcs = gp_funcs;
		}

	}


	template <>
		MatDense<double,Cpu> FaustGPU<double>::get_product(const bool transpose /* = false */, const bool conjugate /* = false */)
		{
			gm_Op op = OP_NOTRANSP;
			if(transpose)
				if(conjugate)
					op = OP_CONJTRANSP;
				else
					op = OP_TRANSP;
			gm_DenseMatFunc_double* dsm_funcs = (gm_DenseMatFunc_double*) this->dsm_funcs;
			gm_MatArrayFunc_double* marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
			double one;
			set_one<double>(&one);
			auto gpu_prod_mat_dense = marr_funcs->chain_matmul(gpu_mat_arr, one, op);
			int32_t nrows, ncols;
			dsm_funcs->info(gpu_prod_mat_dense, &nrows, &ncols);
			MatDense<double, Cpu> gpu2cpu_mat(nrows, ncols);
			dsm_funcs->tocpu(gpu_prod_mat_dense, (double*) reinterpret_cast<double*>(gpu2cpu_mat.getData()));
			dsm_funcs->free(gpu_prod_mat_dense);
			return gpu2cpu_mat;
		}

	template <>
		Vect<double, Cpu> FaustGPU<double>::multiply(const Vect<double,Cpu>& v, const bool transpose, const bool conjugate)
		{
			std::cout << "FaustGPU::multiply(Vect)" << std::endl;
			double one;
			set_one<double>(&one);
			int32_t out_size = this->ncols; // default is transpose here

			gm_Op op;
			if(transpose && conjugate)
				op = OP_CONJTRANSP;
			else if(transpose)
				op = OP_TRANSP;
			else
			{
				op = OP_NOTRANSP;
				out_size = this->nrows;
			}

			Vect<double, Cpu> out_vec(out_size);

			gm_MatArrayFunc_double* marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
			marr_funcs->chain_matmul_by_cpu_dsm_tocpu(gpu_mat_arr, one, op, (double*) reinterpret_cast<double*>(const_cast<double*>(v.getData())), v.size(), 1, (double*) reinterpret_cast<double*>(out_vec.getData()));

			return out_vec;
		}

	template <>
		MatDense<double, Cpu> FaustGPU<double>::multiply(const MatGeneric<double,Cpu>* A, const bool transpose, const bool conjugate)
		{
			std::cout << "FaustGPU::multiply(MatGeneric)" << std::endl;
			const MatSparse<double, Cpu>* sp_mat;
			const MatDense<double, Cpu>* ds_mat;
			int32_t out_nrows;
			double one;
			set_one<double>(&one);
			gm_Op op;
			if(transpose && conjugate)
				op = OP_CONJTRANSP;
			else if(transpose)
				op = OP_TRANSP;
			else
				op = OP_NOTRANSP;

			if(transpose)
				out_nrows = this->ncols;
			else
				out_nrows = this->nrows;

			MatDense<double, Cpu> out_mat(out_nrows, A->getNbCol());

			gm_MatArrayFunc_double* marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
			if(sp_mat = dynamic_cast<const MatSparse<double,Cpu>*>(A))
			{

				throw std::runtime_error("FaustGPU::multiply() by MatSparse isn't yet impl.");
				//		std::cout << "FaustGPU::multiply(MatSparse): " << sp_mat->getNbRow() << " " << sp_mat->getNbCol()<< " " << sp_mat->getNonZeros()<< std::endl;
				marr_funcs->togpu_spm(gpu_mat_arr, sp_mat->getNbRow(), sp_mat->getNbCol(), sp_mat->getNonZeros(), (int32_t *) sp_mat->getOuterIndexPtr(), (int32_t*) sp_mat->getInnerIndexPtr(), (double*) reinterpret_cast<double*>(const_cast<double*>(sp_mat->getValuePtr())));
			}
			else if(ds_mat = dynamic_cast<const MatDense<double,Cpu>*>(A))
			{

//				std::cout << "FaustGPU::multiply(MatDense): " << ds_mat->getNbRow() << " " << ds_mat->getNbCol()<< " " << ds_mat->getNonZeros()<< std::endl;
				marr_funcs->chain_matmul_by_cpu_dsm_tocpu(gpu_mat_arr, one, op, (double*) reinterpret_cast<double*>(const_cast<double*>(ds_mat->getData())), ds_mat->getNbRow(), ds_mat->getNbCol(), (double*) reinterpret_cast<double*>(out_mat.getData()));
			}

			return out_mat;
		}

	template <>
	void FaustGPU<double>::pop_front()
	{
		if(cpu_mat_ptrs.size() > 0)
		{
			if(use_ref_man)
				ref_man.release(*(cpu_mat_ptrs.begin()));
			auto marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
			marr_funcs->erase_at(gpu_mat_arr, 0, !use_ref_man);
			cpu_mat_ptrs.erase(cpu_mat_ptrs.begin());
		}
	}

	template <>
	void FaustGPU<double>::pop_back()
	{
		if(cpu_mat_ptrs.size() > 0)
		{
			if(use_ref_man)
				ref_man.release(*(cpu_mat_ptrs.end()-1));
			auto marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
			marr_funcs->erase_at(gpu_mat_arr, cpu_mat_ptrs.size()-1, !use_ref_man);
			cpu_mat_ptrs.erase(cpu_mat_ptrs.end()-1);
		}
	}

	template <>
	void FaustGPU<double>::push_back(MatGeneric<double,Cpu>* M)
	{
//		std::cout << "push_back M=" << M << " size: " << cpu_mat_ptrs.size() << " this:" << this << std::endl;
		MatSparse<double, Cpu>* sp_mat;
		MatDense<double, Cpu>* ds_mat;
		void* gpu_ref; //sp or ds mat

		auto dsm_funcs = (gm_DenseMatFunc_double*) this->dsm_funcs;
		auto marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
		auto gp_funcs = (gm_GenPurposeFunc_double*) this->gp_funcs;


		if(use_ref_man && cpu_gpu_map.find(M) != cpu_gpu_map.end())
		{
			// already known cpu, gpu mats

			// add the gpu matrix to gpu mat list
			if(dynamic_cast<MatDense<double,Cpu>*>(M))
			{
//				std::cout << "add the gpu dense matrix to gpu mat list" << std::endl;
				marr_funcs->addgpu_dsm(gpu_mat_arr, cpu_gpu_map[M]);
			}
			else
			{
//				std::cout << "add the gpu sparse matrix to gpu mat list" << std::endl;
				// M is sparse
				marr_funcs->addgpu_spm(gpu_mat_arr, cpu_gpu_map[M]);
			}

			cpu_mat_ptrs.push_back(M);

			ref_man.acquire(M);

			return;
		}

		if(sp_mat = dynamic_cast<MatSparse<double,Cpu>*>(M))
		{

//			std::cout << "FaustGPU(): " << sp_mat->getNbRow() << " " << sp_mat->getNbCol()<< " " << sp_mat->getNonZeros()<< std::endl;
			gpu_ref = marr_funcs->togpu_spm(gpu_mat_arr, sp_mat->getNbRow(), sp_mat->getNbCol(), sp_mat->getNonZeros(), sp_mat->getOuterIndexPtr(), sp_mat->getInnerIndexPtr(), (double*) reinterpret_cast<double*>(sp_mat->getValuePtr()));
//			std::cout << "after togpu_spm" << std::endl;
		}
		else if(ds_mat = dynamic_cast<MatDense<double,Cpu>*>(M))
		{
//			std::cout << "FaustGPU(): " << ds_mat->getNbRow() << " " << ds_mat->getNbCol()<< " " << ds_mat->getNonZeros()<< std::endl;
			gpu_ref = marr_funcs->togpu_dsm(gpu_mat_arr, ds_mat->getNbRow(), ds_mat->getNbCol(), (double*) reinterpret_cast<double*>(ds_mat->getData()));
//			std::cout << "after togpu_dsm" << std::endl;
		}

		cpu_mat_ptrs.push_back(M);
		if(use_ref_man)
		{
			cpu_gpu_map[M] = gpu_ref;
			ref_man.acquire(M);
		}
	}

	template <>
		FaustGPU<double>::~FaustGPU()
		{
			gm_users--;
			gm_MatArrayFunc_double* marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
			if(use_ref_man)
			{
				// release all gpu mats
				for(auto m: cpu_mat_ptrs)
				{
//					if(ref_man.contains(m)) // useless assuming we release only the acquired factors
						ref_man.release(m);
				}
			}
			marr_funcs->free(gpu_mat_arr, ! use_ref_man); // delete used mats only if it doesn't use ref_man
			if(gm_users <= 0)
			{
				gm_close_lib(gm_handle);
				delete (gm_MatArrayFunc_double*) marr_funcs;
//				gm_GenPurposeFunc_double* marr_funcs = (gm_GenPurposeFunc_double*) this->gp_funcs;
				delete (gm_DenseMatFunc_double*) dsm_funcs;
//				delete (gm_GenPurposeFunc_double*) gp_funcs;
			}
		}

	template<>
	Faust::RefManager FaustGPU<double>::ref_man([](void *fact)
		{
			gm_GenPurposeFunc_double* gp_funcs = (gm_GenPurposeFunc_double*) Faust::FaustGPU<double>::gp_funcs;
			//normally cpu_gpu_map must contains a the key fac if ref_man knew it (see ctor)
			gp_funcs->free_mat(Faust::FaustGPU<double>::cpu_gpu_map[fact]);
			Faust::FaustGPU<double>::cpu_gpu_map.erase(fact);
		});

	template<>
		FaustGPU<double>::FaustGPU(const std::vector<MatGeneric<double,Cpu>*>& factors) : use_ref_man(true)
	{
		//	std::cout << "FaustGPU<double>::FaustGPU()" << " marr_funcs:" << marr_funcs << std::endl;
		check_gpu_mod_loaded();
		this->load_gm_functions(); //lazy instantiation

		gm_DenseMatFunc_double* dsm_funcs = (gm_DenseMatFunc_double*) this->dsm_funcs;
		gm_MatArrayFunc_double* marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
		gm_GenPurposeFunc_double* gp_funcs = (gm_GenPurposeFunc_double*) this->gp_funcs;

		//	std::cout << "FaustGPU<double>::FaustGPU() factors size: " << factors.size() << " marr_funcs:" << marr_funcs << std::endl;
		gpu_mat_arr = marr_funcs->create();
		//	std::cout << "FaustGPU<double>::FaustGPU() factors size: " << factors.size() << std::endl;
		nrows = factors[0]->getNbRow();
		ncols = (*(factors.end()-1))->getNbCol();
		for(auto m: factors)
		{
			push_back(m);
		}
		//			std::cout << "FaustGPU() factors size, matarray size: " << factors.size() << " " << marr_funcs->size(gpu_mat_arr) << std::endl;
	}


	template <>
	void FaustGPU<double>::update(const Faust::MatGeneric<double,Cpu>* M, int32_t id)
	{
//		std::cout << "update M=" << M << " id=" << id << " this:" << this << std::endl;
		MatGeneric<double,Cpu>* M_ = const_cast<MatGeneric<double,Cpu>*>(M);
		// I promise I won't touch M_ data!
		if(M != cpu_mat_ptrs[id])
		{
//			std::cout << "M: " << M << ", cpu_mat_ptrs[id]: " << cpu_mat_ptrs[id] << std::endl;
			throw std::runtime_error("It's not authorized to update from another cpu matrix than the original one.");
		}

		gm_MatArrayFunc_double* marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
		MatSparse<double, Cpu>* sp_mat;
		MatDense<double, Cpu>* ds_mat;
		void* gpu_ref;

		// if the dims are not equal between M_ and the gpu mat, an exception will be raised by gpu_mod
		if(sp_mat = dynamic_cast<MatSparse<double,Cpu>*>(M_))
		{

//			std::cout << "FaustGPU::update() sparse: " << sp_mat->getNbRow() << " " << sp_mat->getNbCol()<< " " << sp_mat->getNonZeros()<< std::endl;
			gpu_ref = marr_funcs->cpu_set_spm_at(gpu_mat_arr, sp_mat->getNbRow(), sp_mat->getNbCol(), sp_mat->getNonZeros(), sp_mat->getOuterIndexPtr(), sp_mat->getInnerIndexPtr(), (double*) reinterpret_cast<double*>(sp_mat->getValuePtr()), id);
		}
		else if(ds_mat = dynamic_cast<MatDense<double,Cpu>*>(M_))
		{

//			std::cout << "FaustGPU::update() dense: " << ds_mat->getNbRow() << " " << ds_mat->getNbCol()<< " " << ds_mat->getNonZeros()<< std::endl;
			gpu_ref =  marr_funcs->cpu_set_dsm_at(gpu_mat_arr, ds_mat->getNbRow(), ds_mat->getNbCol(), (double*) reinterpret_cast<double*>(ds_mat->getData()), id);
		}
		// gpu_ref is not recorded because this is an assignment of data but the pointers don't change

	}

	template <>
	void FaustGPU<double>::insert(const Faust::MatGeneric<double,Cpu>* M, int32_t id)
	{
		MatGeneric<double,Cpu>* M_ = const_cast<MatGeneric<double,Cpu>*>(M);
		gm_MatArrayFunc_double* marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
		MatSparse<double, Cpu>* sp_mat;
		MatDense<double, Cpu>* ds_mat;
		void* gpu_ref;
		if(use_ref_man && cpu_gpu_map.find(M_) != cpu_gpu_map.end())
		{
			// already known cpu, gpu mats
			// add the gpu matrix to gpu mat list
//			std::cout << "insert at " << id << " the gpu matrix to gpu mat list reused M=" << M << " this:" << this << std::endl;
			if(dynamic_cast<MatDense<double,Cpu>*>(M_))
				marr_funcs->insert_dsm(gpu_mat_arr, cpu_gpu_map[M_], id);
			else
				// M is sparse
				marr_funcs->insert_spm(gpu_mat_arr, cpu_gpu_map[M_], id);

			cpu_mat_ptrs.insert(cpu_mat_ptrs.begin()+id,M_);

			ref_man.acquire(M_);

			return;
		}

		// if the dims are not equal between M_ and the gpu mat, an exception will be raised by gpu_mod
		if(sp_mat = dynamic_cast<MatSparse<double,Cpu>*>(M_))
		{

			gpu_ref = marr_funcs->togpu_insert_spm(gpu_mat_arr, sp_mat->getNbRow(), sp_mat->getNbCol(), sp_mat->getNonZeros(), sp_mat->getOuterIndexPtr(), sp_mat->getInnerIndexPtr(), (double*) reinterpret_cast<double*>(sp_mat->getValuePtr()), id);
		}
		else if(ds_mat = dynamic_cast<MatDense<double,Cpu>*>(M_))
		{

			gpu_ref =  marr_funcs->togpu_insert_dsm(gpu_mat_arr, ds_mat->getNbRow(), ds_mat->getNbCol(), (double*) reinterpret_cast<double*>(ds_mat->getData()), id);
		}

		cpu_mat_ptrs.insert(cpu_mat_ptrs.begin()+id, M_);


//		std::cout << "insert at " << id << " the gpu matrix to gpu mat list M=" << M << " this:" << this << std::endl;

		if(use_ref_man)
		{
			cpu_gpu_map[M_] = gpu_ref;
			ref_man.acquire(M_);
		}

	}

	template<>
		Real<double> FaustGPU<double>::spectral_norm(int32_t max_iter, Real<double> threshold)
		{
			std::cout << "FaustGPU::spectral_norm" << std::endl;
			gm_MatArrayFunc_double* marr_funcs = (gm_MatArrayFunc_double*) this->marr_funcs;
			return marr_funcs->spectral_norm(gpu_mat_arr, (float) threshold, max_iter);
		}
}
