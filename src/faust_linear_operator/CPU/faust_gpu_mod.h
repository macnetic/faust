#ifndef __FAUST_GPU__
#define __FAUST_GPU__
#include "faust_RefManager.h"
#include "faust_gpu_mod_utils.h"
#include <exception>
#include <iostream>
#include <vector>
#include <map>
#include "faust_MatGeneric.h"
#include "faust_Transform.h"

namespace Faust
{

	template<typename FPP>
	class FaustGPU
	{

		static RefManager ref_man;
		public:
#ifdef _MSC_VER
		public: // should not be public but Visual Studio 14 (and only it) can't access private members from lambda exp (error C2248) // cf. ref_man ini.
				//TODO: set back to private later (maybe with a more recent version)
#endif
#ifdef _MSC_VER
		private:
#endif
		gm_MatArray_t gpu_mat_arr;
		std::vector<void*> cpu_mat_ptrs; // cpu mat addresses for which gpu copies are stored in gpu_mat_arr

		// this map is used to retrieve a cpu mat addr from a gpu mat addr
		// this is needed because the RefManager works on cpu mat (to use gpu mat, bijection relationship)
#ifdef _MSC_VER
		public: // should not be public but Visual Studio 14 (and only it) can't access private members from lambda exp (error C2248) // cf. ref_man ini.
		// TODO: set back to private later (mayber with a more recent version)
#endif
		static std::map<void*,void*> cpu_gpu_map;
#ifdef _MSC_VER
		private:
#endif
		int32_t nrows;
		int32_t ncols;

		bool use_ref_man; // default to true

		public:
		FaustGPU(const std::vector<Faust::MatGeneric<FPP,Cpu>*>&);
//		FaustGPU(const Transform<FPP,Cpu>*);
		~FaustGPU();


		MatDense<FPP, Cpu> get_product(const bool transpose = false, const bool conjugate = false);
		MatDense<FPP, Cpu> multiply(const Faust::MatGeneric<FPP,Cpu>*, const bool transpose = false, const bool conjugate = false);

		Vect<FPP, Cpu> multiply(const Faust::Vect<FPP,Cpu>&, const bool transpose = false, const bool conjugate = false);

		void pop_front();

		void pop_back();
		void push_back(MatGeneric<FPP,Cpu>* M);

		void insert(const Faust::MatGeneric<FPP,Cpu>* M, int32_t id);

		/** TODO: this function should be deleted because it definitely slower than the version 2 */
		FPP power_iteration(int32_t max_iter, Real<FPP> threshold, int& flag);

		/** This version differs from power_iteration in the way it computes F'Fv in this order F'(Fv)
		 * so that it doesn't need to load all factors of F'F into GPU memory but only F.
		 */
		FPP power_iteration2(int32_t max_iter, Real<FPP> threshold, int& flag);

		Real<FPP> spectral_norm(int32_t max_iter, Real<FPP> threshold);

		/* Update on gpu the copy matrix of M (M must have already been loaded otherwise an exception is raised) */
		void update(const Faust::MatGeneric<FPP,Cpu>* M, int32_t id);

		/* Returns true if the matrix has already been loaded on GPU, false otherwise */
		static bool is_cpu_mat_known(const Faust::MatGeneric<FPP,Cpu>*);
		/* Returns true if is_cpu_mat_known returns true for each matrix in mats, false otherwise */
		static bool are_cpu_mat_all_known(const std::vector<MatGeneric<FPP,Cpu>*> mats);

	};

	template<typename T>
		void set_one(T* scal);

}

#include "faust_gpu_mod.hpp"

#endif
