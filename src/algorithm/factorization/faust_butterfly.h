#ifndef __FAUST_BUTTERFLY__
#define __FAUST_BUTTERFLY__
namespace Faust
{
	enum ButterflyFactDir
	{
		LEFT,
		RIGHT,
		BALANCED
	};

	enum FactMeth
	{
		DTO,
		LIFTING
	};

	template<typename FPP>
		void retrieveCEC(const Faust::MatDense<FPP, Cpu>& s1, const Faust::MatDense<FPP, Cpu>& s2, vector<vector<faust_unsigned_int>*> &cec, vector<vector<faust_unsigned_int>*> &noncec);

	template<typename FPP>
		void solveDTO(const Faust::MatDense<FPP, Cpu>& A, const Faust::MatDense<FPP, Cpu>& s1, const Faust::MatDense<FPP, Cpu>& s2, Faust::MatDense<FPP, Cpu>& X, Faust::MatDense<FPP, Cpu>& Y);

	template<typename FPP>
		void lifting_two_layers_factorization(const Faust::MatDense<FPP, Cpu>& A, const Faust::MatSparse<FPP, Cpu>& s1, const Faust::MatSparse<FPP, Cpu>& s2, Faust::MatDense<FPP, Cpu>& X, Faust::MatDense<FPP, Cpu>& Y);

	template<typename FPP>
		void lifting_two_layers_factorization(const Faust::MatDense<FPP, Cpu>& A, const Faust::MatDense<FPP, Cpu>& s1, const Faust::MatDense<FPP, Cpu>& s2, Faust::MatDense<FPP, Cpu>& X, Faust::MatDense<FPP, Cpu>& Y);

	template<typename FPP>
		void lifting_two_layers_factorization(const Faust::MatGeneric<FPP, Cpu>& A, const Faust::MatSparse<FPP, Cpu>& s1, const Faust::MatSparse<FPP, Cpu>& s2, Faust::MatSparse<FPP, Cpu>& X, Faust::MatSparse<FPP, Cpu>& Y);

	template<typename FPP>
		void butterfly_support(int nfactors, std::vector<Faust::MatSparse<FPP, Cpu>*> out);

	template<typename FPP>
		std::vector<Faust::MatSparse<FPP, Cpu>*> support_DFT(int nfactors);

	template<typename FPP>
		void bit_reversal_factor(int nfactors, std::vector<Faust::MatSparse<FPP, Cpu>*>& out);

	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical(const Faust::MatDense<FPP, Cpu>& A, const std::vector<Faust::MatSparse<FPP, Cpu>*> &supports, const ButterflyFactDir& dir=RIGHT);

	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical_balanced(const Faust::MatDense<FPP, Cpu>& A, const std::vector<Faust::MatSparse<FPP, Cpu>*> &supports, FactMeth meth = LIFTING);

	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical(const Faust::MatDense<FPP, Cpu>& A, const ButterflyFactDir &dir=RIGHT);


};
#include "faust_butterfly.hpp"
#endif
