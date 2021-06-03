#include "faust_MatDense.h"
#include "faust_TransformHelper.h"
#include <vector>
#include "omp.h" //TODO: ifdef to avoid when BUILD_MULTITHREAD is OFF
#include <cmath>
#include <limits>

namespace Faust
{
	void print_classes(vector<vector<faust_unsigned_int>>& classes)
	{
		for(auto c: classes)
		{
			std::cout << "{";
			for(auto k: c)
			{
				std::cout << k << ", ";
			}
			std::cout << "}" << std::endl;
		}
	}

	void print_classes(vector<vector<faust_unsigned_int>*>& classes)
	{
		for(auto c: classes)
		{
			std::cout << "{";
			for(auto k: *c)
			{
				std::cout << k << ", ";
			}
			std::cout << "}" << " ";
		}
		std::cout << std::endl;
	}
	template<typename FPP>
		void retrieveCEC(const Faust::MatDense<FPP, Cpu>& s1, const Faust::MatDense<FPP, Cpu>& s2, vector<vector<faust_unsigned_int>*> &cec, vector<vector<faust_unsigned_int>*> &noncec)
		{
			faust_unsigned_int r = s1.getNbCol(); // == s2.getNbCol()
			//	Vect<FPP, Cpu> c1(r), c2(r), r1(r), r2(r);
			vector<bool> processed_ids(r, false);
			auto eps = 1e-6;
			int nthreads = 8;
			auto th_class = new vector<faust_unsigned_int>[nthreads];
			for(int i=0;i<r;i++)
			{
				if(! processed_ids[i])
				{
					// create a new equiv. class
					auto class_ = new vector<faust_unsigned_int>();
					class_->push_back(i);
#pragma omp parallel for num_threads(nthreads)
					for(int j=i+1;j<r;j++)
					{
						if(! processed_ids[j])
						{
							if(s1.eq_cols(s2, i, j, eps) && s1.eq_rows(s2, i, j, eps))
							{
								th_class[omp_get_thread_num()].push_back(j);
								processed_ids[j] = true;
							}
						}
					}
					for(int i=0;i<nthreads;i++)
					{
						for(auto k: th_class[i])
							class_->push_back(k);
						th_class[i].clear();
					}
					// class is in cec or noncec
					if(min(Faust::fabs(s1.sum_row(i)), Faust::fabs(s1.sum_col(i))) <= class_->size())
					{
						cec.push_back(class_);
					}
					else
					{
						noncec.push_back(class_);
					}
					processed_ids[i] = true;
				}
			}
			delete[] th_class;
		}

	template<typename FPP>
		void solveDTO(const Faust::MatDense<FPP, Cpu>& A, const Faust::MatDense<FPP, Cpu>& s1, const Faust::MatDense<FPP, Cpu>& s2, Faust::MatDense<FPP, Cpu>& X, Faust::MatDense<FPP, Cpu>& Y)
		{
			vector<vector<faust_unsigned_int>*> cec, noncec;
			Faust::MatDense<FPP, Cpu> submat, bestx, besty;
			X.resize(s1.getNbRow(), s1.getNbCol());
			Y.resize(s2.getNbRow(), s2.getNbCol());
			X.setZeros();
			Y.setZeros();
			retrieveCEC(s1, s2, cec, noncec);
			for(auto ce: noncec)
			{
				auto r = (*ce)[0];
				//		std::cout << "r:" << r << std::endl;
				auto RP = s1.col_nonzero_inds(r);
				//		for(int i=0;i<RP.size();i++)
				//			std::cout << RP[i] << " ";
				//		std::cout << std::endl;
				auto CP = s2.row_nonzero_inds(r);
				//		for(int i=0;i<CP.size();i++)
				//			std::cout << CP[i] << " ";
				//		std::cout << std::endl;
				A.submatrix(RP, CP, submat);
				submat.best_low_rank(ce->size(), bestx, besty);
				//		std::cout << "bestx.norm():" << bestx.norm() << std::endl;
				//		std::cout << bestx.getNbRow() << "," << bestx.getNbCol() << std::endl;
				//		std::cout << "besty.norm():" << besty.norm() << std::endl;
				//		std::cout << besty.getNbRow() << "," << besty.getNbCol() << std::endl;
				//		std::cout << "submat.norm():" << submat.norm() << std::endl;
				//		submat.Display();
				//TODO: OpenMP ?
				for(int i=0;i< ce->size();i++)
				{

					auto col_id = (*ce)[i];
					X.set_col_coeffs(col_id, RP, bestx, i);
					//			cout << "X" << endl;
					//			X.Display();
					//			std::cout << "colx=" << col_id << std::endl;
					//			std::cout << "rowx=";
					//			for(auto x: RP)
					//				std::cout << x << ", ";
					//			std::cout << std::endl;
					auto row_id = col_id;
					Y.set_row_coeffs(row_id, CP, besty, i);
					//			std::cout << "coly=";
					//			for(auto x: CP)
					//				std::cout << x << ", ";
					//			std::cout << std::endl;
					//			std::cout << "rowy=";
					//			std::cout << row_id << std::endl;

					//			std::cout << "X.norm():" << X.norm() << std::endl;
					//			std::cout << "Y.norm():" << Y.norm() << std::endl;
					//			auto X_nnz_inds = X.nonzeros_indices();
					//			std::cout << " X nonzeros:" << std::endl;
					//			for(auto c: X_nnz_inds)
					//				std::cout << c.first << " " << c.second << std::endl;
					//			auto Y_nnz_inds = Y.nonzeros_indices();
					//			std::cout << " Y nonzeros:" << std::endl;
					//			for(auto c: Y_nnz_inds)
					//				std::cout << c.first << " " << c.second << std::endl;
				}
			}
		}

	template<typename FPP>
		void butterfly_support(int size, Faust::MatDense<FPP, Cpu> & out)
		{
			assert(size % 2 == 0);
			Faust::MatDense<FPP, Cpu> ones(2, 2);
			auto s = size / 2;
			Faust::MatDense<FPP, Cpu> id(s, s);
			ones.setOnes();
			id.setEyes();
			Faust::kron(ones, id, out);
		}

	template<typename FPP>
		std::vector<Faust::MatSparse<FPP, Cpu>*> support_DFT(int nfactors)
		{
			std::vector<Faust::MatSparse<FPP, Cpu>*> out;
			auto size = 1 << nfactors;
			Faust::MatDense<FPP, Cpu> kernel, id;
			Faust::MatDense<FPP, Cpu> support;
			Faust::MatSparse<FPP, Cpu> * sp_support;
			for(int i=0;i < nfactors; i++)
			{
				butterfly_support(1 << (nfactors-i), kernel);
				id.resize(1 << i, 1 << i);
				id.setEyes();
				Faust::kron(id, kernel, support);
				sp_support = new Faust::MatSparse<FPP, Cpu>(support);
				out.push_back(sp_support);
			}
			return out;
		}

	template<typename FPP>
Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical(const Faust::MatDense<FPP, Cpu>& A, const std::vector<Faust::MatSparse<FPP, Cpu>*> &supports)
{
	Faust::TransformHelper<FPP, Cpu>* th = new TransformHelper<FPP, Cpu>();
	Faust::MatDense<FPP, Cpu> s2;
	Faust::MatDense<FPP, Cpu> X, Y;
	Faust::MatDense<FPP, Cpu> mat = A;
	for(int i=0;i<supports.size()-1;i++)
	{
		auto s1 = Faust::MatDense<FPP, Cpu>(*supports[i]);
		s2.resize(s1.getNbRow(), s1.getNbCol());
		s2.setEyes();
		for(int j=supports.size()-1;j>=i+1;j--)
		{
			supports[j]->multiply(s2, 'N');
		}
		s2.setNZtoOne();
		solveDTO(mat, s1, s2, X, Y);
		mat = Y;
		th->push_back(&X);
//		std::cout << X.norm() << std::endl;
		std::cout << "factorization #" << i << std::endl;
	}
	th->push_back(&mat);
	return th;
}

	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical(const Faust::MatDense<FPP, Cpu>& A)
	{
		double log2_size = log2((double)A.getNbRow());
		if(A.getNbRow() != A.getNbCol())
			throw std::runtime_error("The matrix to factorize must be square.");
		if(log2_size - int(log2_size) > std::numeric_limits<Real<FPP>>::epsilon())
			throw std::runtime_error("The matrix to factorize must be of a size equal to a power of two.");
		auto support = support_DFT<FPP>((int) log2_size);
		//	std::cout << "support norms" << std::endl;
		//	for(auto s: support)
		//		std::cout << s->norm() << std::endl;
		auto th = butterfly_hierarchical(A, support);
		return th;
	}

};
