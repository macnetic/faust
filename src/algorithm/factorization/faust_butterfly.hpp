#include "faust_MatDense.h"
#include "faust_TransformHelper.h"
#include <vector>
#include "faust_openmp.h"
#ifdef OMP_ENABLED
#include "omp.h"
#endif
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
#ifdef OMP_ENABLED
			int nthreads = 8;
#else
			int nthreads = 1;
#endif
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
#ifdef OMP_ENABLED
								th_class[omp_get_thread_num()].push_back(j);
#else
								th_class[0].push_back(j);
#endif
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
			for(auto ce: cec)
			{
				auto r = (*ce)[0];
				auto RP = s1.col_nonzero_inds(r);
				auto CP = s2.row_nonzero_inds(r);
				if(RP.size() == 0 || CP.size() == 0)
					continue;
				//TODO: maybe some refactoring is possible between the for loops
				if(ce->size() == RP.size() || ce->size() == RP.size())
				{
					noncec.push_back(ce);
					continue;
				}
				A.submatrix(RP, CP, submat);
				if(ce->size() >= RP.size())
				{
					MatDense<FPP, Cpu> eye(RP.size(), ce->size());
					eye.setEyes();
					for(int i=0;i< ce->size();i++)
					{
						auto col_id = (*ce)[i];
						X.set_col_coeffs(col_id, RP, eye, i);
						if(i < RP.size())
						{
							auto row_id = (*ce)[i];
							Y.set_row_coeffs(row_id, CP, submat, i);
						}
					}
				}
				else
				{
					MatDense<FPP, Cpu> eye(ce->size(), CP.size());
					eye.setEyes();
					for(int i=0;i< ce->size();i++)
					{

						if(i < CP.size())
						{
							auto col_id = (*ce)[i];
							X.set_col_coeffs(col_id, RP, eye, i);
						}
						auto row_id = (*ce)[i];
						Y.set_row_coeffs(row_id, CP, submat, i);
					}

				}
			}
			for(auto ce: noncec)
			{
				auto r = (*ce)[0];
				auto RP = s1.col_nonzero_inds(r);
				auto CP = s2.row_nonzero_inds(r);
				A.submatrix(RP, CP, submat);
				submat.best_low_rank(ce->size(), bestx, besty);
				//TODO: OpenMP ?
				for(int i=0;i< ce->size();i++)
				{

					auto col_id = (*ce)[i];
					X.set_col_coeffs(col_id, RP, bestx, i);
					auto row_id = (*ce)[i];
					Y.set_row_coeffs(row_id, CP, besty, i);
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
			if(! std::is_same<FPP, Real<FPP>>::value)
			{
				bit_reversal_factor(nfactors, out);
			}
			return out;
		}

	template<>
		void bit_reversal_factor<double>(int nfactors, std::vector<Faust::MatSparse<double, Cpu>*> &out)
		{
			//nothing to do for double matrices
		}

	template<>
		void bit_reversal_factor<std::complex<double>>(int nfactors, std::vector<Faust::MatSparse<std::complex<double>, Cpu>*>& out)
		{
			// bit reversal permutation factor
			unsigned int dim_size = 1u << nfactors, L, r, L_over_2, L_times_2;
			unsigned int* index = new unsigned int[dim_size];
			unsigned int* new_index = new unsigned int[dim_size];
			for(unsigned int i = 0; i < dim_size; i++)
				index[i] = i;
			memcpy(new_index, index, sizeof(unsigned int)*dim_size);
			bit_rev_permu(nfactors, new_index, false);
			std::vector<std::complex<Real<double>>> ones(dim_size);
			for(typename std::vector<std::complex<double>>::iterator it=ones.begin(); it != ones.end(); it++)
				*it = std::complex<double>(1.0);
			MatSparse<std::complex<double>,Cpu> *P = new MatSparse<std::complex<double>,Cpu>(index, new_index, ones, dim_size, dim_size);
			out.push_back(P);
		}

	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical(const Faust::MatDense<FPP, Cpu>& A, const std::vector<Faust::MatSparse<FPP, Cpu>*> &supports, const ButterflyFactDir& dir/*=RIGHT*/)
		{
			using ID_FUNC = std::function<int(int)>;
			auto th = new TransformHelper<FPP, Cpu>();
			int i, j;
			Faust::MatDense<FPP, Cpu> s2;
			Faust::MatDense<FPP, Cpu> X, Y;
			Faust::MatDense<FPP, Cpu> mat = A;
			assert(dir == RIGHT || dir == LEFT);
			ID_FUNC next_s1_id_left = [&supports](int i){if(i > 1) return i-1; else return -1;};
			ID_FUNC next_s1_id_right = [&supports](int i){if(i < supports.size()-2) return i+1; else return -1;};
			ID_FUNC next_s1_id = dir == RIGHT?next_s1_id_right:next_s1_id_left;
			ID_FUNC next_s2_id_left = [&supports](int j){if(j>0) return j-1; else return -1;};
			ID_FUNC next_s2_id_right = [&supports, &i](int j){if(j>i+1) return j-1; else return -1;};
			ID_FUNC next_s2_id = RIGHT==dir?next_s2_id_right:next_s2_id_left;
			i = dir == RIGHT?0:supports.size()-1;
			do
			{
				auto s1 = Faust::MatDense<FPP, Cpu>(*supports[i]);
				s2.resize(s1.getNbRow(), s1.getNbCol());
				s2.setEyes();
				j = dir == RIGHT?supports.size()-1:i-1;
				do
				{
					supports[j]->multiply(s2, 'N');
				}
				while((j = next_s2_id(j)) > -1);
				s2.setNZtoOne();
				if(dir == RIGHT)
				{
					solveDTO(mat, s1, s2, X, Y);
					th->push_back(&X);
				}
				else
				{
					solveDTO(mat, s2, s1, Y, X);
					th->push_first(&X);
				}
				mat = Y;
				std::cout << "factorization #" << (dir==RIGHT?i+1:supports.size()-i) << std::endl;
			}
			while((i = next_s1_id(i)) > -1);
			if(dir == RIGHT)
				th->push_back(&mat);
			else
				th->push_first(&mat);
			return th;
		}

	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical(const Faust::MatDense<FPP, Cpu>& A, const ButterflyFactDir &dir/*=RIGHT*/)
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
			auto th = butterfly_hierarchical(A, support, dir);
			return th;
		}

};
