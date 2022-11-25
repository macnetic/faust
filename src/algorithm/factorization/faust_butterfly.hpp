#include <Eigen/SVD>
#include "faust_MatDense.h"
#include "faust_TransformHelper.h"
#include <vector>
#include "faust_openmp.h"
#ifdef OMP_ENABLED
#include "omp.h"
#endif
#include <cmath>
#include <limits>
#ifdef USE_GPU_MOD
#include "faust_MatDense_gpu.h"
#endif

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
			//TODO: this function should be deleted because lifting_two_layers_factorization is faster
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
		void lifting_two_layers_factorization(const Faust::MatDense<FPP, Cpu>& A, const Faust::MatSparse<FPP, Cpu>& s1, const Faust::MatSparse<FPP, Cpu>& s2, Faust::MatDense<FPP, Cpu>& X, Faust::MatDense<FPP, Cpu>& Y)
		{
			// TODO: this function should be deleted because lifting_two_layers_factorization(MatGeneric, ...) def is faster
			int r = s1.getNbCol();
			Faust::MatDense<FPP, Cpu> subA, u, v;
			X.resize(s1.getNbRow(), s1.getNbCol());
			Y.resize(s2.getNbRow(), s2.getNbCol());
			X.setZeros();
			Y.setZeros();
			std::vector<MatDense<FPP, Cpu>> U(r), V(r);
			std::vector<std::vector<int>> rows(r), cols(r);
			#pragma omp parallel for schedule(dynamic) private(subA)
			for(int t=0;t<r;t++)
			{
				rows[t] = s1.col_nonzero_inds(t);
				cols[t] = s2.row_nonzero_inds(t);
				A.submatrix(rows[t], cols[t], subA);
				//#define BUTTERFLY_APPROX_RANK1
#ifdef BUTTERFLY_APPROX_RANK1
				subA.approx_rank1(U[t], V[t]);
#else
				if(A.getNbRow() >= (1 << 14) && (std::is_same<FPP, double>::value || std::is_same<FPP, std::complex<double>>::value))
				{
					// use jacobi svd as a workaround to this bug https://gitlab.com/libeigen/eigen/-/issues/1723
					// TODO: remove when it'll be fixed
					Eigen::JacobiSVD<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> svd;
					subA.initJacobiSVD(svd);
					subA.best_low_rank(1, U[t], V[t], svd);
				}
				else
					subA.best_low_rank(1, U[t], V[t]);
#endif
			}
			#pragma omp parallel for schedule(dynamic)
			for(int t=0;t<r;t++)
			{
				X.set_col_coeffs(t, rows[t], U[t], 0);
				Y.set_row_coeffs(t, cols[t], V[t], 0);
			}
		}

	template<typename FPP>
		void lifting_two_layers_factorization(const Faust::MatDense<FPP, Cpu>& A, const Faust::MatDense<FPP, Cpu>& s1, const Faust::MatDense<FPP, Cpu>& s2, Faust::MatDense<FPP, Cpu>& X, Faust::MatDense<FPP, Cpu>& Y)
		{
			// TODO: this function should be deleted because lifting_two_layers_factorization(MatGeneric, ...) def is faster
			int r = s1.getNbCol();
			Faust::MatDense<FPP, Cpu> subA, u, v;
			X.resize(s1.getNbRow(), s1.getNbCol());
			Y.resize(s2.getNbRow(), s2.getNbCol());
			X.setZeros();
			Y.setZeros();
			std::vector<MatDense<FPP, Cpu>> U(r), V(r);
			std::vector<std::vector<int>> rows(r), cols(r);
			#pragma omp parallel for schedule(dynamic) private(subA)
			for(int t=0;t<r;t++)
			{
				rows[t] = s1.col_nonzero_inds(t);
				cols[t] = s2.row_nonzero_inds(t);
				A.submatrix(rows[t], cols[t], subA);
				//#define BUTTERFLY_APPROX_RANK1
#ifdef BUTTERFLY_APPROX_RANK1
				subA.approx_rank1(U[t], V[t]);
#else
				subA.best_low_rank(1, U[t], V[t]);
#endif
			}
			#pragma omp parallel for schedule(dynamic)
			for(int t=0;t<r;t++)
			{
				X.set_col_coeffs(t, rows[t], U[t], 0);
				Y.set_row_coeffs(t, cols[t], V[t], 0);
			}
		}

	template<typename FPP>
		void lifting_two_layers_factorization(const Faust::MatGeneric<FPP, Cpu>& A, const Faust::MatSparse<FPP, Cpu>& s1, const Faust::MatSparse<FPP, Cpu>& s2, Faust::MatSparse<FPP, Cpu>& X, Faust::MatSparse<FPP, Cpu>& Y)
		{
			const Faust::MatDense<FPP, Cpu>* dsA = dynamic_cast<const Faust::MatDense<FPP, Cpu>*>(&A);
			const Faust::MatSparse<FPP, Cpu>* spA;
			if(dsA == nullptr)
			{
				spA = dynamic_cast<const Faust::MatSparse<FPP, Cpu>*>(&A);
				if(spA == nullptr)
					throw std::runtime_error("lifting_two_layers_factorization A must be MatDense or MatSparse.");
			}
			int r = s1.getNbCol();
			MatDense<FPP, Cpu> subA;
			X.resize(s1.getNbRow(), s1.getNbCol());
			Y.resize(s2.getNbRow(), s2.getNbCol());
			X.setZeros();
			Y.setZeros();

			bool svd_on_gpu = false;
			std::vector<std::vector<int>> rows(r), cols(r);

			#pragma omp parallel for schedule(dynamic)
			for(int t=0;t<r;t++)
			{
				rows[t] = s1.col_nonzero_inds(t);
				cols[t] = s2.row_nonzero_inds(t);
			}


			if(svd_on_gpu && rows[0].size() <= 32 && cols[0].size() <= 32) // GPU batched svd is limited to 32 nrows/ncols max
			{
#ifdef USE_GPU_MOD
				compute_XY_on_gpu(A, s1, s2, X, Y, rows, cols);
#else
				throw std::runtime_error("GPU SVD is not enabled in this version of FAµST.");
#endif
			}
			else
			{
				std::vector<MatDense<FPP, Cpu>> U(r), V(r);
#pragma omp parallel for schedule(dynamic) private(subA)
				for(int t=0;t<r;t++)
				{
					if(dsA)
						dsA->submatrix(rows[t], cols[t], subA);
					else
						spA->submatrix(rows[t], cols[t], subA);
					//#define BUTTERFLY_APPROX_RANK1
#ifdef BUTTERFLY_APPROX_RANK1
					subA.approx_rank1(U[t], V[t]);
#else
					if(A.getNbRow() >= (1 << 14) && (std::is_same<FPP, double>::value || std::is_same<FPP, std::complex<double>>::value))
					{
						// use jacobi svd as a workaround to this bug https://gitlab.com/libeigen/eigen/-/issues/1723
						// TODO: remove when it'll be fixed
						Eigen::JacobiSVD<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> svd;
						subA.initJacobiSVD(svd);
						subA.best_low_rank(1, U[t], V[t], svd);
					}
					else
						subA.best_low_rank(1, U[t], V[t]);
#endif
				}

				std::vector<Eigen::Triplet<FPP>> XtripletList;
				std::vector<Eigen::Triplet<FPP>> YtripletList;
				for(int t=0;t<r;t++)
				{
					for(int i=0;i<rows[t].size();i++)
						XtripletList.push_back(Eigen::Triplet<FPP>(rows[t][i], t, U[t].getData()[i]));
					for(int i=0;i<cols[t].size();i++)
						YtripletList.push_back(Eigen::Triplet<FPP>(t, cols[t][i], V[t](0,i)));
				}
				X = MatSparse<FPP, Cpu>(XtripletList, s1.getNbRow(), s1.getNbCol());
				Y = MatSparse<FPP, Cpu>(YtripletList, s2.getNbRow(), s2.getNbCol());
			}
		}

#ifdef USE_GPU_MOD
	template<typename FPP>
		void compute_XY_on_gpu(const Faust::MatGeneric<FPP, Cpu>& A, const Faust::MatSparse<FPP, Cpu>& s1, const Faust::MatSparse<FPP, Cpu>& s2, Faust::MatSparse<FPP, Cpu>& X, Faust::MatSparse<FPP, Cpu>& Y, std::vector<std::vector<int>> &rows, std::vector<std::vector<int>> &cols)
		{
			const Faust::MatDense<FPP, Cpu>* dsA = dynamic_cast<const Faust::MatDense<FPP, Cpu>*>(&A);
			const Faust::MatSparse<FPP, Cpu>* spA;
			if(dsA == nullptr)
			{
				spA = dynamic_cast<const Faust::MatSparse<FPP, Cpu>*>(&A);
				if(spA == nullptr)
					throw std::runtime_error("lifting_two_layers_factorization A must be MatDense or MatSparse.");
			}
			int r = s1.getNbCol();
			int m = rows[0].size(), n = cols[0].size();
			MatDense<FPP, Cpu> Us(m, r), Vs(n, r);
			MatDense<Real<FPP>, Cpu> Ss(r);
			MatDense<FPP, Cpu> subAs(m, r * n);
			for(int i=0; i < r; i++)
				if(dsA)
					dsA->submatrix(rows[i], cols[i], subAs.getData() + i * m * n);
				else
					spA->submatrix(rows[i], cols[i], subAs.getData() + i * m * n);

			Faust::batched_svd(subAs, r /* batch_sz */, Us, Vs, Ss, /* rank */ 1);


			using Map = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;

			for(int i=0; i < r; i++)
			{
				Map U(Us.getData() + i * m, m, 1);
				U *= Ss(i);
				Vs.adjoint();
			}

			std::vector<Eigen::Triplet<FPP>> XtripletList;
			std::vector<Eigen::Triplet<FPP>> YtripletList;
			for(int t=0;t<r;t++)
			{
				Map U_(Us.getData() + t * m, m, 1);
				Map V_(Vs.getData() + t * n, 1, n);
				for(int i=0;i<rows[t].size();i++)
					XtripletList.push_back(Eigen::Triplet<FPP>(rows[t][i], t, U_(i, 0)));
				for(int i=0;i<cols[t].size();i++)
					YtripletList.push_back(Eigen::Triplet<FPP>(t, cols[t][i], V_(0, i)));
			}
			X = MatSparse<FPP, Cpu>(XtripletList, s1.getNbRow(), s1.getNbCol());
			Y = MatSparse<FPP, Cpu>(YtripletList, s2.getNbRow(), s2.getNbCol());
		}
#endif

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
//				Faust::MatDense<FPP, Cpu> submat;
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
			#pragma omp parallel for
			for(int i=0;i<noncec.size();i++)
			{
				Faust::MatDense<FPP, Cpu> submat, bestx, besty;
				auto ce = noncec[i];
//			for(auto ce: noncec)
//			{
				auto r = (*ce)[0];
				auto RP = s1.col_nonzero_inds(r);
				auto CP = s2.row_nonzero_inds(r);
				A.submatrix(RP, CP, submat);
#ifdef BUTTERFLY_APPROX_RANK1
				submat.approx_rank1(bestx, besty);
#else
				submat.best_low_rank(ce->size(), bestx, besty);
#endif
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
		std::vector<Faust::MatSparse<FPP, Cpu>*> butterfly_support(int nfactors, Faust::MatSparse<FPP, Cpu>* P/*=nullptr*/)
		{
			// WARNING: free the MatSparse after calling this function
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
			if(P != nullptr)
				out.push_back(P); // the caller is responsible to delete the matrix P
			return out;
		}

	// TODO: specializations should be in .cpp file
	template<>
		void bit_reversal_factor<double>(int nfactors, Faust::MatSparse<double, Cpu>*& out)
		{
			//nothing to do for double matrices
		}

	template<>
		void bit_reversal_factor<float>(int nfactors, Faust::MatSparse<float, Cpu>*& out)
		{
			//nothing to do for float matrices
		}

	template<>
		void bit_reversal_factor<std::complex<double>>(int nfactors, Faust::MatSparse<std::complex<double>, Cpu>*& out)
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
			out = new MatSparse<std::complex<double>,Cpu>(index, new_index, ones, dim_size, dim_size);
		}

	template<>
		void bit_reversal_factor<double>(int nfactors, std::vector<Faust::MatSparse<double, Cpu>*> &out)
		{
			//nothing to do for double matrices
		}

	template<>
		void bit_reversal_factor<float>(int nfactors, std::vector<Faust::MatSparse<float, Cpu>*> &out)
		{
			//nothing to do for float matrices
		}

	template<>
		void bit_reversal_factor<std::complex<double>>(int nfactors, std::vector<Faust::MatSparse<std::complex<double>, Cpu>*>& out)
		{
			MatSparse<std::complex<double>,Cpu> *P;
			bit_reversal_factor(nfactors, P);
			out.push_back(P);
		}

	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical(const Faust::MatDense<FPP, Cpu>& A, const std::vector<Faust::MatSparse<FPP, Cpu>*> &supports, const ButterflyFactDir& dir/*=RIGHT*/)
		{
			using ID_FUNC = std::function<int(int)>;
			auto th = new TransformHelper<FPP, Cpu>();
			int i, j;
			Faust::MatDense<FPP, Cpu> s2;
			std::vector<Faust::MatDense<FPP, Cpu>> s2_vec(supports.size());
			Faust::MatDense<FPP, Cpu> X, Y;
			Faust::MatDense<FPP, Cpu> mat = A;
			assert(dir == RIGHT || dir == LEFT);
			ID_FUNC next_s1_id_left = [&supports](int i){if(i > 1) return i-1; else return -1;};
			ID_FUNC next_s1_id_right = [&supports](int i){if(i < supports.size()-2) return i+1; else return -1;};
			ID_FUNC next_s1_id = dir == RIGHT?next_s1_id_right:next_s1_id_left;
			ID_FUNC next_s2_id_left = [&supports, &i](int j){if(j<i-1) return j+1; else return -1;};
			ID_FUNC next_s2_id_right = [&supports, &i](int j){if(j>i+1) return j-1; else return -1;};
			ID_FUNC next_s2_id = RIGHT==dir?next_s2_id_right:next_s2_id_left;

			s2.resize(A.getNbRow(), A.getNbCol());
			s2.setEyes();
			i = dir == RIGHT?0:supports.size()-1;
			j = dir == RIGHT?supports.size()-1:0;
			do
			{
				supports[j]->multiply(s2, 'N');
				s2.setNZtoOne();
				s2_vec[dir == RIGHT?j-1:j+1] = s2;
			}
			while((j = next_s2_id(j)) > -1);
			do
			{
				auto s1 = Faust::MatDense<FPP, Cpu>(*supports[i]);
				if(dir == RIGHT)
				{
					solveDTO(mat, s1, s2_vec[i], X, Y);
					th->push_back(&X);
				}
				else
				{
					solveDTO(mat, s2_vec[i], s1, Y, X);
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
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical_balanced(const Faust::MatDense<FPP, Cpu>& A, const std::vector<Faust::MatSparse<FPP, Cpu>*> &supports, const FactMeth meth/*=LIFTING*/)
		{
			//TODO: refactor with other prototypes
			if(meth != LIFTING)
				throw std::runtime_error("Only FactMeth::LIFTING is supported here, use the other prototype for DTO.");
			auto th = new TransformHelper<FPP, Cpu>();
			double l2_nfacts = std::log2((double)supports.size());
			size_t d = ceil(l2_nfacts); // depth of factorization tree
			/** initialize the support tree structure (the size of each vector -- one per level of tree) **/
			std::vector<std::vector<MatSparse<FPP, Cpu>>> tree_supports(d);
			// set supports for each tree level
			size_t i=0;
			while(i<d-1)
			{
				tree_supports[i] = std::vector<MatSparse<FPP, Cpu>>(1 << (i+1));
//				std::cout << "tree_supports["<<i<<"].size():" << tree_supports[i].size() << std::endl;
				i++;
			}
			// the number of leafs of tree is not necessarily a power of two (remaining factorization(s))
			size_t n = supports.size();
			if(d > 1)
				tree_supports[d-1] = std::vector<MatSparse<FPP, Cpu>>(2*(n - tree_supports[d-2].size()));
			else // factorization in two factors only
				tree_supports[d-1] = std::vector<MatSparse<FPP, Cpu>>(2);
//			std::cout << "tree_supports[" << d-1 << "].size():" << tree_supports[d-1].size() << std::endl;
			/** initialize the supports in the tree starting from the leafs -- that are the supports arguments **/
			for(int i=d-1;i>=0;i--)
			{
				#pragma omp parallel for schedule(dynamic)
				for(int j = 0;j < tree_supports[i].size(); j++)
				{
					if(i == d - 1)
					{
						// tree_supports[i][j] is a leaf on the last level
//						std::cout << "tree_supports[" << i << "," << j << "] received support[" << j << "]" << std::endl;
						tree_supports[i][j] = *supports[j];
					}
					else if (2*(j+1) > tree_supports[i+1].size())
					{
						// tree_supports[i][j] is a leaf
//						std::cout << "tree_supports[" << i << "," << j << "] received support[" << j+tree_supports[i+1].size()/2 << "]" << std::endl;
						tree_supports[i][j] = *supports[j+tree_supports[i+1].size()/2];
					}
					else
					{
						// tree_supports[i][j] is not a leaf so it's calculated using child nodes
						// tree_supports[i][j] = tree_supports[i+1][2*j] * tree_supports[i+1][2*j+]
//						std::cout << "tree_supports[" << i << "," << j << "] received tree_supports[" << i+1 << "," << 2*j << "]*tree_supports[" << i+1 << "," << 2*j+1 << "]" << std::endl;
						MatDense<FPP, Cpu> s;
						gemm_gen(tree_supports[i+1][2*j], tree_supports[i+1][2*j+1], s, FPP(1.0), FPP(0), 'N', 'N');
						tree_supports[i][j] = s;
					}

				}
			}
			// the supports are now set for all factorization tree nodes
			// factorize according to the tree of supports, level by level
			std::vector<MatSparse<FPP, Cpu>> next_lvl_facs;
			std::vector<MatSparse<FPP, Cpu>> input_facs;
			// first factorization
			next_lvl_facs.resize(2);
			lifting_two_layers_factorization(A, tree_supports[0][0], tree_supports[0][1], next_lvl_facs[0], next_lvl_facs[1]);
#ifdef BUTTERFLY_BALANCED_CHECK_ERRORS
			std::cout << "computing error... lvl1" << std::endl;
			MatDense<FPP, Cpu> L = next_lvl_facs[0];
			MatDense<FPP, Cpu> R = next_lvl_facs[1];
			L.multiplyRight(R);
			L -= A;
			std::cout << "err:" << L.norm()/A.norm() << std::endl;
#endif
			input_facs = next_lvl_facs;

			for(size_t i=0;i<d-1;i++)
			{
				next_lvl_facs.resize(tree_supports[i+1].size());
				int j_bound = std::min(tree_supports[i].size(), tree_supports[i+1].size()/2);
				#pragma omp parallel for schedule(static)
				for(int j=0;j < j_bound;j++)
					lifting_two_layers_factorization(input_facs[j], tree_supports[i+1][2*j], tree_supports[i+1][2*j+1], next_lvl_facs[j*2], next_lvl_facs[j*2+1]);
#ifdef BUTTERFLY_BALANCED_CHECK_ERRORS
				for(int j=0;j < j_bound;j++)
				{
					std::cout << "computing error... lvl" << (i+2) << " factorization of index" << j << " factor" << std::endl;
					MatDense<FPP, Cpu> A = input_facs[j];
					MatDense<FPP, Cpu> L = next_lvl_facs[j*2];
					MatDense<FPP, Cpu> R = next_lvl_facs[j*2+1];
					L.multiplyRight(R);
					L -= A;
					std::cout << "err:" << L.norm()/A.norm() << std::endl;
				}
#endif

				if(i < d - 2)
				{
					input_facs = next_lvl_facs;
				}
			}
			// build the resulting Faust
			for(auto f: next_lvl_facs)
			{
				th->push_back(new MatSparse<FPP, Cpu>(f), false, false);
			}
			if(next_lvl_facs.size() < supports.size())
				for(int i=next_lvl_facs.size()/2;i<input_facs.size();i++)
				{
					th->push_back(new MatSparse<FPP, Cpu>(input_facs[i]), false, false);
				}
			return th;
		}


	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical_balanced_dense(const Faust::MatDense<FPP, Cpu>& A, const std::vector<Faust::MatSparse<FPP, Cpu>*> &supports, const FactMeth meth/*=LIFTING*/)
		{
			//TODO: this function should be deleted because it has shown to be very inefficient relatively to its MatSparse counterpart
			auto th = new TransformHelper<FPP, Cpu>();
			double l2_nfacts = std::log2((double)supports.size());
			size_t d = ceil(l2_nfacts); // depth of factorization tree
			/** initialize the support tree structure (the size of each vector -- one per level of tree) **/
			std::vector<std::vector<MatDense<FPP, Cpu>>> tree_supports(d);
			// set supports for each tree level
			size_t i=0;
			while(i<d-1)
			{
				tree_supports[i] = std::vector<MatDense<FPP, Cpu>>(1 << (i+1));
//				std::cout << "tree_supports["<<i<<"].size():" << tree_supports[i].size() << std::endl;
				i++;
			}
			// the number of leafs of tree is not necessarily a power of two (remaining factorization(s))
			size_t n = supports.size();
			if(d > 1)
				tree_supports[d-1] = std::vector<MatDense<FPP, Cpu>>(2*(n - tree_supports[d-2].size()));
			else // factorization in two factors only
				tree_supports[d-1] = std::vector<MatDense<FPP, Cpu>>(2);
//			std::cout << "tree_supports[" << d-1 << "].size():" << tree_supports[d-1].size() << std::endl;
			/** initialize the supports in the tree starting from the leafs -- that are the supports arguments **/
			for(int i=d-1;i>=0;i--)
			{
				#pragma omp parallel for schedule(dynamic)
				for(int j = 0;j < tree_supports[i].size(); j++)
				{
					if(i == d - 1)
					{
						// tree_supports[i][j] is a leaf on the last level
//						std::cout << "tree_supports[" << i << "," << j << "] received support[" << j << "]" << std::endl;
						tree_supports[i][j] = *supports[j];
					}
					else if (2*(j+1) > tree_supports[i+1].size())
					{
						// tree_supports[i][j] is a leaf
//						std::cout << "tree_supports[" << i << "," << j << "] received support[" << j+tree_supports[i+1].size()/2 << "]" << std::endl;
						tree_supports[i][j] = *supports[j+tree_supports[i+1].size()/2];
					}
					else
					{
						// tree_supports[i][j] is not a leaf so it's calculated using child nodes
						// tree_supports[i][j] = tree_supports[i+1][2*j] * tree_supports[i+1][2*j+]
//						std::cout << "tree_supports[" << i << "," << j << "] received tree_supports[" << i+1 << "," << 2*j << "]*tree_supports[" << i+1 << "," << 2*j+1 << "]" << std::endl;
						gemm(tree_supports[i+1][2*j], tree_supports[i+1][2*j+1], tree_supports[i][j], FPP(1.0), FPP(0), 'N', 'N');
					}

				}
			}
			// the supports are now set for all factorization tree nodes
			// factorize according to the tree of supports, level by level
			std::vector<MatDense<FPP, Cpu>> next_lvl_facs;
			std::vector<MatDense<FPP, Cpu>> input_facs;
			// first factorization
			next_lvl_facs.resize(2);
			if(meth == DTO)
				solveDTO(A, tree_supports[0][0], tree_supports[0][1], next_lvl_facs[0], next_lvl_facs[1]);
			else // meth == ButterflyFactDir
				lifting_two_layers_factorization(A, tree_supports[0][0], tree_supports[0][1], next_lvl_facs[0], next_lvl_facs[1]);
			input_facs = next_lvl_facs;

			for(size_t i=0;i<d-1;i++)
			{
				next_lvl_facs.resize(tree_supports[i+1].size());
				int j_bound = std::min(tree_supports[i].size(), tree_supports[i+1].size()/2);
				#pragma omp parallel for schedule(static)
				for(int j=0;j < j_bound;j++)
				{
					if(meth == DTO)
						solveDTO(input_facs[j], tree_supports[i+1][2*j], tree_supports[i+1][2*j+1], next_lvl_facs[j*2], next_lvl_facs[j*2+1]);
					else // meth == LIFTING
						lifting_two_layers_factorization(input_facs[j], tree_supports[i+1][2*j], tree_supports[i+1][2*j+1], next_lvl_facs[j*2], next_lvl_facs[j*2+1]);
				}
				if(i < d - 2)
				{
					input_facs = next_lvl_facs;
				}
			}
			// build the resulting Faust
			for(auto f: next_lvl_facs)
			{
				th->push_back(new MatSparse<FPP, Cpu>(f), false, false);
			}
			if(next_lvl_facs.size() < supports.size())
				for(int i=next_lvl_facs.size()/2;i<input_facs.size();i++)
				{
					th->push_back(new MatSparse<FPP, Cpu>(input_facs[i]), false, false);
				}
			return th;
		}

	template<typename FPP>
		Faust::TransformHelper<FPP, Cpu>* butterfly_hierarchical(const Faust::MatDense<FPP, Cpu>& A, const ButterflyFactDir &dir/*=RIGHT*/, Faust::MatSparse<FPP, Cpu>* P/*=nullptr*/, const bool mul_perm /*= true*/)
		{
			double log2_size = log2((double)A.getNbRow());
			if(A.getNbRow() != A.getNbCol())
				throw std::runtime_error("The matrix to factorize must be square.");
			if(log2_size - int(log2_size) > std::numeric_limits<Real<FPP>>::epsilon())
				throw std::runtime_error("The matrix to factorize must be of a size equal to a power of two.");
			// let the permutation P argument to nullptr (default value) because we don't want to add it to the factorization supports
			// we rather multiply the matrix by P in the next of the function
			auto support = butterfly_support<FPP>((int) log2_size);

			TransformHelper<FPP, Cpu>* th = nullptr;
			MatSparse<FPP, Cpu> *Pt = nullptr;
			auto A_ = A;
			if (P != nullptr)
			{
				Pt = new MatSparse<FPP, Cpu>(*P);
				Pt->transpose();
				A_.multiplyRight(*Pt);
			}
			if(dir == BALANCED)
				th = butterfly_hierarchical_balanced(A_, support);
			else
				th = butterfly_hierarchical(A_, support, dir);
			for(int i = 0; i < support.size();i++)
				delete support[i];
			if(P != nullptr)
			{
				//TODO: maybe a swapcols on th would be wiser/quicker
				const MatSparse<FPP, Cpu>* sp_last_fac;
				const MatDense<FPP, Cpu>* ds_last_fac;
				auto new_last_fac = new MatSparse<FPP, Cpu>(*P);
				if(mul_perm)
				{
					if(sp_last_fac = dynamic_cast<const MatSparse<FPP, Cpu>*>(th->get_gen_fact(th->size()-1)))
					{
						sp_last_fac->multiply(*new_last_fac, 'N');
					}
					else if(ds_last_fac = dynamic_cast<const MatDense<FPP, Cpu>*>(th->get_gen_fact(th->size()-1)))
					{
						ds_last_fac->multiply(*new_last_fac, 'N');
					}
					else
						// it can't happen but still
						throw std::runtime_error("butterfly supports only MatSparse and MatDense matrix");
					th->replace(new_last_fac, th->size()-1);
				}
				else
				{
					// the permutation is not multiplied
					th->push_back(new_last_fac);
				}
				if(Pt != nullptr)
					delete Pt;
			}
			return th;
		}

};
