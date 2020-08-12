#include "faust_prox.h"
namespace Faust
{
	template<>
		void prox_skperm<double>(Faust::MatDense<double, Cpu> & M_, const unsigned int k, const bool normalized /* default to true*/, const bool pos/* default to false */)
		{
			//	template<typename FPP> using DenseMat = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
			//	using DenseMatInt = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
			// use preproc constants instead of nice C++ template aliases because they are only used here and can't be local to a function
#define DenseMatFPP Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
#define DenseMatRealFPP Eigen::Matrix<Real<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
#define DenseMatInt Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
			// M is a DenseMat<double>&, val is a double
#define set_all_to_self_minus(M, val) \
			for(int i = 0; i < M.rows(); i++) \
			for(int j = 0; j < M.cols(); j++) \
			M(i,j) -= val
			// use a Map to avoid data buffer copying
			Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> M(M_.getData(), M_.getNbRow(), M_.getNbCol());
			assert(M.rows() == M.cols());
			assert(k > 0);
			unsigned long shape[2] = {(unsigned long)M.rows(), (unsigned long)M.cols()};
			unsigned long size = shape[0];
			int side, vertex, v;
			DenseMatInt matching(M.rows(), M.cols()), visited(2,size);
			DenseMatInt degree(2, size), parent(2, size);
			DenseMatFPP result(M.rows(), M.cols());
			DenseMatFPP potential(2, size), slack(2, size);
			DenseMatRealFPP C(M.rows(), M.cols());

			C = (- square(M.array())).real(); // M.square() doesn't work... why ?

			auto inf = numeric_limits<double>::infinity();
			auto eps = numeric_limits<double>::epsilon();
			eps = 1e-7;
			inf = 1e10;

			matching.setZero();
			potential.setZero();
			potential.block(0,0, 1, size).setConstant(C.minCoeff()); //C.rowwise().minCoeff();
			degree.setZero();
			for(int i=0;i<k*size;i++)
			{
				queue<pair<int,int>> q;
				parent.setConstant(-1);
				visited.setZero();
				slack.setConstant(inf);
				//		slack.setConstant(1e-10);
				int unmatched_vertex = -1;

				for(int v=0; v < size; v++)
				{
					if(degree(0,v) < k)
					{
						q.push(std::make_pair(0,v));
						visited(0,v) = 1;
						assert(visited(0,v) == 1);
					}
				}

				while(true)
				{
					if(q.empty())
					{
						double delta = inf;
						for(int side=0;side<2;side++)
						{
							for(int vertex=0; vertex < size; vertex++)
								if(! visited(side, vertex))
									delta = std::min(delta, slack(side, vertex));
						}
						set_all_to_self_minus(slack, double(delta));
						for(int vertex=0; vertex < size; vertex++)
						{
							if(visited(0,vertex))
								potential(0, vertex) += delta;
							if(visited(1, vertex))
								potential(1, vertex) -= delta;
						}

						for(int side = 0; side < 2; side++)
						{
							for(int vertex = 0; vertex < size; vertex++)
							{
								if(abs(slack(side,vertex)) < eps && ! visited(side, vertex))
								{
									visited(side, vertex) = 1;
									q.push(make_pair(side, vertex));
								}
							}
						}
					}

					side = q.front().first;
					vertex = q.front().second;
					q.pop();

					if(side == 1 && degree(side, vertex) < k)
					{
						unmatched_vertex = vertex;
						degree(1, unmatched_vertex) += 1;
						break;
					}

					DenseMatRealFPP weight;
					DenseMatInt connected;

					if(side == 0)
					{
						weight.resize(1, C.cols());
						weight = C.block(vertex, 0, 1, C.cols());
						connected.resize(1, matching.cols());
						connected = matching.block(vertex, 0, 1, matching.cols());
					}
					else
					{
						weight.resize(C.rows(), 1);
						weight = C.block(0, vertex, C.rows(), 1);
						connected.resize(matching.rows(), 1);
						connected = matching.block(0, vertex, matching.rows(), 1);
						set_all_to_self_minus(connected, 1);
						connected *= -1;
					}

					for(int u=0; u < size; u++)
					{
						double p_diff = double(weight(u)) - potential(side, vertex) - potential(1 - side, u);
						if(std::abs(p_diff) > eps)
						{
							if(side == 0 && ! visited(1 - side, u) && (p_diff) > 0 && (p_diff) < (slack(1-side, u)))
							{
								slack(1, u) = p_diff;
								parent(1, u) = vertex;
							}
							if(side == 1 && ! visited(1 - side, u) && (p_diff) < 0 && - (p_diff) < (slack(1-side, u)))
							{
								slack(0, u) = - p_diff;
								parent(0, u) = vertex;
							}
							continue;
						}

						if(visited(1 - side, u) || connected(u) == 1)
							continue;

						q.push(make_pair(1 - side, u));
						parent(1 - side, u) = vertex;
						visited(1 - side, u) = 1;

					}


				}

				v = unmatched_vertex;

				while(true)
				{
					auto u = parent(1, v);
					auto p = parent(0, u);
					matching(u, v) = 1;
					if(p == -1)
					{
						degree(0, u) += 1;
						break;
					}
					else
					{
						matching(u, p) = 0;
						v = p;
					}
				}
			}
			M = (matching.array() > 0).select(M, 0);
		}

	//TODO: duplicate of specialization (with double) is to remove when a general impl. (for real and complex) will be available
	template<>
		void prox_skperm<float>(Faust::MatDense<float, Cpu> & M_, const unsigned int k, const bool normalized /* default to true*/, const bool pos/* default to false */)
		{
			//	template<typename FPP> using DenseMat = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
			//	using DenseMatInt = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
			// use preproc constants instead of nice C++ template aliases because they are only used here and can't be local to a function
#define DenseMatFPP Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
#define DenseMatRealFPP Eigen::Matrix<Real<float>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
#define DenseMatInt Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
			// M is a DenseMat<float>&, val is a float
#define set_all_to_self_minus(M, val) \
			for(int i = 0; i < M.rows(); i++) \
			for(int j = 0; j < M.cols(); j++) \
			M(i,j) -= val
			// use a Map to avoid data buffer copying
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> M(M_.getData(), M_.getNbRow(), M_.getNbCol());
			assert(M.rows() == M.cols());
			assert(k > 0);
			unsigned long shape[2] = {(unsigned long)M.rows(), (unsigned long)M.cols()};
			unsigned long size = shape[0];
			int side, vertex, v;
			DenseMatInt matching(M.rows(), M.cols()), visited(2,size);
			DenseMatInt degree(2, size), parent(2, size);
			DenseMatFPP result(M.rows(), M.cols());
			DenseMatFPP potential(2, size), slack(2, size);
			DenseMatRealFPP C(M.rows(), M.cols());

			C = (- square(M.array())).real(); // M.square() doesn't work... why ?

			auto inf = numeric_limits<float>::infinity();
			auto eps = numeric_limits<float>::epsilon();
			eps = 1e-7;
			inf = 1e10;

			matching.setZero();
			potential.setZero();
			potential.block(0,0, 1, size).setConstant(C.minCoeff()); //C.rowwise().minCoeff();
			degree.setZero();
			for(int i=0;i<k*size;i++)
			{
				queue<pair<int,int>> q;
				parent.setConstant(-1);
				visited.setZero();
				slack.setConstant(inf);
				//		slack.setConstant(1e-10);
				int unmatched_vertex = -1;

				for(int v=0; v < size; v++)
				{
					if(degree(0,v) < k)
					{
						q.push(std::make_pair(0,v));
						visited(0,v) = 1;
						assert(visited(0,v) == 1);
					}
				}

				while(true)
				{
					if(q.empty())
					{
						float delta = inf;
						for(int side=0;side<2;side++)
						{
							for(int vertex=0; vertex < size; vertex++)
								if(! visited(side, vertex))
									delta = std::min(delta, slack(side, vertex));
						}
						set_all_to_self_minus(slack, float(delta));
						for(int vertex=0; vertex < size; vertex++)
						{
							if(visited(0,vertex))
								potential(0, vertex) += delta;
							if(visited(1, vertex))
								potential(1, vertex) -= delta;
						}

						for(int side = 0; side < 2; side++)
						{
							for(int vertex = 0; vertex < size; vertex++)
							{
								if(abs(slack(side,vertex)) < eps && ! visited(side, vertex))
								{
									visited(side, vertex) = 1;
									q.push(make_pair(side, vertex));
								}
							}
						}
					}

					side = q.front().first;
					vertex = q.front().second;
					q.pop();

					if(side == 1 && degree(side, vertex) < k)
					{
						unmatched_vertex = vertex;
						degree(1, unmatched_vertex) += 1;
						break;
					}

					DenseMatRealFPP weight;
					DenseMatInt connected;

					if(side == 0)
					{
						weight.resize(1, C.cols());
						weight = C.block(vertex, 0, 1, C.cols());
						connected.resize(1, matching.cols());
						connected = matching.block(vertex, 0, 1, matching.cols());
					}
					else
					{
						weight.resize(C.rows(), 1);
						weight = C.block(0, vertex, C.rows(), 1);
						connected.resize(matching.rows(), 1);
						connected = matching.block(0, vertex, matching.rows(), 1);
						set_all_to_self_minus(connected, 1);
						connected *= -1;
					}

					for(int u=0; u < size; u++)
					{
						float p_diff = float(weight(u)) - potential(side, vertex) - potential(1 - side, u);
						if(std::abs(p_diff) > eps)
						{
							if(side == 0 && ! visited(1 - side, u) && (p_diff) > 0 && (p_diff) < (slack(1-side, u)))
							{
								slack(1, u) = p_diff;
								parent(1, u) = vertex;
							}
							if(side == 1 && ! visited(1 - side, u) && (p_diff) < 0 && - (p_diff) < (slack(1-side, u)))
							{
								slack(0, u) = - p_diff;
								parent(0, u) = vertex;
							}
							continue;
						}

						if(visited(1 - side, u) || connected(u) == 1)
							continue;

						q.push(make_pair(1 - side, u));
						parent(1 - side, u) = vertex;
						visited(1 - side, u) = 1;

					}


				}

				v = unmatched_vertex;

				while(true)
				{
					auto u = parent(1, v);
					auto p = parent(0, u);
					matching(u, v) = 1;
					if(p == -1)
					{
						degree(0, u) += 1;
						break;
					}
					else
					{
						matching(u, p) = 0;
						v = p;
					}
				}
			}
			M = (matching.array() > 0).select(M, 0);
		}

	template<>
		void prox_skperm<complex<double>>(Faust::MatDense<complex<double>, Cpu> & M_, const unsigned int k, const bool normalized /* default to true*/, const bool pos/* default to false */)
		{
			throw std::runtime_error("Complex prox_skperm is not yet implemented.");
		}

	//TODO: the complex impl. below doesn't work: a complex matrix with a nul imaginary part doesn't project to the same image as the real proj impl return
//	template<> void prox_skperm<complex<double>>(Faust::MatDense<complex<double>, Cpu> & M_, const unsigned int k, const bool normalized /* default to true*/, const bool pos/* default to false */)
//		{
//			//	template<typename FPP> using DenseMat = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
//			//	using DenseMatInt = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
//			// use preproc constants instead of nice C++ template aliases because they are only used here and can't be local to a function
//#define DenseMatFPP Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
//#define DenseMatRealFPP Eigen::Matrix<Real<complex<double>>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
//#define DenseMatInt Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
//			// M is a DenseMat<FPP>&, val is a FPP
//#define set_all_to_self_minus(M, val) \
//			for(int i = 0; i < M.rows(); i++) \
//			for(int j = 0; j < M.cols(); j++) \
//			M(i,j) -= val
//			// use a Map to avoid data buffer copying
//			Eigen::Map<Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>> M(M_.getData(), M_.getNbRow(), M_.getNbCol());
//			assert(M.rows() == M.cols());
//			assert(k > 0);
//			unsigned long shape[2] = {M.rows(), M.cols()};
//			unsigned long size = shape[0];
//			int side, vertex, v;
//			DenseMatInt matching(M.rows(), M.cols()), visited(2,size);
//			DenseMatInt degree(2, size), parent(2, size);
//			DenseMatFPP result(M.rows(), M.cols());
//			DenseMatFPP potential(2, size), slack(2, size);
//			DenseMatRealFPP C(M.rows(), M.cols());
//
//			C = (- square(M.array())).real(); // M.square() doesn't work... why ?
//
//			auto inf = numeric_limits<double>::infinity();
//			auto eps = numeric_limits<double>::epsilon();
//			eps = 1e-7;
//			inf = 1e10;
//
//			matching.setZero();
//			potential.setZero();
//			potential.block(0,0, 1, size).setConstant(C.minCoeff()); //C.rowwise().minCoeff();
//			degree.setZero();
//			for(int i=0;i<k*size;i++)
//			{
//				queue<pair<int,int>> q;
//				parent.setConstant(-1);
//				visited.setZero();
//				slack.setConstant(inf);
//				//		slack.setConstant(1e-10);
//				int unmatched_vertex = -1;
//
//				for(int v=0; v < size; v++)
//				{
//					if(degree(0,v) < k)
//					{
//						q.push(std::make_pair(0,v));
//						visited(0,v) = 1;
//						assert(visited(0,v) == 1);
//					}
//				}
//
//				while(true)
//				{
//					if(q.empty())
//					{
//						complex<double> delta = inf;
//						for(int side=0;side<2;side++)
//						{
//							for(int vertex=0; vertex < size; vertex++)
//								if(! visited(side, vertex))
//									delta = Faust::fabs(delta) < Faust::fabs(slack(side, vertex))?delta:slack(side,vertex);
//						}
//						set_all_to_self_minus(slack, delta);
//						for(int vertex=0; vertex < size; vertex++)
//						{
//							if(visited(0,vertex))
//								potential(0, vertex) += delta;
//							if(visited(1, vertex))
//								potential(1, vertex) -= delta;
//						}
//
//						for(int side = 0; side < 2; side++)
//						{
//							for(int vertex = 0; vertex < size; vertex++)
//							{
//								if(abs(slack(side,vertex)) < eps && ! visited(side, vertex))
//								{
//									visited(side, vertex) = 1;
//									q.push(make_pair(side, vertex));
//								}
//							}
//						}
//					}
//
//					side = q.front().first;
//					vertex = q.front().second;
//					q.pop();
//
//					if(side == 1 && degree(side, vertex) < k)
//					{
//						unmatched_vertex = vertex;
//						degree(1, unmatched_vertex) += 1;
//						break;
//					}
//
//					DenseMatRealFPP weight;
//					DenseMatInt connected;
//
//					if(side == 0)
//					{
//						weight.resize(1, C.cols());
//						weight = C.block(vertex, 0, 1, C.cols());
//						connected.resize(1, matching.cols());
//						connected = matching.block(vertex, 0, 1, matching.cols());
//					}
//					else
//					{
//						weight.resize(C.rows(), 1);
//						weight = C.block(0, vertex, C.rows(), 1);
//						connected.resize(matching.rows(), 1);
//						connected = matching.block(0, vertex, matching.rows(), 1);
//						set_all_to_self_minus(connected, 1);
//						connected *= -1;
//					}
//
//					for(int u=0; u < size; u++)
//					{
//						auto p_diff = complex<double>(weight(u)) - potential(side, vertex) - potential(1 - side, u);
//						if(Faust::fabs(p_diff) > eps)
//						{
//							if(side == 0 && ! visited(1 - side, u) && Faust::fabs(p_diff) > 0 && Faust::fabs(p_diff) < Faust::fabs(slack(1-side, u)))
//							{
//								slack(1, u) = p_diff;
//								parent(1, u) = vertex;
//							}
//							if(side == 1 && ! visited(1 - side, u) && Faust::fabs(p_diff) < 0 && - Faust::fabs(p_diff) < Faust::fabs(slack(1-side, u)))
//							{
//								slack(0, u) = - p_diff;
//								parent(0, u) = vertex;
//							}
//							continue;
//						}
//
//						if(visited(1 - side, u) || connected(u) == 1)
//							continue;
//
//						q.push(make_pair(1 - side, u));
//						parent(1 - side, u) = vertex;
//						visited(1 - side, u) = 1;
//
//					}
//
//
//				}
//
//				v = unmatched_vertex;
//
//				while(true)
//				{
//					auto u = parent(1, v);
//					auto p = parent(0, u);
//					matching(u, v) = 1;
//					if(p == -1)
//					{
//						degree(0, u) += 1;
//						break;
//					}
//					else
//					{
//						matching(u, p) = 0;
//						v = p;
//					}
//				}
//			}
//			M = (matching.array() > 0).select(M, 0);
//		}

}
