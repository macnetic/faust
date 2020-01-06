#include "faust_MatDense.h"
#include "faust_MatSparse.h"

namespace Faust {
	template<typename FPP>
		void wht_factors(unsigned int n, std::vector<MatGeneric<FPP,Cpu>*>&  factors, const bool cloning_fact, const bool norma)
		{
			if(n == 0)
			{
				factors.resize(1);
				MatDense<FPP,Cpu>* mat = new MatDense<FPP,Cpu>(1,1);
				(*mat)[0] = 1;
				factors[0] = mat;
			}
			else
			{
				factors.resize(n);
				unsigned int order = 1ull << n;
				std::vector<int> col_ids(order), row_ids(order);
				std::vector<FPP> vals(order);
				unsigned int order_over_2 = order >> 1;
				unsigned int order_times_2 = order << 1;
				unsigned int i_times_2;

				//			// init the permutation matrix

				for(unsigned int i=0; i < order_over_2; i++)
				{
					i_times_2 = i << 1;
					row_ids[i] = i_times_2;
					row_ids[i+order_over_2] = i_times_2+1;
					col_ids[i] = i;
					col_ids[i+order_over_2] = i + order_over_2;
					vals[i] = vals[i+order_over_2] = 1;
				}
				//			cout << row_ids.size() << endl;
				//			cout << col_ids.size() << endl;
				//			cout << vals.size() << endl;
				MatSparse<FPP,Cpu> P(row_ids, col_ids, vals, order, order);
				P.update_dim();

				// init the base matrix
				int *row_ptr = new int[order+1];
				row_ptr[0] = 0;
				int *bcol_ids = new int[order_times_2];
				FPP *bvalues = new FPP[order_times_2];

				//			cout << "row_ptr: ";
				for(int i = 1; i < order+1;i++)
				{
					row_ptr[i] = 2+row_ptr[i-1];
					//				cout << row_ptr[i] << ",";
				}
				//			cout << endl;

				bool parity = true; //row parity

				int col_id = 0;
				//			cout << "bvalues: ";
				//			cout << "bcol_ids: ";
				for(unsigned int i=0; i < order_times_2; i+=2)
				{

					if(parity) //row index is pair
						bvalues[i] = bvalues[i+1] = 1;
					else
					{
						bvalues[i+1] = -1;
						bvalues[i] = 1;
					}
					//				cout << bvalues[i] << " " << bvalues[i+1];
					parity = ! parity;
					bcol_ids[i] = col_id;
					bcol_ids[i+1] = col_id+1;
					//				cout << bcol_ids[i] << " " << bcol_ids[i+1] << " ";
					if(((i + 1) & 3u) == 3u) col_id+=2; // i+1 mod 4 == 0
				}
				//			cout << endl;
				MatSparse<FPP,Cpu> B(order_times_2, order, order, bvalues, row_ptr, bcol_ids, false);
				//			cout << "Faust::TransformHelper::hadamardFaust(), B:" << endl;
				//			B.Display();
				delete[] bvalues;
				delete[] bcol_ids;
				delete[] row_ptr;
				MatSparse<FPP,Cpu>* factor = new MatSparse<FPP,Cpu>(order,order);
				factor->mat = B.mat*P.mat;
				factor->update_dim();

				factors[0] = factor;
				for(int i=1; i < n; i++)
					if(cloning_fact)
						factors[i] = factor->Clone();
					else
						factors[i] = factor;

				if(norma)
				{
					factors[0] = factor->Clone();
					*factors[0] *= static_cast<FPP>(1.0/sqrt((float)order));
				}

			}
		}
}
