#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include "faust_bit_rev_permu.h"

namespace Faust {
	template<typename FPP>
void fft_factors(unsigned int n, std::vector<MatGeneric<complex<FPP>,Cpu>*>&  v)
{
	//TODO: clean this code
	//Cooley-Tukey
	// number of factors to set : n+1
	v.resize(n+1);
	//TODO: rename index and new_index
	unsigned int dim_size = 1u << n, L, r, L_over_2, L_times_2;
	unsigned int* index = new unsigned int[dim_size];
	unsigned int* new_index = new unsigned int[dim_size];
	complex<FPP>* fac_data = new complex<FPP>[dim_size*2];
	int* fac_row_ptr = new int[dim_size+1]; //int and not uint because of MatSparse
	// we should change it later when MatSparse will be updated
	int* fac_col_ind = new int[dim_size*2];
	complex<FPP> wL;
	const FPP pi = std::acos(-1);
	const std::complex<double> cplx_i(0, 1); //TODO: rename to i (np for loop vars because they have their own scope)
	for(unsigned int i = 0; i < dim_size; i++)
		index[i] = i;
	memcpy(new_index, index, sizeof(unsigned int)*dim_size);
	bit_rev_permu(n, new_index, false);
	vector<complex<FPP>> ones(dim_size);
	for(typename vector<complex<FPP>>::iterator it=ones.begin(); it != ones.end(); it++)
		*it = complex<FPP>(1.0);
	MatSparse<complex<FPP>,Cpu> *P = new MatSparse<complex<FPP>,Cpu>(index, new_index, ones, dim_size, dim_size);
	MatSparse<complex<FPP>,Cpu> * factor;
	//			cout << "fft_factors() P:" << endl;
	//			P->Display();
	v[n] = P;
	for(unsigned int q=0; q < n; q++) //n+1 factors if counting P
	{
		L_over_2 = 1u << q;
		L = L_over_2 << 1;
		L_times_2 = L << 1;
		r = dim_size >> (q+1); // dim_size/L
		wL = std::exp(FPP(-2.0)*cplx_i*pi/FPP(L));
		//TODO: infer j from i and delete j
		//TODO: there is a potential for OpenMP opt. here and in the next loops
		for(unsigned int i=0, j=0; i < L_over_2; i++,j+=2)
		{
			fac_data[j] = 1; //identity diag
			fac_data[j+1] = std::pow(wL, index[i]);
			fac_data[j+L] = 1;
			fac_data[j+L+1] = -fac_data[j+1];
			//					cout << "fac_data[j+1]: " << fac_data[j+1] << endl;
		}
		//Rhs == fac_data[0], ..., fac_data[2*L-1]
		for(unsigned int i=1;i<r;i++)
			//TODO: var L_times_2
			memcpy(fac_data+i*(L<<1), fac_data, sizeof(complex<FPP>)*(L<<1));
		//ok for the factor data, now the indices (CSR format)
		//				cout << "fact_data[i]:" << endl;
		//				for(int i = 0; i < dim_size*2; i++)
		//					cout << fac_data[i] << " ";
		//				cout << endl;
		//row_ptr
		// Rhs (see CooleyTukeyFact.m)
		for(int i=0,j=0;i<L+1;i+=2,j++)
			fac_row_ptr[j] = i;
		//the first L_over_2+1 eles of fac_row_ptr are set
		// set now the second part of Rhs
		// (ignoring the first 0 when copying)
		memcpy(fac_row_ptr+L_over_2+1, fac_row_ptr+1, sizeof(int)*L_over_2);
		//shift by L the second part
		for(int i=L_over_2+1;i<L+1;i++)
			fac_row_ptr[i] += L;
		//				for(int i=0;i<L+1;i++)
		//					cout << "fac_row_ptr[i]=" << fac_row_ptr[i] << " ";
		//				cout << endl;
		// ok we have the first L+1 elements of fac_row_ptr
		// we need r*L+1 elements for Aq (whole factor)
		for(int i=1;i<r;i++)
		{
			//ignoring again the first ele. == 0
			memcpy(fac_row_ptr+i*L+1, fac_row_ptr+1, sizeof(int)*L);
			for(int j=i*L+1;j<(i+1)*L+1;j++)
				fac_row_ptr[j] += fac_row_ptr[i*L];
		}
		//				for(int i=0;i<dim_size+1;i++)
		//					cout << "fac_row_ptr[i]=" << fac_row_ptr[i] << " ";
		//				cout << endl;				//fac_row_ptr is fully ready to set the factor
		// now set the column indices
		// starting from Rhs indices
		for(int i=0;i<L_over_2;i++)
		{
			fac_col_ind[(i<<1)] = i;
			fac_col_ind[(i<<1)+1] = L_over_2+i;
		}
		//copy the first half height of the matrix columns into the second buf part
		// because they are the same
		memcpy(fac_col_ind+L, fac_col_ind, sizeof(int)*L);
		// we need r times the Rhs columns
		for(int i=1;i<r;i++)
		{
			memcpy(fac_col_ind+i*L_times_2,fac_col_ind,sizeof(int)*L_times_2);
			// shift because the Rhs are diagonal blocks in factor (Aq)
			for(int j=i*L_times_2;j<(i+1)*L_times_2;j++)
				fac_col_ind[j] += i*L;
		}
		//				for(int i=0;i<dim_size*2;i++)
		//					cout << "fac_col_ind[i]=" << fac_col_ind[i] << " ";
		//				cout << endl;
		// init the MatSparse with its buffers
		factor = new MatSparse<complex<FPP>,Cpu>(dim_size*2, dim_size, dim_size, fac_data,
				fac_row_ptr, fac_col_ind);
		//				cout << "fft_factors() factor:" <<endl;
		//				MatDense<complex<FPP>,Cpu>(*factor).Display();
		//				cout << "id: " << n-q-1 << endl;
		//				cout << factor << endl;
		//				factor->Display();
		v[n-q-1] = factor; //TODO: directly affects v with the MatSparse, delete factor

	}
	delete[] index;
	delete[] new_index;
	delete[] fac_data;
}

	template<typename FPP>
void fft_factors(unsigned int n, vector<MatGeneric<FPP,Cpu>*>&  v)
{
	handleError("Faust::TransformHelper<FPP,Cpu>", "Can't get fft factors for real matrices.");
}

}
