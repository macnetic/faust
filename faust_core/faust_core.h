#ifndef __FAUST_CORE_H__
#define __FAUST_CORE_H__

#include <vector>
#include "faust_spmat.h"

class faust_mat;
//class faust_spmat;
class faust_vec;
class faust_params;

class faust_core
{
	public:
		faust_core();
		faust_core(const std::vector<faust_spmat>& facts, const faust_real lambda_ = (faust_real)1.0);
		faust_core(const faust_params& params);
		faust_core(const faust_core & A);
		void get_facts(std::vector<faust_spmat>& sparse_facts)const{sparse_facts = data;}; 
		int size()const{return data.size();} 
                faust_mat get_product();
		faust_spmat get_fact(int id) const;		
		int getNbRow() const;
		int getNbCol() const;
		long long int get_total_nnz()const{return totalNonZeros;}
		void clear(){data.resize(0);totalNonZeros=0;}
		void push_back(const faust_spmat& S);
		void push_first(const faust_spmat& S);
		void pop_back(faust_spmat& S);
		void pop_first(faust_spmat& S);
		void pop_first(faust_spmat& S) const;
		void Display()const;
		void transpose();
		//(*this) = (*this) * A
		void multiply(const faust_core & A);
		//(*this) = A * (*this)
		void multiplyLeft(const faust_core & A);	
		void scalarMultiply(const faust_real scalar){data[0]*=scalar;}
		faust_real spectralNorm(const int nbr_iter_max, faust_real threshold, int &flag) const;
		~faust_core(){}


	public:
		void operator=(const faust_core&  f){data=f.data;totalNonZeros=f.totalNonZeros;}
		// add all of the sparse matrices from f.data to this->data
		void operator*=(const faust_real  scalar){scalarMultiply(scalar);};
		void operator*=(const faust_core&  f){multiply(f);};
		// add the sparse matrix S to this->data
		void operator*=(const faust_spmat&  S){push_back(S);totalNonZeros+=S.getNonZeros();}




	private:
		std::vector<faust_spmat> data;
		long long int totalNonZeros;

		

	friend faust_vec operator*(const faust_core& f, const faust_vec& v);
	friend faust_mat operator*(const faust_core& f, const faust_mat& M);
	friend void multiply(const faust_core & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, char typeA, char typeMult);
		
};



#endif
