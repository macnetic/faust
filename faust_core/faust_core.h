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
		faust_core(const std::vector<faust_spmat>& facts, const faust_real lambda_ = 1.0);
		faust_core(const faust_params& params);
		
		void get_facts(std::vector<faust_spmat>& sparse_facts)const; 
		int size()const{return data.size();} 
                faust_mat get_product();

		long long int get_total_nnz()const;

		~faust_core(){}


	public:
		// add all of the sparse matrices from f.data to this->data
		void operator*=(const faust_core&  f);
		// add the sparse matrix S to this->data
		void operator*=(const faust_spmat&  S);




	private:
		std::vector<faust_spmat> data;
		bool isDataInit;
		long long int totalNonZeros;

	friend faust_vec operator*(const faust_core& f, const faust_vec& v);
	friend faust_mat operator*(const faust_core& f, const faust_mat& M);
		
};



#endif
