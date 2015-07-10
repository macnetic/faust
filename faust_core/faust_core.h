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
		faust_core(const std::vector<faust_spmat>& facts, const faust_real lambda_);
		faust_core(const faust_params& params);
		
		void get_facts(std::vector<faust_spmat>& sparse_facts)const; 
		faust_real get_lambda()const;

		~faust_core(){}


	public:
		//void operator*=(const faust_core&  f);
		//faust_mat& operator*=(const faust_mat&   f)const;
		//faust_mat& operator*=(const faust_spmat& f)const;



	private:
		std::vector<faust_spmat> data;
		faust_real lambda;
		// true if all facts have been multiplied to obtain the product matrix factProduct;
		bool isDataInit;
		//faust_mat factProduct;

	friend faust_vec operator*(const faust_core& f, const faust_vec& v);
		
};



#endif
