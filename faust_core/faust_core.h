#ifndef __FAUST_CORE_H__
#define __FAUST_CORE_H__

#include <vector>

class faust_mat;
class faust_spmat;
class faust_vec;

class faust_core
{
	public:
		faust_core();
		faust_core(const std::vector<faust_spmat>& facts);
		
		~faust_core();


	public:
		//void operator*=(const faust_core&  f);
		//faust_mat& operator*=(const faust_mat&   f)const;
		//faust_mat& operator*=(const faust_spmat& f)const;



	private:
		const std::vector<faust_spmat> data;
		// true if all facts have been multiplied to obtain the product matrix factProduct;
		bool factsMultiplied;
		faust_mat factProduct;

	friend faust_vec operator*(const faust_core& f, const faust_vec& v);
		
};



#endif
