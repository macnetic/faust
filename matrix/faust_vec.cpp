#include "faust_vec.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

void faust_vec::setOnes()
{
	vec.setOnes();
}


void faust_vec::Display() const
{
	 std::cout << "dim = " << getDim() << std::endl;
	 std::cout << vec << std::endl;
}


void faust_vec::resize(const int new_dim)
{
	
		if (new_dim <0)
		{
			std::cerr << "ERROR FAUST_VEC resize : les nouvelles dimensions doivent etre strictement positive" << std::endl;
			exit( EXIT_FAILURE);
		}
		else if (dim != new_dim)
		{
			dim = new_dim;
			vec.conservativeResize(dim);
			
		}
}

void faust_vec::operator=(faust_vec const& y)
{
	  vec = y.vec;
	  dim = y.dim;

}

