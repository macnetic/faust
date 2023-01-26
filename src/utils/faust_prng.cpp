#include <cstdlib>
#include <iostream>
#include <random>
#include "faust_prng.h"

using namespace std;
namespace Faust
{

	unsigned int seed(unsigned seed/*=0*/)
	{
		static unsigned int seed_ = 0;
		if(seed == 0)
			return seed_;
		else
		{
			seed_ = seed;
			srand(seed_);
			default_random_engine gen(seed);
			generator(&gen);
		}
		return seed_;
	}

	default_random_engine& generator(default_random_engine *gen)
	{
		static default_random_engine generator;
		if(gen)
			generator = *gen;
		return generator;
	}

}
