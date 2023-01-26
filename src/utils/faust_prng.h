namespace Faust
{

	//TODO: a proper class with seed and generator as attributes and static member functions

	/**
	 * This module is for C PRNG initialization in Faust.
	 */

	/**
	 * \brief Initializes the C PRNG and C++ default_random_engine with a particular seed or returns the current one.
	 *
	 * \param seed: if 0 the current seed is returned (default behaviour) otherwise it is the new seed for the PRNGs. 
	 */
	unsigned int seed(unsigned int seed=0);


	/**
	 * \brief Gets the C++ default_random_engine used in Faust or assigns it if gen is not nullptr.
	 */
	std::default_random_engine& generator(std::default_random_engine *gen = nullptr);
}

