/**
 *
 *  Init v to (1,..,2^n-1).
 *  And proceeds with a bit reversal permutation.
 *  E.g. : if n = 3, v=(0, 1, 2, 3, 4, 5, 6, 7, 8) and becomes (0, 4, 2, 6, 1, 5, 3, 7).
 */
void bit_rev_permu(unsigned int n, unsigned int* v, const bool initv=true);

