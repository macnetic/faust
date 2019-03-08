namespace Faust {
	void bit_rev_permu(unsigned int n, unsigned int* v, const bool initv)
	{
		unsigned int size = 1<<n;
		if(initv)
			for(unsigned int i=0;i<size;i++)
				v[i] = i;
		unsigned int lower_mask = 1;
		unsigned int upper_mask = 1<<(n-1);
		unsigned int shift=0;
		unsigned short lobit, hibit;

		while(lower_mask < upper_mask)
		{
			for(unsigned int i=0;i<size;i++){
				lobit = (v[i] & lower_mask) >>(shift);
				hibit = (v[i] & upper_mask)>>(n-shift-1);
				if(lobit > hibit){
					//lobit == 1, hibit == 0
					v[i] ^= lower_mask;
					v[i] |= upper_mask;
				}
				else if (lobit < hibit){
					//swap lobit == 0, hibit == 1
					v[i] |= lower_mask;
					v[i] ^= upper_mask;
				}
			}
			lower_mask <<= 1;
			upper_mask >>= 1;
			shift++;
		}
	}
}
