
#include "faust_constant.h"
#ifndef __FAUST_SLICE__
#define __FAUST_SLICE__
namespace Faust
{


	struct Slice
	{

		faust_unsigned_int start_id; // first index of slice
		faust_unsigned_int end_id; // end index of slice (excluded)

		Slice(faust_unsigned_int start_id, faust_unsigned_int end_id);
		Slice(const Slice&);
		Slice();

		bool belong_to(faust_unsigned_int min_i, faust_unsigned_int max_i);
		bool belong_to(Slice& s);
		void copy(const Slice& s);
		Slice& operator=(const Slice&);

		void display() const;
		static void swap(Slice& s1, Slice& s2);
		size_t size() const;

	};

}
#endif
