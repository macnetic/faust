#include "faust_Slice.h"

using namespace Faust;


Slice::Slice(faust_unsigned_int start_id, faust_unsigned_int end_id)
{
	//enforce ascendant order: this->start_id <= this->end_id
	if(start_id <= end_id) {
		this->start_id = start_id;
		this->end_id = end_id;
	}
	else {
		this->end_id = start_id;
		this->start_id = end_id;
	}
}

Slice::Slice(Slice& s)
{
	this->copy(s);
}

Slice::Slice() :  start_id(0), end_id(0)
{
}

bool Slice::belong_to(faust_unsigned_int min_id, faust_unsigned_int max_id)
{
	Slice s(min_id, max_id);
	return belong_to(s);
}

bool Slice::belong_to(Slice& s)
{
	return this->start_id >= s.start_id && this->end_id <= s.end_id;
}

void Slice::copy(Slice& s)
{
	this->start_id = s.start_id;
	this->end_id = s.end_id;
}

void Slice::swap(Slice& s1, Slice& s2)
{
	Slice tmp(s1.start_id, s1.end_id);
	s1.start_id = s2.start_id;
	s1.end_id = s2.end_id;
	s2.start_id = tmp.start_id;
	s2.end_id = tmp.end_id;
}

