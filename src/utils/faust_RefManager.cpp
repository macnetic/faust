#include "faust_RefManager.h"

void Faust::RefManager::acquire(void* ref)
{
	if(refCounts.find(ref) != refCounts.end())
	{
		refCounts[ref]++;
	}
	else
	{
		refCounts.insert(make_pair(ref,1));
	}
}

void Faust::RefManager::release(void* ref)
{
	if(refCounts.find(ref) != refCounts.end())
	{
		if(refCounts[ref] <= 1)
		{
			this->cb(ref);
			refCounts.erase(ref);
		}
		else
			refCounts[ref]--;
	}
}

void Faust::RefManager::set_free_cb(void(*cb)(void*))
{
	this->cb = cb;
}
