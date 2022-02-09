#include "faust_RefManager.h"
#include <exception>

using namespace std;
void Faust::RefManager::acquire(void* ref)
{
#ifdef DEBUG
	cout << "Faust::RefManager::acquire ref: " << ref << endl;
#endif
	if(refCounts.find(ref) != refCounts.end())
	{
		refCounts[ref]++;
#ifdef DEBUG
		cout << "Faust::RefManager::acquire refCount=" << refCounts[ref] << endl;
#endif
	}
	else
	{
		refCounts[ref] = 1;
#ifdef DEBUG
		cout << "Faust::RefManager::acquire refCount=" << refCounts[ref] << endl;
#endif
	refCounts.insert(make_pair(ref,1));
	}
}

void Faust::RefManager::release(void* ref)
{
#ifdef DEBUG
	cout << "Faust::RefManager::release" << endl;
#endif
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
	else
		throw runtime_error("RefManager error: unknown pointer reference asked.");
}

void Faust::RefManager::set_free_cb(void(*cb)(void*))
{
	this->cb = cb;
}

bool Faust::RefManager::contains(void* ref) const
{

	return refCounts.find(ref) != refCounts.end();
}

Faust::RefManager::RefManager()
{
}

Faust::RefManager::RefManager(void(*cb)(void*))
{
	set_free_cb(cb);
}

void Faust::RefManager::print_refs()
{
	std::cout << "number of refs managed: " << refCounts.size() << std::endl;
	for(map<void*, unsigned int>::const_iterator it = refCounts.begin();
			it != refCounts.end(); ++it)
	{
		std::cout << "address: " << it->first << " ref counts:" << it->second << std::endl;
	}
}
