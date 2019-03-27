#include <map>
#include <iostream>
using namespace std;

#ifndef __FAUST_REF_MANAGER__
#define __FAUST_REF_MANAGER__
namespace Faust {
	class RefManager {
		map<void*,unsigned int> refCounts;
		typedef void (RefManager::*free_cb)(void*);
		void(*cb)(void*);

		public:
		void acquire(void* ref);

		void release(void* ref);

		void set_free_cb(void (*cb)(void*));
	};
}
#endif

