#include <map>
#include <iostream>

#ifndef __FAUST_REF_MANAGER__
#define __FAUST_REF_MANAGER__
namespace Faust {
	class RefManager {
		std::map<void*,unsigned int> refCounts;
		void(*cb)(void*);

		public:
		void acquire(void* ref);

		void release(void* ref);

		bool contains(void* ref) const;

		void set_free_cb(void (*cb)(void*));

		/**
		 * Prints all managed addresses with their references count.
		 */
		void print_refs();

		RefManager();
		RefManager(void(*cb)(void*));
	};
}
#endif

