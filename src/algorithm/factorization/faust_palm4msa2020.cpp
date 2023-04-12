#include <string>
#include <iostream>
#include "faust_palm4msa2020.h"

using namespace std;

namespace Faust {

	volatile bool palm4msa2_interrupted;

	void palm4msa2_signal_handler(int signal)
	{
		if(signal == SIGINT)
		{
			cerr << string(15, '=') << " INTERRUPT PALM4MSA" << endl;
			palm4msa2_interrupted = true;
		}
	}

	bool is_palm4msa2_interrupted()
	{
		return palm4msa2_interrupted;
	}

	void init_palm4msa2_interrupt()
	{
		palm4msa2_interrupted = false;
	}
}
