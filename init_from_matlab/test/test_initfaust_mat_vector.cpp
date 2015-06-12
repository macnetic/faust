#include "faust_mat.h"
#include <cstring>
#include <sstream>
#include <iostream>
#include <vector>
#include "faust_init_from_matio.h"


using namespace std;

int main(int argc, char* argv[])
{	
	string fileName = "cas_test3.mat";
	string variableName = "facts";
	vector<faust_mat> facts;
	
	
	cout<<"avant init"<<endl;
	init_faust_mat_vector_from_matiofile(facts,fileName.c_str(),variableName.c_str());
	cout<<"apres init"<<endl;

	for (int i=0;i<facts.size();i++)
	{
		facts[i].Display();
	}
	
	
	
}