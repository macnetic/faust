#include <stdlib.h>
#include "faust_params.h"
#include <iostream>
#include <vector>
#include<string>
#include "faust_init_params_from_xml.h"
#include "hierarchical_fact.h"


typedef float FPP;
using namespace std;

int main(int argc, char* argv[]) 
{
	
	faust_params<FPP> params;
	string config_filename(argv[1]);
	init_params_from_xml(config_filename.c_str(),params);
	hierarchical_fact<FPP> hier_fact(params);
	
	
	hier_fact.compute_facts();
	cout<<"lambda="<<std::setprecision(20)<<hier_fact.get_lambda()<<endl;
	
}