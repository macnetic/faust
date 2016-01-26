#include <stdlib.h>
#include "faust_params.h"
#include <iostream>
#include <vector>
#include<string>
#include "faust_init_params_from_xml.h"
#include "hierarchical_fact.h"
#include "faust_core.h"


typedef double FPP;
using namespace std;

int main(int argc, char* argv[]) 
{
	
	faust_params<FPP> params;
	string config_filename(argv[1]); // filename of xml config hierarchial_fact file
	size_t ind = config_filename.find_last_of(".");
	if(ind<=0 || ind>= config_filename.size())
	{
		cerr << "Le nom du fichier est incorrect" << endl;
		exit(EXIT_FAILURE);
	}
	
	string config_filename_extension(config_filename, ind);
	if(config_filename_extension.compare(".xml") != 0)
	{
		cerr << "Le nom du fichier doit se terminer par \".xml\"" << endl;
		exit(EXIT_FAILURE);
	}
	string config_filename_body(config_filename, 0, ind); 
	string outputfilename = config_filename_body + "_FACT.txt";
	
	
	
	init_params_from_xml(config_filename.c_str(),params);
	hierarchical_fact<FPP> hier_fact(params);
	hier_fact.compute_facts();
	
	
	cout<<"lambda="<<std::setprecision(20)<<hier_fact.get_lambda()<<endl;
	
	
	faust_core<FPP> faust_facts;
	hier_fact.get_facts(faust_facts);
	faust_facts.scalarMultiply(hier_fact.get_lambda());
	std::cout<<"writing factorization into "<< outputfilename <<endl;
	faust_facts.print_file(outputfilename.c_str());
	
	
	
}