#include <stdlib.h>
#include "faust_params_palm.h"
#include <iostream>
#include <vector>
#include<string>
#include "faust_init_params_from_xml.h"
#include "palm4MSA.h"
#include "faust_core.h"


typedef float faust_real;
using namespace std;

int main(int argc, char* argv[]) 
{
	
	faust_params_palm<faust_real> params;
	string config_filename(argv[1]);
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
	
	
	init_palm_params_from_xml(config_filename.c_str(),params);
	cout<<"params initialisé"<<endl;
	params.Display();
	palm4MSA<faust_real> palm(params);
	
	
	palm.compute_facts();
	cout<<"RMSE : "<<palm.get_RMSE()<<endl;
	cout<<"lambda="<<std::setprecision(20)<<palm.get_lambda()<<endl;
	
	faust_core<faust_real> faust_facts;
	palm.get_facts(faust_facts);
	std::cout<<"faust_fact size : "<<faust_facts.size()<<endl;
	faust_facts.scalarMultiply(palm.get_lambda());
	std::cout<<"writing factorization into "<< outputfilename <<endl;
	faust_facts.print_file(outputfilename.c_str());
	
}
