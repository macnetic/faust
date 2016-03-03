#include <stdlib.h>
#include "faust_params.h"
#include <iostream>
#include <vector>
#include <string>
#include "faust_init_params_from_xml.h"
#include "hierarchical_fact.h"
#include "faust_core.h"

/// Definition of Floating Point Precision
typedef double FPP; 
using namespace std;

/*! \brief Compute the hierarchical factorization of a given data matrix A in cmdline mode.<br>
* Projet name is "launch_hierarchical_fact". It is available in the /wrapper/cmd_line/src/*.cpp <br> 
* \param fileName1 : An xml configuration file specifying the different parameters of the hierarchical_fact algorithm <br>
  \param fileName2 : The data text file where the matrix A that will be factorized is stored. <br> 
*		The first line of the file contains 2 integer : the number of row and the number of column.<br>
*		All the other line contains one coefficient in ColMajor access of the matrix).<br>
* \param character : (Optionnal) Character specifying if the matrix stored correspond to the data matrix or its transposed <br>
*	'N' (default value) --> the two matrix are the same <br>
*       'T'                 --> the data matrix is the transposed of the matrix stored in the text file.
*/
int main(int argc, char* argv[]) 
{
	faust_params<FPP> params;
	cout<<"argc : "<<argc<<endl; 
	if (argc < 3)
	{
		cerr << "incorrect number of argument (at least to argument must be specified) : 1st argument is the data filename\n second is the configuration filename  " << endl;
		exit(EXIT_FAILURE);
	}

	string config_filename(argv[1]); // filename of xml config hierarchial_fact file
	string data_filename(argv[2]); // filename of the data (matrix) filename

	char operator_data='N';
	if (argc >= 4)
	{	
		
		char* operator_data_tmp(argv[3]);
		if ((strcmp(operator_data_tmp,"N")!=0) && (strcmp(operator_data_tmp,"T")!=0))
		{
			cerr << " the optionnal third argument must be 'T' if the data is the transpose of the matrix stored in  data filename else it is 'N'" << endl;
			exit(EXIT_FAILURE);
		}	
		operator_data = operator_data_tmp[0];
	}
		
	size_t ind = config_filename.find_last_of(".");
	if(ind<=0 || ind>= config_filename.size())
	{
		cerr << "incorrect configuration filename" << endl;
		exit(EXIT_FAILURE);
	}
	
	string config_filename_extension(config_filename, ind);
	if(config_filename_extension.compare(".xml") != 0)
	{
		cerr << "extension of the configuration file must be  \".xml\"" << endl;
		exit(EXIT_FAILURE);
	}
	string config_filename_body(config_filename, 0, ind); 
	string outputfilename = config_filename_body + "_FAUST.txt";
	
	
	// initialization 
	init_params_from_xml(config_filename.c_str(),params);
	params.data.init_from_file(data_filename.c_str());
	if (operator_data=='T')
		params.data.transpose();
	
	params.check_constraint_validity();
	params.Display();
	hierarchical_fact<FPP> hier_fact(params);
	
	hier_fact.compute_facts();
	
	cout<<"lambda="<<std::setprecision(20)<<hier_fact.get_lambda()<<endl;
	
	
	faust_core<FPP> faust_facts;
	hier_fact.get_facts(faust_facts);
	faust_facts.scalarMultiply(hier_fact.get_lambda());
	std::cout<<"writing factorization into "<< outputfilename <<endl;
	faust_facts.print_file(outputfilename.c_str());
		
	
}
