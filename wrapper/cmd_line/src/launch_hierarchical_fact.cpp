#include "faust_MatSparse.h"
#include "faust_HierarchicalFact.h"
#include "faust_Timer.h"
#include "faust_Transform.h"
//#include "faust_init_from_matio_params.h"
//#include "faust_init_from_matio_core.h"
#include <string>
#include <sstream>
#include "faust_BlasHandle.h"
#include "faust_SpBlasHandle.h"
#include "faust_init_params_from_xml.h"

#include "faust_Timer.h"


/// Definition of Floating Point Precision
typedef double FPP;
using namespace std;


/*! \brief Compute the hierarchical factorization of a given data matrix A in cmdline mode.<br>
* Projet name is "launch_hierarchical_fact". It is available in the /wrapper/cmd_line/src/*.cpp <br>
* \param fileName1 : An xml configuration file specifying the different parameters of the HierarchicalFact algorithm <br>
  \param fileName2 : The data text file where the matrix A that will be factorized is stored. <br>
*		The first line of the file contains 2 integer : the number of row and the number of column.<br>
*		All the other line contains one coefficient in ColMajor access of the matrix).<br>
* \param character : (Optionnal) Character specifying if the matrix stored correspond to the data matrix or its transposed <br>
*	'N' (default value) --> the two matrix are the same <br>
*       'T'                 --> the data matrix is the transposed of the matrix stored in the text file. <br>
*
* \param nbiter_comp_time : (Optionnal) number of iteration for time comparison betwwen matrix-vector product and faust-vector product  <br>
* default value 0 (no time comparison)
*/
int main(int argc, char* argv[])
{
	Faust::Params<FPP,Cpu> params;
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

	// time comparison between matrix-vector product
	int niter_time_comp = 0;
	if (argc >= 5)
	{
	 	niter_time_comp=atoi(argv[4]);
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
	Faust::MatDense<FPP,Cpu> data_matrix;
	data_matrix.init_from_file(data_filename.c_str());
	if (operator_data=='T')
		data_matrix.transpose();
	params.data=data_matrix;

	params.check_constraint_validity();
	Faust::BlasHandle<Cpu> blas_handle;
	Faust::SpBlasHandle<Cpu> spblas_handle;
	std::cout<<"**************** PARAMETER OF HIERARCHICAL_FACT **************** "<<std::endl;
	params.Display();
	Faust::HierarchicalFact<FPP,Cpu> hier_fact(params,blas_handle,spblas_handle);

	std::cout<<"****************  FACTORIZATION **************** "<<std::endl;
	hier_fact.compute_facts();

	cout<<"lambda="<<std::setprecision(20)<<hier_fact.get_lambda()<<endl;


	Faust::Transform<FPP,Cpu> faust_facts;
	hier_fact.get_facts(faust_facts);
	faust_facts.scalarMultiply(hier_fact.get_lambda());


	std::cout<<"writing faust (factorisation) into  "<< outputfilename <<endl;
	faust_facts.print_file(outputfilename.c_str());
	std::cout<<std::endl;
		std::cout<<"**************** RELATIVE ERROR BETWEEN FAUST AND DATA MATRIX **************** "<<std::endl;

	//relative_error
	Faust::MatDense<FPP,Cpu> faust_product;
	faust_product=faust_facts.get_product();
	faust_product-=data_matrix;
	FPP relative_error = faust_product.norm()/data_matrix.norm();

	std::cout<<"		"<<relative_error<<std::endl<<std::endl;

	//time comparison between matrix vector product and faust-vector product
	if (niter_time_comp > 0)
	{

		Faust::Timer tdense;
		Faust::Timer tfaust;
		Faust::Vect<FPP,Cpu> x(data_matrix.getNbCol());
		Faust::Vect<FPP,Cpu> ydense(data_matrix.getNbRow());
		Faust::Vect<FPP,Cpu> yfaust(faust_facts.getNbRow());
		for (int i=0;i<niter_time_comp;i++)
		{
			//random initilisation of vector x
	 		for (int j=0;j<x.size();j++)
			{
				x[j]=std::rand()*2.0/RAND_MAX-1.0;
			}

			tdense.start();
				ydense = data_matrix * x;
			tdense.stop();

	 		tfaust.start();
				yfaust = faust_facts * x;
			tfaust.stop();

	 	}
		std::cout<<std::endl;


		std::cout<<"**************** TIME COMPARISON MATRIX VECTOR PRODUCT **************** "<<std::endl;
		std::cout<<"	TIME  SPEED-UP : "<<tdense.get_time()/tfaust.get_time()<<std::endl;
		std::cout<<"	MEAN TIME dense : "<<tdense.get_time()/((float) niter_time_comp)<<std::endl;
		std::cout<<"	MEAN TIME faust : "<<tfaust.get_time()/((float) niter_time_comp)<<std::endl;
	}




}
