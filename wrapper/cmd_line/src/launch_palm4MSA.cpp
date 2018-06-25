/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#include <stdlib.h>
#include "faust_ParamsPalm.h"
#include <iostream>
#include <vector>
#include<string>
#include "faust_init_params_from_xml.h"
#include "faust_Palm4MSA.h"
#include "faust_Transform.h"

/// Definition of Floating Point Precision
typedef float FPP;
using namespace std;

/*! \brief Palm4MSA consists of "Proximal Alternating Linearized Minimization for Multi-layer Sparse Approximation". <br>
* This command-line Palm4MSA factorizes an input matrix corresponding to an operator of interest into J sparse factors and converges to a stationary point. <br>
* Projet name is launch_Palm4MSA. It is available in the ./wrapper/cmd_line/src/ <br>
* \param fileName1 : An xml configuration file specifying the different parameters of the Palm4MSA algorithm <br>
  \param fileName2 : The data text file where the matrix A that will be factorized is stored. <br>
*/
int main(int argc, char* argv[])
{

	if (argc < 2)
	{
		cerr << "incorrect number of argument (at least to argument must be specified) : 1st argument is the data filename\n second is the configuration filename  " << endl;
		exit(EXIT_FAILURE);
	}

	Faust::ParamsPalm<FPP,Cpu> params;


	string config_filename(argv[1]);
	string data_filename(argv[2]); // filename of the data (matrix) filename

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
	string outputfilename = config_filename_body + "_FAUST.txt";


	init_palm_params_from_xml(config_filename.c_str(),params);
	params.data.init_from_file(data_filename.c_str());
	params.check_constraint_validity();
	cout<<"params initialisé"<<endl;
	params.Display();
	Faust::BlasHandle<Cpu> blas_handle;
	Faust::Palm4MSA<FPP,Cpu> palm(params,blas_handle);


	palm.compute_facts();
	cout<<"RMSE : "<<palm.get_RMSE()<<endl;
	cout<<"lambda="<<std::setprecision(20)<<palm.get_lambda()<<endl;

	Faust::Transform<FPP,Cpu> faust_facts;
	palm.get_facts(faust_facts);
	std::cout<<"faust_fact size : "<<faust_facts.size()<<endl;
	faust_facts.scalarMultiply(palm.get_lambda());
	std::cout<<"writing factorization into "<< outputfilename <<endl;
	faust_facts.print_file(outputfilename.c_str());

}
