#include "faust_core.h"
#include "faust_init_from_matio_core.h"
#include "faust_init_from_matio_mat.h"
#include "faust_timer.h"
#include "LinAlgebra.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>



// use this script to make runtime comparison between Faust-vector product and Dense matrix-vector product , this script takes to arguments :
// - the inputfile Faust_example.mat where all the Fausts with different Density, Dimensions and number of factor are stored. This file is generated by test/gen_artificial_faust.m 
// - NB_run the number of test
//
//The different computed times are stored in 2 different files :
// temps_faust.dat and temps_dense.dat
// Three headers files are generated that recaps the different Dimensions, theoretical RCGs and number of factors that are tested :
// - DIMS.dat
// - RCGS.dat
// - FACTORS.dat
 
using namespace std;
typedef FFPP faust_real; 

int main(int argc, char* argv[])
{



	

	if(argc < 2 )
	{
		cerr << "nombre d'arguments incorrect. 1 argument doit etre renseigne (nom du fichier (terminant par \".mat\") contenant les Faust et les matrices dense (fichier genere par test/gen_artificial_faust.m))"<<endl;
		
		cerr<<"second argument optionnel le nombre de RUN" << endl;
		
		exit(EXIT_FAILURE);
	}
	
	
	
	
	
	//string mat_file ="/udd/nbellot/Devel/Faust/trunk/devcpp/interface_matlab/Faust_example.mat";
	string mat_file(argv[1]);
	string prefix_output(argv[2]);
	size_t ind = mat_file.find_last_of(".");
	
	if(ind<=0 || ind>= mat_file.size())
	{
		cerr << "Le nom du fichier est incorrect" << endl;
		exit(EXIT_FAILURE);
	}
	string mat_file_extension(mat_file, ind);
	if(mat_file_extension.compare(".mat") != 0)
	{
		cerr << "Le nom du fichier doit se terminer par \".mat\"" << endl;
		exit(EXIT_FAILURE);
	}


	string mat_file_body_tmp(mat_file, 0, ind);


	string mat_file_body_dir, mat_file_body_file;

	ind = mat_file_body_tmp.find_last_of("/");
	if(ind<=0 || ind>= mat_file_body_tmp.size())
	{
		mat_file_body_dir = string("");
		mat_file_body_file = mat_file_body_tmp;
	}
	else
	{
		mat_file_body_dir = string(mat_file_body_tmp, 0, ind+1);
		mat_file_body_file = string(mat_file_body_tmp, ind+1);
	}
	
	int nb_run_tmp = 10;
	int id_file = 0;
	if (argc >= 4) nb_run_tmp = atoi(argv[3]);
	if (argc >= 5) id_file = atoi(argv[4]);
	
	const int NB_RUN = nb_run_tmp;
	
	int nb_fact;
	int RCG;
	int DIM;
	
	faust_mat<faust_real> DIMS;
	faust_mat<faust_real> RCGS;
	faust_mat<faust_real> NB_FACTS;
	faust_mat<faust_real> Dense_mat;
	faust_core<faust_real> fc;

	
	//void init_faust_mat_from_matio(faust_mat& M, const char* fileName, const char* variableName);
	init_faust_mat_from_matio(DIMS,mat_file.c_str(),"Dims");
	init_faust_mat_from_matio(RCGS,mat_file.c_str(),"RCGs");
	init_faust_mat_from_matio(NB_FACTS,mat_file.c_str(),"nb_facts");
	
	
	
	
	
	DIMS.Display();
	RCGS.Display();
	
	int nDIMS=DIMS.getNbRow()*DIMS.getNbCol();
	int nRCGS=RCGS.getNbRow()*RCGS.getNbCol();
	int nNB_FACTS=NB_FACTS.getNbRow()*NB_FACTS.getNbCol();
	
	vector<faust_mat<faust_real> > Dense_matS(nDIMS);
	//vector<vector<faust_core> > FcoreS(nRCGS,vector<faust_core>(nDIMS));
	faust_core<faust_real> FcoreS[nNB_FACTS][nRCGS][nDIMS];
	
	stringstream ssInputVarFaust;
	stringstream ssInputVarDense;
	string OuputDir="../output/";
	
	
	cout << "loading of dense matrix en faust" << endl;
	
	for (int k=0;k<nNB_FACTS;k++)
	{
		nb_fact = NB_FACTS[k];
		cout<<"  NB_FACT : "<<nb_fact<<endl;	
		for (int j=0;j<nDIMS;j++)
		{
		
		
			DIM = DIMS[j];
			cout<<"  DIM : "<<DIM<<endl;
			ssInputVarDense.str("");
			ssInputVarDense<<"DMat_Dim_"<<DIM;
			init_faust_mat_from_matio(Dense_mat,mat_file.c_str(),ssInputVarDense.str().c_str());
		
			Dense_matS[j]=Dense_mat;
		
			cout<<"      RCG : ";
			for (int i=0;i<nRCGS;i++)
			{
				RCG = RCGS[i];
				cout<<RCG<<" ";
				ssInputVarFaust.str("");
				ssInputVarFaust<<"Faust_nfact_"<<nb_fact<<"_RCG_"<<RCG<<"_Dim_"<<DIM;
				//init_faust_core_from_matiofile(faust_core& core, const char* fileName, const char* variableName)
				init_faust_core_from_matiofile(fc,mat_file.c_str(),ssInputVarFaust.str().c_str());
				FcoreS[k][i][j]=fc;
			}
			cout<<endl;
		}
	}	
	cout<<endl<<endl<<endl;
	
	
	// vector<vector<vector<float> > t_dense(vector<nRCGS, vector<float>(nDIMS));
	// vector<vector<float> > t_faust(nRCGS, vector<float>(nDIMS));
	
	
	float t_dense[NB_RUN][nNB_FACTS][nRCGS][nDIMS];
	float t_faust[NB_RUN][nNB_FACTS][nRCGS][nDIMS];
	
	faust_timer timer_faust;
	faust_timer timer_dense;




	cout << "calculus computing" << endl;
	faust_vec<faust_real> x_tmp(Dense_matS[0].getNbCol());
	for (int j=0 ; j<x_tmp.size() ;j++)
		x_tmp[j] = std::rand()*2.0/RAND_MAX-1.0;
	faust_vec<faust_real> y_dense_tmp(Dense_matS[0].getNbRow());
	y_dense_tmp = Dense_matS[0] * x_tmp;
	for (int k=0;k<NB_RUN;k++)
	{
		cout<<"RUN "<<k<<endl;
		for (int n=0;n<nNB_FACTS;n++)
		{
		
			for (int i=0;i<nRCGS;i++)
			{	
				RCG = RCGS[i];
				for (int j=0;j<nDIMS;j++)
				{
					DIM = DIMS[j];
					faust_vec<faust_real> x(DIM);
					Dense_mat = Dense_matS[j];
					faust_core<faust_real> fc1(FcoreS[n][i][j]);
					for (int ii=0;ii<DIM;ii++)x[ii]=std::rand()*2.0/RAND_MAX-1;
					faust_vec<faust_real> y_dense(DIM);
					faust_vec<faust_real> y_faust(DIM);

					timer_dense.reset();
					timer_dense.start();
					y_dense = Dense_mat * x;
					timer_dense.stop();
			
					timer_faust.reset();
					timer_faust.start();
					y_faust = fc1 * x;
					timer_faust.stop();
			
					t_faust[k][n][i][j] = timer_faust.get_time();
					t_dense[k][n][i][j] = timer_dense.get_time();
			
				}
			}
		}
	}

	ofstream fichier;
	stringstream ss;
	
	ss.str("");
	ss <<OuputDir<<prefix_output<<"temps_dense"<<id_file<<".dat";
	cout << "fichier output 1 : "<<ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int k=0 ; k <NB_RUN;k++)
	{	
		for (int n=0;n<nNB_FACTS;n++)
		{
			for (int j=0 ; j<nDIMS ; j++)
			{
				for (int i=0 ; i < nRCGS ; i++)
				{
					fichier << setprecision(20) << t_dense[k][n][i][j] << " ";
				}
			}
		}	
		fichier << endl;
	}
	
	fichier.close();

	ss.str("");
	ss <<OuputDir<<prefix_output<<"temps_faust"<<id_file<<".dat";
	cout <<"fichier ouput2  : "<<ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int k=0 ; k <NB_RUN;k++)
	{
		for (int n=0;n<nNB_FACTS;n++)
		{
			for (int j=0 ; j<nDIMS; j++)
			{
				for (int i=0 ; i < nRCGS ; i++)
				{
					fichier << setprecision(20) << t_faust[k][n][i][j] << " ";
				}
			}
		}
		fichier << endl;
	}
	fichier.close();
	
	
	//ecriture des headers
	ss.str("");
	ss <<OuputDir<<prefix_output<<"DIMS"<<id_file<<".dat";
	cout <<"fichier header DIMS  : "<<ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int j=0 ; j<nDIMS; j++)fichier << setprecision(20) << DIMS[j] << " ";
	fichier.close();
	
	
	ss.str("");
	ss <<OuputDir<<prefix_output<<"RCGS"<<id_file<<".dat";
	cout <<"fichier header RCGS  : "<<ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int j=0 ; j<nRCGS; j++)fichier << setprecision(20) << RCGS[j] << " ";
	fichier.close();

	
	ss.str("");
	ss <<OuputDir<<prefix_output<<"NB_FACTS"<<id_file<<".dat";
	cout <<"fichier header NB_FACTS  : "<<ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int j=0 ; j<nNB_FACTS; j++)
	fichier << setprecision(20) << NB_FACTS[j] << " ";
	fichier.close();
	
	
	
	return 0;
   
}
