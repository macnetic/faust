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


using namespace std;

int main(int argc, char* argv[])
{



	



	const int NB_RUN = 10;

	
	
	
	string mat_file ="/udd/nbellot/Devel/Faust/trunk/devcpp/interface_matlab/Faust_example.mat";
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
	
	int nb_fact = 2;
	int RCG = 2;
	int DIM = 512;
	
	faust_mat DIMS;
	faust_mat RCGS;
	faust_mat Dense_mat;
	faust_core fc;

	
	//void init_faust_mat_from_matio(faust_mat& M, const char* fileName, const char* variableName);
	init_faust_mat_from_matio(DIMS,mat_file.c_str(),"Dims");
	init_faust_mat_from_matio(RCGS,mat_file.c_str(),"RCGs");

	
	
	
	
	DIMS.Display();
	RCGS.Display();
	
	int nDIMS=DIMS.getNbRow()*DIMS.getNbCol();
	int nRCGS=RCGS.getNbRow()*RCGS.getNbCol();
	
	vector<faust_mat> Dense_matS(nDIMS);
	vector<vector<faust_core> > FcoreS(nRCGS,vector<faust_core>(nDIMS));
	
	stringstream ssInputVarFaust;
	stringstream ssInputVarDense;
	cout << "loading of dense matrix en faust" << endl;
	for (int j=0;j<nDIMS;j++)
	{
		
		
		DIM = DIMS[j];
		cout<<"DIM : "<<DIM<<endl;
		ssInputVarDense.str("");
		ssInputVarDense<<"DMat_Dim_"<<DIM;
		init_faust_mat_from_matio(Dense_mat,mat_file.c_str(),ssInputVarDense.str().c_str());
		
		Dense_matS[j]=Dense_mat;
		
		cout<<"  RCG : ";
		for (int i=0;i<nRCGS;i++)
		{
			RCG = RCGS[i];
			cout<<RCG<<" ";
			ssInputVarFaust.str("");
			ssInputVarFaust<<"Faust_nfact_"<<nb_fact<<"_RCG_"<<RCG<<"_Dim_"<<DIM;
			//init_faust_core_from_matiofile(faust_core& core, const char* fileName, const char* variableName)
			init_faust_core_from_matiofile(fc,mat_file.c_str(),ssInputVarFaust.str().c_str());
			FcoreS[i][j]=fc;
		}
		cout<<endl;
		
	}	
	cout<<endl<<endl<<endl;
	
	
	// vector<vector<vector<float> > t_dense(vector<nRCGS, vector<float>(nDIMS));
	// vector<vector<float> > t_faust(nRCGS, vector<float>(nDIMS));
	
	
	float t_dense[NB_RUN][nRCGS][nDIMS];
	float t_faust[NB_RUN][nRCGS][nDIMS];
	
	faust_timer timer_faust;
	faust_timer timer_dense;




	cout << "calculus computing" << endl;
	faust_vec x_tmp(Dense_matS[0].getNbCol());
	for (int j=0 ; j<x_tmp.size() ;j++)
		x_tmp[j] = std::rand()*2.0/RAND_MAX-1.0;
	faust_vec y_dense_tmp(Dense_matS[0].getNbRow());
	y_dense_tmp = Dense_matS[0] * x_tmp;
for (int k=0;k<NB_RUN;k++)
{
	cout<<"RUN "<<k<<endl;
	for (int i=0;i<nRCGS;i++)
	{	
		RCG = RCGS[i];
		for (int j=0;j<nDIMS;j++)
		{
			DIM = DIMS[j];
			faust_vec x(DIM);
			Dense_mat = Dense_matS[j];
			faust_core fc1(FcoreS[i][j]);
			for (int ii=0;ii<DIM;ii++)x[ii]=std::rand()*2.0/RAND_MAX-1;
			faust_vec y_dense(DIM);
			faust_vec y_faust(DIM);

			timer_dense.reset();
			timer_dense.start();
			y_dense = Dense_mat * x;
			timer_dense.stop();
			
			timer_faust.reset();
			timer_faust.start();
			y_faust = fc1 * x;
			timer_faust.stop();
			
			t_faust[k][i][j] = timer_faust.get_time();
			t_dense[k][i][j] = timer_dense.get_time();
			
		}
	}
}

	ofstream fichier;
	stringstream ss;
	
	ss.str("");
	ss << mat_file_body_dir <<"temps_dense_" << mat_file_body_file<< ".dat";
	cout << "fichier output 1 : "<<ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int k=0 ; k <NB_RUN;k++)
	{	
		for (int j=0 ; j<nDIMS ; j++)
		{
			for (int i=0 ; i < nRCGS ; i++)
			{
				fichier << setprecision(20) << t_dense[k][i][j] << " ";
			}
		}
		fichier << endl;
	}
	
	fichier.close();

	ss.str("");
	ss << mat_file_body_dir <<"temps_faust_" << mat_file_body_file.c_str() << ".dat";
	cout <<"fichier ouput2  : "<<ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int k=0 ; k <NB_RUN;k++)
	{
		for (int j=0 ; j<nDIMS; j++)
		{
			for (int i=0 ; i < nRCGS ; i++)
			{
				fichier << setprecision(20) << t_faust[k][i][j] << " ";
			}
		}
		fichier << endl;
	}
	fichier.close();

	
	
/*
//debut affichage des matrices
	for(int i=0 ; i<core.size() ;i++)
	{
		cout<<"core["<<i<<"] contains " << core[i].size() << " sparse matrices : "<<endl<<endl;
		core[i].Display();
		cout<<endl;
	}	
        cout<<endl;
	for(int i=0 ; i<dense_mat.size() ; i++)
	{
		if(dense_mat[i].getNbCol()>20)
			cout << "dense_mat[" << i << "] = " << dense_mat[i].getNbRow() <<"x"<< dense_mat[i].getNbCol()<<" dense matrix" << endl;
		else
		{
			cout << "dense_mat[" << i << "] = " << endl;
			dense_mat[i].Display();
			cout<<endl;
		}
	}	
//fin affichage des matrices
*/

/*


	std::srand(std::time(0));

	vector<vector<float> > t_mean_faust(NB_RUN, vector<float>(core.size()));
	vector<vector<float> > t_faust(NB_RUN, vector<float>(core.size()));

	faust_timer timer_dense_tmp;
	faust_timer timer_faust_tmp;



	faust_vec x_tmp(dense_mat[0].getNbCol());
	for (int j=0 ; j<x_tmp.size() ;j++)
		x_tmp[j] = std::rand()*2.0/RAND_MAX-1.0;
	faust_vec y_dense_tmp(dense_mat[0].getNbRow());
	y_dense_tmp = dense_mat[0] * x_tmp;


	for (int run=0 ; run < NB_RUN ; run++)
	{
		faust_timer t_run;
		//cout << "run " << run+1 <<"/"<<NB_RUN<<endl;
		t_run.start();
		for (int i=core.size()-1 ; i>=0 ; i--)
		{
			if (status[i]==0)
			{
				faust_vec x(dense_mat[i].getNbCol());
				for (int j=0 ; j<x.size() ;j++)
					x[j] = std::rand()*2.0/RAND_MAX-1.0;
				faust_vec y_dense(dense_mat[i].getNbRow());
				faust_vec y_faust(dense_mat[i].getNbRow());
				
				timer_dense_tmp.reset();
				timer_dense_tmp.start();
				y_dense = dense_mat[i] * x;
				timer_dense_tmp.stop();
				t_dense[run][i] = timer_dense_tmp.get_time();
	
				timer_faust_tmp.reset();
				timer_faust_tmp.start();
				y_faust = core[i] * x;
				timer_faust_tmp.stop();
				t_faust[run][i] = timer_faust_tmp.get_time();
	
				faust_real err_rel = y_faust.mean_relative_error(y_dense);
				//cout<<"err relative = " << err_rel << endl;
			}
			else
			{
				t_dense[run][i] = nan("");
				t_faust[run][i] = nan("");
			}
		}
		t_run.stop();
		//cout << "temps run = " << t_run.get_time() << " s" <<endl;
	}




	
// debut ecriture des resultats dans fichier
	char filename_dense[150];
	char filename_faust[150];

	ofstream fichier;
	stringstream ss;
	
	ss.str("");
	ss << mat_file_body_dir <<"temps_dense_" << mat_file_body_file.c_str() << ".dat";
	cout << ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int run=0 ; run < NB_RUN ; run++)
	{
		for (int i=0 ; i<core.size() ; i++)
		{
			fichier << setprecision(20) << t_dense[run][i] << " ";
		}
		fichier << endl;
	}
	fichier.close();

	ss.str("");
	ss << mat_file_body_dir <<"temps_faust_" << mat_file_body_file.c_str() << ".dat";
	cout << ss.str().c_str() <<endl;
	fichier.open(ss.str().c_str());
	for (int run=0 ; run < NB_RUN ; run++)
	{
		for (int i=0 ; i<core.size() ; i++)
		{
			fichier << setprecision(15) << t_faust[run][i] << " ";
		}
		fichier << endl;
	}
	fichier.close();
// fin ecriture des resultats dans fichier
*/

	return 0;
   
}
