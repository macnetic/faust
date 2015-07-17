#include "faust_core.h"
#include "faust_init_from_matio_core.h"
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

// debut gestion arguement (nom de fichier)	
	if(argc != 2 )
	{
		cerr << "nombre d'arguments incorrect. Un argument doit etre renseigne (nom du fichier (terminant par \".mat\") contenant la cellule mat_cells)" << endl;
		exit(EXIT_FAILURE);
	}
	string mat_file(argv[1]);
	
	size_t ind = mat_file.find_last_of(".");
	
	string mat_file_body(mat_file, 0, ind);
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
// fin gestion arguement (nom de fichier)	




	const int NB_RUN = 1000;

	vector<faust_core>* vec_core = new vector<faust_core>();
	vector<faust_mat>* vec_dense_mat = new vector<faust_mat>();
	init_faust_data_from_matiofile(*vec_dense_mat, *vec_core, argv[1], "mat_cells");

	const vector<faust_core> core(*vec_core);
	const vector<faust_mat> dense_mat(*vec_dense_mat);

	delete vec_core; vec_core=NULL;
	delete vec_dense_mat; vec_dense_mat=NULL;

	
	if(core.size()!=dense_mat.size())
	{
		cerr << "probleme de chargement des matrices a partir de Matlab"<<endl;
		exit(EXIT_FAILURE);
	}




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





	std::srand(std::time(0));

	vector<vector<float> > t_dense(NB_RUN, vector<float>(core.size()));
	vector<vector<float> > t_faust(NB_RUN, vector<float>(core.size()));

	faust_timer timer_dense_tmp;
	faust_timer timer_faust_tmp;

	for (int run=0 ; run < NB_RUN ; run++)
	{
		faust_timer t_run;
		cout << "run " << run+1 <<"/"<<NB_RUN<<endl;
		t_run.start();
		for (int i=core.size()-1 ; i>=0 ; i--)
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
			cout<<"err relative = " << err_rel << endl;
		}
		t_run.stop();
		cout << "temps run = " << t_run.get_time() << " s" <<endl;
	}




	
// debut ecriture des resultats dans fichier
	char filename_dense[150];
	char filename_faust[150];

	ofstream fichier;
	stringstream ss;
	
	ss.str("");
	ss << "temps_dense_" << mat_file_body.c_str() << ".dat";
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
	ss << "temps_faust_" << mat_file_body.c_str() << ".dat";
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

   
}
