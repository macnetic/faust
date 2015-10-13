#include "faust_core.h"
#include "faust_spmat.h"
#include "faust_mat.h"
#include "faust_vec.h"
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
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include "faust_constant.h"       /* time */
using namespace std;
vector<int> gen_unique_rand(int nbr, int range, bool UpperSide);


int main(int argc, char* argv[])
{

// debut gestion arguement (nom de fichier)	
	/*if(argc != 2 )
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
	*/


	int DIM;
	DIM=atoi(argv[1]);
	const int NB_RUN = 10000;
	const int NB_SAMPLE = 8;
	 vector<int> suite_nbr_factor={1,2,4,8};
	
	 
	int nnz;
	int MAX_RAND=9000;
	int nb_mat = 1;
	double density_max = 0.5;
	
	
	
	vector<vector<faust_core>> vec_core(NB_SAMPLE,vector<faust_core>(suite_nbr_factor.size()));
	vector<faust_mat> vec_mat(nb_mat);
	vector<int> global_id,id_row,id_col;
	vector<faust_real> values;
	faust_spmat spmat;
	faust_mat mat(DIM,DIM);
	faust_core f;
	

	cout<<"generation des faust_core"<<endl;
	std::srand(std::time(0));
	for (int id_factor=0;id_factor<suite_nbr_factor.size();id_factor++)
	{	
		
		int nbr_factor = suite_nbr_factor[id_factor];
		cout<<"nbr_factor :"<< nbr_factor<<endl;
		
		for (int i=1;i<=NB_SAMPLE;i++)
		{		
		
			f.clear();
		
			nnz =(int) ((double)(i*DIM*DIM/NB_SAMPLE/nbr_factor))*density_max;
			for (int k=0;k<nbr_factor;k++)
			{	
				values.resize(nnz),
				id_row.resize(nnz);
				id_col.resize(nnz);
				global_id =  gen_unique_rand(nnz,DIM*DIM,true);
		
				for (int j=0; j<nnz;j++)
				{
					id_row[j]=global_id[j]%DIM;
					id_col[j]=(int) (global_id[j]/DIM);
					values[j] = (faust_real)(std::rand()%MAX_RAND);			
				}
			//void init(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<faust_real>& values, const int dim1_, const int dim2_);
				spmat.init(id_row,id_col,values,DIM,DIM);
				//cout<<spmat.getNonZeros()<<endl;
				f.push_back(spmat);
			}
			
			//cout<<"ecriture faust_core " << f.size()<<" densite "<<(double)((double) f.get_total_nnz())/((double) f.getNbRow()*f.getNbCol())<<endl;
			vec_core[i-1][id_factor]=f;
			cout<<"ecriture faust_core " << vec_core[i-1][id_factor].size()<<" densite "<<(double)((double) vec_core[i-1][id_factor].get_total_nnz())/((double) vec_core[i-1][id_factor].getNbRow()*vec_core[i-1][id_factor].getNbCol())<<endl;
			cout<<"fin ecriture faust_core:"<< nbr_factor<<endl;
			f.clear();
		}
	}
	
	cout<<"GEN_DENSE_MAT "<<endl;
	for (int i=0;i<nb_mat;i++)
	{
		for (int k=0;k< (mat.getNbRow() * mat.getNbCol());k++)
		{	
		mat[k]=	std::rand()%MAX_RAND;
		}
		vec_mat[i]=mat;
	}

	for (int j=0;j<suite_nbr_factor.size();j++)
	{	
		for (int i=0;i<NB_SAMPLE;i++)
		{
			cout<<"FAUST_CORE "<<i<<j<<" : "<<endl;
			cout<<"nb fact : "<<vec_core[i][j].size()<<endl;
			cout<<"dim : "<<vec_core[i][j].getNbRow()<<"  "<<vec_core[i][j].getNbCol()<<endl;
			cout<<"nnz : "<<vec_core[i][j].get_total_nnz()<<endl;
			cout<<"density : "<<(double)((double) vec_core[i][j].get_total_nnz())/((double) vec_core[i][j].getNbRow()*vec_core[i][j].getNbCol())<<endl;
			cout<<endl<<endl;
		}
	}
	cout<<"FAUST_MAT_VEC "<<" : "<<endl;
	for (int i=0;i<nb_mat;i++)
	{
		cout<<"dim : "<<vec_mat[i].getNbRow()<<"  "<<vec_mat[i].getNbCol()<<endl;
	}

	
	faust_mat t_dense(NB_RUN,NB_SAMPLE);
	faust_mat t_faust(NB_RUN,NB_SAMPLE);
	
	t_dense.setZeros();
	t_faust.setZeros();
	vector<faust_mat> T_dense;
	vector<faust_mat> T_faust;
	for (int i=0;i<suite_nbr_factor.size();i++)
	{
		T_dense.push_back(t_dense);
		T_faust.push_back(t_faust);
	}
	
	faust_timer timer_dense_tmp;
	faust_timer timer_faust_tmp;

	for (int run=0 ; run < NB_RUN ; run++)
	{
		faust_timer t_run;
		cout << "run " << run+1 <<"/"<<NB_RUN<<endl;
		t_run.start();
		for (int j=0;j<suite_nbr_factor.size();j++)
		{	
			for (int i=0; i<NB_SAMPLE ; i++)
			{
				faust_vec x(DIM);
				for (int ll=0 ; ll<x.size() ;ll++)
					x[ll] = std::rand()*2.0/RAND_MAX-1.0;
				faust_vec y_dense(DIM);
				faust_vec y_faust(DIM);
				timer_dense_tmp.reset();
				timer_dense_tmp.start();
				y_dense = vec_mat[i%nb_mat] * x;
				timer_dense_tmp.stop();
				//t_dense[run][i] = timer_dense_tmp.get_time();
				T_dense[j].setCoeff(timer_dense_tmp.get_time(),run,i);

				timer_faust_tmp.reset();
				timer_faust_tmp.start();
				y_faust = vec_core[i][j] * x;
				timer_faust_tmp.stop();
			//t_faust[run][i] = timer_faust_tmp.get_time();
				T_faust[j].setCoeff(timer_faust_tmp.get_time(),run,i);

				
			}
		}
		t_run.stop();
		cout << "temps run = " << t_run.get_time() << " s" <<endl;
	}
	
	
	
	
	
	// debut ecriture des resultats dans fichier

	char filename_dense[150];
	char filename_faust[150];
	const char* path = "../output/";


	stringstream ss_outputfile;
	
	ss_outputfile.str("");
	ss_outputfile<<path;
	ss_outputfile << "time_comp_" <<"NB_FACT_";
	
	for (int i=0;i<suite_nbr_factor.size();i++)
	{
	ss_outputfile<<suite_nbr_factor[i]<<"_";	
	}	
	ss_outputfile<<"DIM_"<<DIM<<"_density_max_"<<density_max<< ".mat";
	
	cout<< "Writing results into " <<ss_outputfile.str()<< " file "<<endl;
	
	for (int i=0;i<suite_nbr_factor.size();i++)
	{
		stringstream ss_var_dense,ss_var_faust;
		ss_var_dense.str("");
		ss_var_faust.str("");
		ss_var_dense<<"t_dense_"<<suite_nbr_factor[i];
		ss_var_faust<<"t_fact_"<<suite_nbr_factor[i];
		
		write_faust_mat_into_matfile(T_dense[i],ss_outputfile.str().c_str(),ss_var_dense.str().c_str());
		write_faust_mat_into_matfile(T_faust[i],ss_outputfile.str().c_str(),ss_var_faust.str().c_str());
	}
	
	
	
	
	
}

vector<int> gen_unique_rand(int nbr, int range,bool UpperSide)
{	
	if (nbr>range)
	{
		cerr<<"gen_unique_rand : it's impossible to gen nbr  random different number in range if nbr > range"<<endl;
	}
	srand (time(NULL));
	int range_case = (int) range/nbr;
	vector<int> int_numbers(nbr);
	int lower_bound = 0;
	if (!UpperSide)
	{	
		for (int i=0; i< nbr ; i++)
		{
			int_numbers[i] = rand() % range_case + i*range_case;
		}
	}else
	{
		for (int i=0; i< nbr ; i++)
		{
			int_numbers[i] = range-1 - (rand() % range_case + i*range_case);
		}
	}	
	return int_numbers;
	
	
}



	

	
