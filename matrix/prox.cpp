#include "prox.h"
#include <vector>
#include "faust_constant.h"
#include <iostream>

void prox_sp(faust_mat & M,int k)
{	
	int dim1 = M.getNbRow();
	int dim2 = M.getNbCol();
	if (k>dim1*dim2)
	{
		std::cerr << "ERROR prox_sp : k > number of elements of M " << std::endl;
		exit( EXIT_FAILURE);	
	}
	faust_mat new_M(dim1,dim2);
	faust_mat M_abs(dim1,dim2);
	M_abs = M;
	M_abs.abs();
	new_M.setZeros();
	
	/*
	std::cout<<"*** M ****"<<std::endl;
	M.Display();
	
	
	std::cout<<"*** M_abs ****"<<std::endl;
	M_abs.Display();
	
	
	std::cout<<"*** new_M ****"<<std::endl;
	new_M.Display();
	*/
	
	int nb_elt_found = 0;



	std::vector<int> current_id_row;
	std::vector<int> current_id_col;
	std::vector<faust_real> current_signed_values;
	faust_real current_value;
	int nbr_new_elt = 0;
	int i, nbr_elt_to_add;
	
	while (nb_elt_found < k)
	{
		current_id_row.resize(0);
		current_id_col.resize(0);
		
		current_value=M_abs.max(current_id_row,current_id_col);
		nbr_new_elt = current_id_row.size();
		
		nbr_elt_to_add=std::min(k-nb_elt_found,nbr_new_elt);
		
		current_id_row.resize(nbr_elt_to_add);
		current_id_col.resize(nbr_elt_to_add);
		current_signed_values.resize(nbr_elt_to_add);
		M.getCoeffs(current_signed_values,current_id_row,current_id_col);
		
		new_M.setCoeffs(current_signed_values,current_id_row,current_id_col);
		nb_elt_found += nbr_elt_to_add;
		M_abs.setCoeffs(0,current_id_row,current_id_col);
	
	/*
	std::cout<<std::endl<<std::endl;	
	std::cout<<"*** ITER *********************** : "<< nb_elt_found <<std::endl;	
	std::cout<<"*** M ****"<<std::endl;
	M.Display();
	
	
	std::cout<<"*** M_abs ****"<<std::endl;
	M_abs.Display();
	
	
	std::cout<<"*** new_M ****"<<std::endl;
	new_M.Display();
	*/
		
		
	}
	
	new_M.scalarMultiply(1/new_M.norm());
	M = new_M; 	
}








void prox_spcol(faust_mat & M,int k)
{
	int dim1 = M.getNbRow();
	int dim2 = M.getNbCol();
	if (k>dim1)
	{
		std::cerr << "ERROR prox_spcol : k > nb_Row of M " << std::endl;
		exit( EXIT_FAILURE);	
	}
	faust_mat new_M(dim1,dim2);
	faust_mat M_abs(dim1,dim2);
	faust_real current_max;
	M_abs = M;
	M_abs.abs();
	new_M.setZeros();
	
	faust_mat current_col(dim1,1);
	std::vector<int> id_row,id_col_mat,id_col,id_row_max,id_col_max;
	std::vector<faust_real> signed_values,values_per_Col;
	id_row.resize(dim1);
	id_col.assign(dim1,0);
	values_per_Col.resize(dim1);
	
	
	int nb_elt_found,nbr_new_elt,nbr_elt_to_add;
	
	for (int i=0;i<dim1;i++)id_row[i]=i;
	
	for (int j=0;j<dim2;j++)
	{	
		nb_elt_found = 0;
		id_col_mat.assign(dim1,j);
		M_abs.getCoeffs(values_per_Col,id_row,id_col_mat);
		current_col.setCoeffs(values_per_Col,id_row,id_col);//copie des coefficents de la matrice dans une matrice column
		std::cout<<"avant boucle while"<<std::endl;
		while (nb_elt_found < k)
		{
			id_row_max.resize(0);
			id_col_max.resize(0);
			
			// calculus of the max of the current column and its index
			current_max=current_col.max(id_row_max,id_col_max);
			nbr_new_elt = id_row_max.size();
		
			nbr_elt_to_add=std::min(k-nb_elt_found,nbr_new_elt);
			id_col_mat.assign(nbr_elt_to_add,j);
			id_row_max.resize(nbr_elt_to_add);
			signed_values.resize(nbr_elt_to_add);
			
			M.getCoeffs(signed_values,id_row_max,id_col_mat);
			new_M.setCoeffs(signed_values,id_row_max,id_col_mat);
			id_col_max.resize(nbr_elt_to_add);
			current_col.setCoeffs(0,id_row_max,id_col_max);
		
			nb_elt_found += nbr_elt_to_add;
		}
		
	}
	M=new_M;
	
}



void prox_splin(faust_mat & M,int k){M.init_from_file("facts0.txt");}
void prox_normcol(faust_mat & M, faust_real k){M.init_from_file("facts1.txt");}
















