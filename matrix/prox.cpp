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

	std::vector<int> id_row;
	std::vector<int> id_col;
	std::vector<int> values;

	id_row.assign(k,-1);
	id_col.assign(k,-1);
	values.assign(k,-1);

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