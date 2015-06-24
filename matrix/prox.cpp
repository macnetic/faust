#include "prox.h"
#include <vector>
#include "faust_constant.h"
#include <iostream>

void prox_sp(faust_mat & M,int k)
{	
	int dim1 = M.getNbRow();
	int dim2 = M.getNbCol();
	if (k<dim1*dim2)
	{	
		faust_mat new_M(dim1,dim2);
		faust_mat M_abs(dim1,dim2);
		M_abs = M;
		M_abs.abs();
		new_M.setZeros();
	
	
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
		M = new_M;
	}
	M.scalarMultiply(1/M.norm());
	

	
}








void prox_spcol(faust_mat & M,int k)
{
	int dim1 = M.getNbRow();
	int dim2 = M.getNbCol();
	if (k<dim1)
	{
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
		
		//memcpy(&((M.getData())[j*dim1]),&(current_col),dim1*sizeof(faust_real));
		//std::cout<<"avant boucle while"<<std::endl;
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
	faust_real normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}
}



void prox_splin(faust_mat & M,int k)
{
	M.transpose();
	prox_spcol(M,k);
	M.transpose();

}

//  normcol of the zero matrix equal to the matrix with all elements equal to s,
//  in order to have the same behaviour as matlab prox
void prox_normcol(faust_mat & M,faust_real s)
{

	int dim1 = M.getNbRow();
	int dim2 = M.getNbCol();
	if (s<0)
	{
		std::cerr << "ERROR prox_normcol : s < 0" << std::endl;
		exit( EXIT_FAILURE);	
	}
	
	
	faust_mat current_col(dim1,1);
	std::vector<int> id_row,id_col_mat,id_col;
	std::vector<faust_real> values_per_Col;
	id_row.resize(dim1);
	id_col.assign(dim1,0);
	values_per_Col.resize(dim1);
	faust_real norm_col;
	
	if (s == 0)
	{
		M.setZeros();
	}else
	{
		
		for (int i=0;i<dim1;i++)id_row[i]=i;
	
		for (int j=0;j<dim2;j++)
		{	
			id_col_mat.assign(dim1,j);
			M.getCoeffs(values_per_Col,id_row,id_col_mat);
			current_col.setCoeffs(values_per_Col,id_row,id_col);//copie des coefficents de la matrice dans une matrice column
			norm_col = current_col.norm();
			
			if (norm_col == 0)
			{
				M.setCoeffs(s,id_row,id_col_mat);
			}else
			{
				for (int k=0;k<values_per_Col.size();k++)
				{
					values_per_Col[k]=values_per_Col[k]/norm_col*s;
				}
				M.setCoeffs(values_per_Col,id_row,id_col_mat);
			}
				
		}
	
	}
	
}


void prox_normlin(faust_mat & M,faust_real s)
{
	M.transpose();
	prox_normcol(M,s);
	M.transpose();
	
}



void prox_supp(faust_mat & M, const faust_mat & supp)
{
	if ( (supp.getNbRow() != M.getNbRow()) || (supp.getNbCol() != M.getNbCol()) )
	{
		std::cerr << "ERROR prox_supp : dimensions of the matrix mismatch " << std::endl;
		exit( EXIT_FAILURE);	
	}
	faust_real value;
	for (int i=0;i<M.getNbRow();i++)
	{
		for (int j=0;j<M.getNbCol();j++)
		{
			value=supp.getCoeff(i,j);
			if (value == 0)
			{
				M.setCoeff(0,i,j);
			}
		}
	}
	
	faust_real normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}
	
}












void prox_sp_pos(faust_mat & M,int k)
{
	//treshold de la matrice
	for (int i=0;i<(M.getNbRow() * M.getNbCol());i++)
	{
		if ((((M.getData()))[i]) < 0)
		{
			(((M.getData()))[i])=0;
		}
		
	}
	
	prox_sp(M,k);
	
}



// M needs to be square and k must divide dimension of M
void prox_blkdiag(faust_mat & M,int k)
{   
	int i,j;
	int Msize = M.getNbRow();
	if (Msize != M.getNbCol())
	{
		std::cerr << "ERROR prox_blkdiag : input matrix must be square " << std::endl;
		exit( EXIT_FAILURE);	
	}
	
	int sizeblock = Msize/k;
	faust_mat new_M(Msize,Msize);
	new_M.setZeros();
	std::cout<<"sizeblock : "<<sizeblock<<std::endl;
	
	for (int ii=0;ii<k;ii++)
	{
		for (i=sizeblock*ii;i<sizeblock*(ii+1);i++)
		{
			for (j=sizeblock*ii;j<sizeblock*(ii+1);j++)
			{	
				//std::cout<<M.getCoeff(i,j)<<std::endl;
				new_M.setCoeff(M.getCoeff(i,j),i,j);	
			}
		
		}
		
	}

	M = new_M;

	faust_real normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}	
		
	
}



void prox_toeplitz(faust_mat & M, int k)
{	
	int nl = M.getNbRow();
	int nc = M.getNbCol();
	faust_mat crit(nl+nc-1,1);
	faust_mat new_M(nl,nc);
	new_M.setZeros();
	std::vector<faust_real> diag_meanS;
	std::vector<int> num_elt_diagS;
	int num_elt_diag;
	int min_dim = std::min(nl,nc);
	int shift_dim = nc-nl;
	faust_real current_mean;
	int id_col,id_row;
	
	diag_meanS.resize(nl+nc-1);
	num_elt_diagS.resize(nl+nc+1);
	
	if (k>nl+nc-1)
	{
		std::cerr << "ERROR prox_toeplitz : k > (nb_Row+nbCol-1) of M " << std::endl;
		exit( EXIT_FAILURE);	
	}
	
	std::cout<<"boucle critere :"<<std::endl;
	for (int i=-nl+1;i<nc;i++)
	{	
		/*
		if ( ((i>=shift_dim) && (i<=0)) || ((i<=(shift_dim)) && (i>=0)) )
		{	std::cout<<"min_dim"<<std::endl;
			num_elt_diag = min_dim; 
		}else
		{	
			std::cout<<"soustraction"<<std::endl;
			num_elt_diag = min_dim - (abs(shift_dim) - abs(i));
		}
		*/
		if (shift_dim < 0)
		{
			if (i<=0)
			{
				if (i<shift_dim)
				{	std::cout<<"cas1"<<std::endl;
					num_elt_diag = min_dim-abs(i+abs(shift_dim));
				}else
				{	
					std::cout<<"cas2"<<std::endl;
					num_elt_diag = min_dim;
				}
			}else
			{
				std::cout<<"cas3"<<std::endl;
				num_elt_diag = min_dim-i;
			}
		}else
		{
			if (i>=0)
			{
				if (i>shift_dim)
				{	std::cout<<"cas4"<<std::endl;
					num_elt_diag = min_dim-abs(i-abs(shift_dim));
				}else
				{	
					std::cout<<"cas5"<<std::endl;
					num_elt_diag = min_dim;
				}
			}else
			{
				std::cout<<"cas6"<<std::endl;
				num_elt_diag = min_dim-abs(i);
			}
		}
		
		
		num_elt_diagS[i+nl-1]=num_elt_diag;
		
		if (i < 0)
		{
			id_row = -i;
			id_col = 0;
		}else
		{
			id_row = 0;
			id_col = i;
		}
		std::cout<<"id_row : "<<id_row<<std::endl;
		std::cout<<"id_col : "<<id_col<<std::endl;
		std::cout<<"nombre diagonal : "<<num_elt_diag<<std::endl;
		
		current_mean = 0;
		std::cout<<"coeff diag"<<std::endl;
		for (int l=0;l<num_elt_diag;l++)
		{	
			std::cout<<M.getCoeff(id_row+l,id_col+l)<<std::endl;
			current_mean = current_mean + M.getCoeff(id_row+l,id_col+l);
		}
		
		current_mean = current_mean/((faust_real) num_elt_diag);
		std::cout<<"mean : "<<current_mean<<std::endl;
		
		diag_meanS[i+nl-1]=current_mean;
		crit.setCoeff(current_mean*current_mean*((faust_real)num_elt_diag),i+nl-1,0);
		
	}
	
	std::cout<<" crit : "<<std::endl;
	crit.Display();
	prox_sp(crit,k);
	int ll;
	std::cout<<"boucle sparse"<<std::endl;
	for (int i=0;i<nl+nc-1;i++)
	{
		if (crit.getCoeff(i,0) != 0)
		{
			ll = i-nl+1;
			if (ll < 0)
			{
				id_row = -ll;
				id_col = 0;
			}else
			{
				id_row = 0;
				id_col = ll;
			}
			current_mean = diag_meanS[i];
			for (int l=0;l<num_elt_diagS[i];l++)
			{
				new_M.setCoeff(current_mean,id_row+l,id_col+l);
			}
		}
		
	}
	
	M = new_M;

	faust_real normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}	
	
	
	
	
}
	
	
	
	
	
	


















