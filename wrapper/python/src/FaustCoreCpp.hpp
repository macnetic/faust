/****************************************************************************/
/*                              Description:                                */
/*  C++ Hpp file wrapping the  Class Faust::Transform,                      */
/*  the methods are quite the same but this works with pointer              */
/*  which are easier to handle from Python                                  */
/*                                                                          */     
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
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

#include "faust_Transform.h"
#include <iostream>
#include <exception>

template<typename FPP>
void FaustCoreCpp<FPP>::push_back(FPP* valueMat, unsigned int nbrow,unsigned int nbcol)
{
	Faust::MatDense<FPP,Cpu> dense_mat(valueMat,nbrow,nbcol);
	Faust::MatSparse<FPP,Cpu> sparse_mat(dense_mat);
	//sparse_mat.Display();
	this->transform.push_back(&sparse_mat);
	
	
	
	
}

template<typename FPP>
void FaustCoreCpp<FPP>::multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x,bool isTranspose)const
{
	
	
	char op = 'N';
	if (isTranspose)
		op = 'T';

	unsigned int nbRowThis,nbColThis;
	
	
	this->setOp(isTranspose,nbRowThis,nbColThis);


	if ( (nbrow_y != nbRowThis) | (nbrow_x != nbColThis) | (nbcol_y != nbcol_x) )
	{	
		std::cout<<"nbRowThis "<<nbRowThis<<" must be equal to nb row y  "<<nbrow_y<<std::endl;
		std::cout<<"nbColThis "<<nbColThis<<" must be equal to nb row x  "<<nbrow_x<<std::endl;
		std::cout<<"nbcol_y "<<nbcol_y<<" must be equal to nbcol_x  "<<nbcol_x<<std::endl;
		handleError("FaustCpp"," multiply : invalid dimension");
	}
	if (nbcol_x == 1)
	{
		Faust::Vect<FPP,Cpu> X(nbrow_x,value_x);
		Faust::Vect<FPP,Cpu> Y;

		
		Y = this->transform.multiply(X,op);

		memcpy(value_y,Y.getData(),sizeof(FPP)*nbrow_y);
	}else
	{
		Faust::MatDense<FPP,Cpu> X(value_x,nbrow_x,nbcol_x);
		Faust::MatDense<FPP,Cpu> Y;

		Y = this->transform.multiply(X,op);

		memcpy(value_y,Y.getData(),sizeof(FPP)*nbrow_y*nbcol_y);
	}
	
	
}


template<typename FPP>
void FaustCoreCpp<FPP>::setOp(const bool isTransposed,unsigned int& nbRowOp, unsigned int& nbColOp)const
{
	char trans_flag('N');
	if (isTransposed)
		trans_flag='T';

	faust_unsigned_int nb_row,nb_col;
	this->transform.setOp(trans_flag,nb_row,nb_col);
	nbRowOp=(unsigned int) nb_row;
	nbColOp=(unsigned int) nb_col;
	




}



template<typename FPP>
unsigned long long FaustCoreCpp<FPP>::nnz() const
{
    return this->transform.get_total_nnz();
}


template<typename FPP>
double FaustCoreCpp<FPP>::norm() const
{
    double precision = 0.001;
    faust_unsigned_int nbr_iter_max = 100;
    int flag; // not used yet
    return this->transform.spectralNorm(nbr_iter_max, precision, flag);
}

template<typename FPP>
double FaustCoreCpp<FPP>::get_nb_factors() const
{
    double nb_fact = this->transform.size();
    return nb_fact;
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::get_fact_nb_rows(unsigned int& i) const
{
    Faust::MatGeneric<FPP,Cpu>* const factor_generic = this->transform.get_fact(i);
    unsigned int nb_rows = factor_generic->getNbRow();
    delete factor_generic;
    return nb_rows;
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::get_fact_nb_cols(unsigned int& i) const
{
    Faust::MatGeneric<FPP,Cpu>* const factor_generic = this->transform.get_fact(i);
    unsigned int nb_cols = factor_generic->getNbCol();
    delete factor_generic;
    return nb_cols;
}

template<typename FPP>
void FaustCoreCpp<FPP>::get_fact(unsigned int& i, FPP* fact_ptr) const
{
    Faust::MatGeneric<FPP,Cpu>* const factor_generic = this->transform.get_fact(i);
    Faust::MatDense<FPP,Cpu> dense_factor;

    switch (factor_generic->getType())
    {
        case Dense :
            {
                Faust::MatDense<FPP,Cpu>* factor_dense_ptr = dynamic_cast<Faust::MatDense<FPP,Cpu>* > (factor_generic);
                dense_factor = (*factor_dense_ptr);
            }
            break;

        case Sparse :
            {
                Faust::MatSparse<FPP,Cpu>* factor_sparse_ptr = dynamic_cast<Faust::MatSparse<FPP,Cpu>* > (factor_generic);
                dense_factor = (*factor_sparse_ptr);
            }
            break;

        default:
            throw std::runtime_error("get_fact : unknown type of the factor matrix.");
    }


    memcpy(fact_ptr, dense_factor.getData(),
            sizeof(FPP)*factor_generic->getNbCol()*factor_generic->getNbRow());

    delete factor_generic;

}

template<typename FPP>
void FaustCoreCpp<FPP>::save_mat_file(const char* filepath, bool transpose_flag) const
{
//    std::cout << "FaustCoreCpp::save_mat_file()" << std::endl;
    this->transform.save_mat_file(filepath, transpose_flag);
}
