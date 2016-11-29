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

template<typename FPP>
void FaustCpp<FPP>::push_back(FPP* valueMat, unsigned int nbrow,unsigned int nbcol)
{
	Faust::MatDense<FPP,Cpu> dense_mat(valueMat,nbrow,nbcol);
	Faust::MatSparse<FPP,Cpu> sparse_mat(dense_mat);
	//sparse_mat.Display();
	this->transform.push_back(sparse_mat);
	
	
	
	
}

template<typename FPP>
void FaustCpp<FPP>::multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x,bool isTranspose)const
{
	if ( (nbrow_y != this->transform.getNbRow()) | (nbrow_x != this->transform.getNbCol()) | (nbcol_y != nbcol_x) )
		handleError("FaustCpp"," invalid dimension");

	char op = 'N';
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
void FaustCpp<FPP>::setOp(const bool isTransposed,unsigned int& nbRowOp, unsigned int& nbColOp)const
{
	char trans_flag('N');
	if (isTransposed)
		trans_flag='T';
	faust_unsigned_int nb_row,nb_col;
	this->transform.setOp(trans_flag,nb_row,nb_col);
	nbRowOp=(unsigned int) nb_row;
	nbColOp=(unsigned int) nb_col;
	




}


