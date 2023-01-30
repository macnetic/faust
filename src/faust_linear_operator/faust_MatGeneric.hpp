/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
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
#include "faust_exception.h"

// useful for optimize with multiplication
#include "faust_Timer.h"
#include "faust_MatBSR.h"
#include "faust_MatPerm.h"
#include "faust_MatButterfly.h"


template<typename FPP,FDevice DEVICE>
void Faust::MatGeneric<FPP,DEVICE>::setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const
{
    if(op == 'N')
    {
        nbRowOp=getNbRow();
        nbColOp=getNbCol();
    }
    else if(op == 'T' || op == 'H')
    {
        nbRowOp=getNbCol();
        nbColOp=getNbRow();
    }
    else
        handleError("Faust::MatGeneric::","setOp : invalid character");
}



//template <typename FPP, FDevice DEVICE>
template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::optimize(Faust::MatDense<FPP,DEV> const & M,Faust::MatSparse<FPP,DEV> const & S)
{
	//std::cout<<"DEBUT OPTIMIZE "<<std::endl;


	if ( (M.getNbCol() != S.getNbCol()) | (M.getNbRow() != S.getNbRow()) )
		handleError("Faust::MatGeneric::", " Faust::optimize : matrix M and S have not the same size");

	Faust::Vect<FPP,DEV> x_dense(M.getNbCol());

	for (int i=0;i<M.getNbCol();i++)
	{
//		x_dense[i]=i*0.005;
		x_dense.set_coeff(i, i*.005);
	}

	Faust::Vect<FPP,DEV> const x(x_dense);
	Faust::Vect<FPP,DEV> x_sparse(x_dense);

	int nb_mult=10;
	Faust::Timer t_dense,t_sparse;
	for (int i=0;i<nb_mult;i++)
	{
		x_sparse=x;
		x_dense=x;
		t_sparse.start();
		S.multiply(x_sparse,'N');
		t_sparse.stop();
		t_dense.start();
		M.multiply(x_dense,'N');
		t_dense.stop();

	}
	//float density = ((float)S.getNonZeros())/((float)(S.getNbRow()*S.getNbCol()));
	float density = S.density();
	//std::cout<<" density "<<density<<std::endl;
	//std::cout<<" tps sparse "<<t_sparse.get_time()<<std::endl;
	//std::cout<<" tps dense "<<t_dense.get_time()<<std::endl;

	//if (M.getNbCol() != M.getNbRow())
	if (t_sparse.get_time() <= t_dense.get_time())
	{
		//std::cout<<" CHOICE SPARSE "<<t_dense.get_time()<<std::endl;
		return new MatSparse<FPP,DEV>(S);
	}else
	{
		//std::cout<<" CHOICE DENSE "<<t_dense.get_time()<<std::endl;
		return new MatDense<FPP,DEV>(M);
	}

	//std::cout<<"FIN OPTIMIZE "<<t_sparse.get_time()<<std::endl;


}

template<typename FPP,FDevice DEVICE>
Faust::MatGeneric<FPP,DEVICE>::~MatGeneric()
{

}

template<typename FPP,FDevice DEVICE>
void Faust::MatGeneric<FPP,DEVICE>::Display() const
{
	std::cout << to_string() << std::endl;
}

template<typename FPP,FDevice DEVICE>
std::string Faust::MatGeneric<FPP,DEVICE>::to_string(MatType type, const bool transpose /* set to false by default */, const bool displaying_small_mat_elts) const
{
	return this->to_string(getNbRow(), getNbCol(), transpose, this->density(), this->getNonZeros(), this->is_identity, type);
}

template<typename FPP,FDevice DEVICE>
std::string Faust::MatGeneric<FPP,DEVICE>::to_string(const bool transpose /* set to false by default */, const bool displaying_small_mat_elts) const
{
	MatType type = Dense;
	if(dynamic_cast<const MatSparse<FPP,DEVICE>*>(this))
		type = Sparse;
	else if (dynamic_cast<const MatDense<FPP,DEVICE>*>(this))
		type = Dense;
	else if (dynamic_cast<const MatDiag<FPP>*>(this))
		type = Diag;
	else if (dynamic_cast<const MatBSR<FPP, DEVICE>*>(this))
		type = BSR;
	else if (dynamic_cast<const MatPerm<FPP, DEVICE>*>(this))
		type = Perm;
	else if (dynamic_cast<const MatButterfly<FPP, DEVICE>*>(this))
		type = Butterfly;
	else
		throw std::runtime_error("Unhandled matrix type in MatGeneric::to_string()"); // shouldn't happen
	return this->to_string(getNbRow(), getNbCol(), transpose, this->density(), this->getNonZeros(), this->is_identity, type);
}

template<typename FPP, FDevice DEVICE>
std::string Faust::MatGeneric<FPP,DEVICE>::get_scalar_type_str()
{
	std::string type_str;
	if (! std::is_same<FPP,Real<FPP>>::value)
		type_str = "complex";
	else
	{
		if(std::is_same<FPP, float>::value)
			type_str = "float";
		else
			type_str = "double";
	}
	return type_str;
}


template<typename FPP,FDevice DEVICE>
std::string Faust::MatGeneric<FPP,DEVICE>::to_string(int32_t nrows, int32_t ncols, bool transpose, Real<FPP> density, int32_t nnz, bool is_identity, MatType type) const
{
	std::ostringstream str;
	str << " (" << MatGeneric<FPP,DEVICE>::get_scalar_type_str() << ") ";
	if(type == Dense)
		str << "DENSE,";
	else if(type == Sparse)
		str << "SPARSE,";
	else if(type == Diag)
		str << "DIAG,";
	else if(type == BSR)
		str << "BSR,";
	else if(type == Perm)
		str << "PERM,";
	else if(type == Butterfly)
		str << "Butterfly,";
	else if(type == None)
		str << "UNKNOWN MATRIX TYPE,";
	else
		throw std::runtime_error("Invalid MatType passed to MatGeneric::to_string()");
	str << " size ";
	if(transpose)
		str << ncols << "x" << nrows;
	else
		str << nrows << "x" << ncols;
	if(type == BSR) //TODO: refactor with above almost same code into a new MatBSR member function
		str << dynamic_cast<const MatBSR<FPP, Cpu>*>(this)->to_string_blocks(transpose);
	str << ", density "<< density <<", nnz "<< nnz; 
	if(type == BSR)
		str <<  " (nnz blocks: " << dynamic_cast<const MatBSR<FPP, Cpu>*>(this)->getNBlocks() << ")";
//	str << ", addr: " << this;
	str << std::endl;
	if (is_identity)
		str <<" identity matrix flag" << std::endl;
	return str.str();
}

