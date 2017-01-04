/****************************************************************************/
/*                              Description:                                */
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
#ifndef __FAUST_Transform_H__
#define __FAUST_Transform_H__

#include <vector>
#include "faust_LinearOperator.h"


#include "faust_transform_algebra.h"

//#include "faust_Vect.h"
//#include "faust_MatDense.h"
//#include "faust_MatSparse.h"

//! \class Faust::Transform
//! \brief contains the data of the FAUST method (sparses matrix, lambda, number of non-zeros...) */


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

	template<typename FPP,Device DEVICE> class LinearOperator;

	template<typename FPP,Device DEVICE> class Transform;
	template<typename FPP,Device DEVICE> class Vect;
	template<typename FPP,Device DEVICE> class MatDense;
	template<typename FPP,Device DEVICE> class MatSparse;
	template<Device DEVICE> class BlasHandle;
	template<Device DEVICE> class SpBlasHandle;


	// forward definition of friend function
	template<typename FPP>
	Faust::Vect<FPP,Cpu> operator*(const Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu>& v);
	template<typename FPP>
	Faust::MatDense<FPP,Cpu> operator*(const Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu>& M);


	template<typename FPP>
	class Transform<FPP,Cpu> : public Faust::LinearOperator<FPP,Cpu>
	{

		public:
		void faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB)const;

        /** \brief Constructor
        * \param data : Vector including sparse matrix
        * \param totalNonZeros : Number of nonzeros Value in the data (all factorized matrix) */
		Transform();

        /** \brief Constructor
        * \param data : Vector including sparse matrix
        * \param lambda : the multiplicative scalar*/
        Transform(const std::vector<Faust::MatSparse<FPP,Cpu> >& facts, const FPP lambda_ = (FPP)1.0);

        Transform(const Transform<FPP,Cpu> & A);
	
		/** \brief 
        * check the factors validity of the faust, if the list of factors represents a valid matrix
        * */		
	void check_factors_validity() const;

        /** \brief Constructor
        * \param facts : Vector including dense matrix*/
        Transform(const std::vector<Faust::MatDense<FPP,Cpu> >&facts);

        //void get_facts(std::vector<Faust::MatSparse<FPP,Cpu> >& sparse_facts)const{sparse_facts = data;}
        //void get_facts(std::vector<Faust::MatDense<FPP,Cpu> >& facts)const;
        
	faust_unsigned_int size()const{return data.size();}
        void size(int size_)const{ data.resize(size_);}


        //template<Device DEVICE> class BlasHandle;
        //template<Device DEVICE> class SpBlasHandle;
        /** \brief Perform the product of all factorized matrix. */
        Faust::MatDense<FPP,Cpu> get_product(const char opThis='N')const;
        Faust::MatDense<FPP,Cpu> get_product(Faust::BlasHandle<Cpu> blas_handle,Faust::SpBlasHandle<Cpu> spblas_handle)const;
        // modif AL AL
        // Faust::MatDense<FPP,Cpu> get_product(Faust::BlasHandle<Cpu> blas_handle,Faust::SpBlasHandle<Cpu> spblas_handle)const
        // {return (*this).get_product();}
        Faust::MatSparse<FPP,Cpu> get_fact(faust_unsigned_int id) const;
	
        faust_unsigned_int getNbRow() const;
        faust_unsigned_int getNbCol() const;
        void print_file(const char* filename) const;
        void init_from_file(const char* filename);
        long long int get_total_nnz()const{return totalNonZeros;}
        void clear(){data.resize(0);totalNonZeros=0;}
	void push_back(const Faust::MatGeneric<FPP,Cpu>* M);
        void push_first(const Faust::MatGeneric<FPP,Cpu>* M);
        void pop_back(Faust::MatGeneric<FPP,Cpu>* M);
        void pop_first(Faust::MatGeneric<FPP,Cpu>* M);
        void pop_first(Faust::MatGeneric<FPP,Cpu>* M) const;
        void pop_first(Faust::MatSparse<FPP,Cpu>& S) const;
        void Display()const;
        void transpose();
        void updateNonZeros();
		void setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const;
        ///(*this) = (*this) * A
        void multiply(const Transform<FPP,Cpu> & A);
        ///(*this) = A * (*this)
        void multiplyLeft(const Transform<FPP,Cpu> & A);
        void scalarMultiply(const FPP scalar);
        FPP spectralNorm(const int nbr_iter_max, FPP threshold, int &flag) const;
        ~Transform(){}
	

	/*!
	* \brief multiplication between with vector x
	*  x = op((*this)) * x
	<br>
    	* op((*this)) = (*this) if opThis='N', op((*this)= = transpose((*this)) if opThis='T'<br>
	*! \tparam  x :  the vector to be multiplied
	*! \tparam opThis : character 
        */                               					
	Faust::Vect<FPP,Cpu> multiply(const Faust::Vect<FPP,Cpu> x,const char opThis='N') const;


	
	Faust::MatDense<FPP,Cpu> multiply(const Faust::MatDense<FPP,Cpu> A,const char opThis='N') const;

       
		void operator=(const Transform<FPP,Cpu>&  f){data=f.data;totalNonZeros=f.totalNonZeros;}
        /// add all of the sparse matrices from f.data to this->data
        void operator*=(const FPP  scalar){scalarMultiply(scalar);};
        void operator*=(const Transform<FPP,Cpu>&  f){multiply(f);};
        /// add the sparse matrix S to this->data
        void operator*=(const Faust::MatSparse<FPP,Cpu>&  S){push_back(S);totalNonZeros+=S.getNonZeros();}
	
	#ifdef __COMPILE_TIMERS__
		void print_timers() const;
	#endif




	private:
	long long int totalNonZeros;
	static const char * m_className;
	std::vector<Faust::MatGeneric<FPP,Cpu>*> data;
	
	#ifdef __COMPILE_TIMERS__
		mutable std::vector<Faust::Timer> t_multiply_vector;
	#endif 


	// friend function
	friend Faust::Vect<FPP,Cpu> Faust::operator*<>(const Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu>& v);
	friend Faust::MatDense<FPP,Cpu> Faust::operator*<>(const Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu>& M);
    };

}



#include "faust_Transform.hpp"


#endif
