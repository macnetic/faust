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
#ifndef FAUST_VECT_H
#define FAUST_VECT_H

#include <Eigen/Dense>
#include "faust_constant.h"
#include <vector>
#include <iterator>
#include "faust_exception.h"
#include <complex>
#include <iostream>


#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif


//! \class Faust::Vect faust_Vect.h
//! \brief Class template representing dense vector <br>
//! This class implements basic linear algebra operation (addition, multiplication, frobenius and spectral norm...) <br>





//template<typename FPP>
//


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

	template<typename FPP>
		FPP fabs(FPP c);

	template<typename FPP>
		FPP fabs(std::complex<FPP> c);

    template<typename FPP,FDevice DEVICE>
		class Vect;

    template<typename FPP,FDevice DEVICE>
		class MatDense;

    template<typename FPP,FDevice DEVICE>
		class MatSparse;


	template<typename FPP>
		class MatDiag;

	template<typename FPP,FDevice DEVICE>
		class MatBSR;

   // friend function of faust_linear_algebra.h
   template<typename FPP>
   void gemv(const MatDense<FPP,Cpu> & A,const Vect<FPP,Cpu> & x,Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);

	template<typename FPP,FDevice DEVICE> class Transform;
	template<typename FPP> class TransformHelperPoly;
	template<typename FPP, FDevice DEVICE> class TransformHelper;



    template<typename FPP>
    class Vect<FPP,Cpu>
    {
        template<class,FDevice> friend class Vect;
//		friend double Transform<FPP,Cpu>::normL1(const bool transpose) const;
//		friend Vect<FPP,Cpu> MatDiag<FPP>::multiply(const Vect<FPP,Cpu> & vec) const;
		friend class Transform<FPP,Cpu>;
		friend class MatDiag<FPP>;
		friend class MatDense<FPP,Cpu>;
		friend class MatSparse<FPP,Cpu>; // TODO: limit friendship to the necessary functions
//		friend void  MatDiag<FPP>::multiply(Vect<FPP,Cpu> & vec, char opThis) const;

		friend class MatBSR<FPP, Cpu>;
		friend TransformHelperPoly<FPP>; // TODO: limit to needed member functions only
		friend TransformHelper<FPP,Cpu>; // TODO: limit to needed member functions only (sliceMultiply at least)

        public :
        Vect() : dim(0), vec() {}
        Vect(const int _dim) : dim(_dim), vec(_dim){}

        // template<typename U>
        // Vect(const Vect<U>& v) : dim(v.dim), vec(v.vec){}
        template<typename FPP1>
        Vect(const Vect<FPP1,Cpu>& v); // copy-convert scalar ctor
		Vect(const Vect<FPP,Cpu>& v); //copy ctor
        Vect(const faust_unsigned_int dim_, const FPP* data_);
		Vect(Vect<FPP,Cpu>&& v); //move ctor

        FPP* getData(){return vec.data();}

	bool isReal() const;

        const FPP* getData() const {return vec.data();}
        void setOnes();
        void Display() const;
        void print_file(const char* filename)const;

        faust_unsigned_int size() const {return dim;}
        void resize(const int new_dim);
        double norm(){return vec.norm();}
        void scalarMultiply(FPP const scalar){vec *= scalar;}
	template<typename FPP1>
	void scalarMultiply(FPP1 const scalar){vec *= scalar;}
		void conjugate();
		FPP normL1() const;
		void normalize();


        // multiply (*this) =  A * (*this)
        // modif AL AL
        //void  multiplyLeft(MatDense<FPP,Cpu> const& A){gemv(A, *this, *this, 1.0, 0.0, 'N');}
        //!  \brief
        //! \brief Vect::multiplyLeft is used to replace this by A * (*this)
        //! \param A is a matrix (dense or sparse).
        void  multiplyLeft(MatDense<FPP,Cpu> const& A); //{gemv(A, *this, *this, 1.0, 0.0, 'N');}
        void  multiplyLeft(MatSparse<FPP,Cpu> const& S,const char TransS='N');

		bool isZero() const {return vec.isZero(0);};
        FPP sum()const{return vec.sum();}
        FPP mean()const{return vec.mean();}
        FPP dot(const Vect<FPP,Cpu> & v)const;
//		FPP min_coeff() const {return vec.minCoeff();};
		FPP min_coeff() const {int index; return this->min_coeff(&index);};
		FPP min_coeff(int *index) const { int col_index; return vec.minCoeff(index, &col_index); }
		FPP max_coeff(int *index) const;


        template<typename FPP1>
        void operator=(Vect<FPP1,Cpu> const& y); // copy assignment from different scalar type
        Vect<FPP,Cpu>& operator=(const Vect<FPP,Cpu>& y); // copy assignment
        Vect<FPP,Cpu>& operator=(Vect<FPP,Cpu>&& y); // move assignment

        void operator*=(const FPP alpha);
        void operator+=(const FPP alpha);
        void operator-=(const FPP alpha);

        void operator+=(const Vect<FPP,Cpu> & v);
        void operator-=(const Vect<FPP,Cpu> & v);
		void operator-=(const FPP* v_data);

        FPP mean_relative_error(const Vect<FPP,Cpu> & v);

        FPP& operator[](faust_unsigned_int i){return vec(i);}
        const FPP& operator[](faust_unsigned_int i)const{return vec(i);}
		void set_coeff(faust_unsigned_int i, const FPP& val) { (*this)[i] = val;};

        const FPP& operator()(faust_unsigned_int i)const{return vec(i);}
        bool equality(Vect<FPP,Cpu> const &x, const Real<FPP> precision) const;
		FPP mean();

		static Vect<FPP, Cpu>* rand(faust_unsigned_int size);
		void setRand();

        // friend algebra
        friend void gemv<>(const MatDense<FPP,Cpu> & A,const Vect<FPP,Cpu> & x,Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);


        private:
        faust_unsigned_int dim;
        static const char * m_className;
        Eigen::Matrix<FPP, Eigen::Dynamic,1> vec;

        #ifdef __COMPILE_TIMERS__
            public:
            Timer t_local_multiplyLeft;
        #endif

    };
}





#include "faust_Vect.hpp"

#endif
