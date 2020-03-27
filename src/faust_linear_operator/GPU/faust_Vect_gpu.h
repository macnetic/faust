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
#ifndef __FAUST_CU_VEC_H__
#define __FAUST_CU_VEC_H__

#include "faust_constant_gpu.h"
#include <vector>
#include <iterator>
#include "cuda.h"
#include "cuda_runtime.h"
#include "faust_BlasHandle_gpu.h"



#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif
#include "faust_exception.h"


#ifdef __COMPILE_SPMAT__
    #include "faust_SpBlasHandle_gpu.h"
    template <typename FPP,FDevice DEVICE> class MatSparse;
#endif

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    template <typename FPP,FDevice DEVICE> class Vect;
    template <typename FPP,FDevice DEVICE> class MatDense;
#ifdef __COMPILE_SPMAT__
    template <typename FPP,FDevice DEVICE> class MatSparse;
#endif

    template <typename FPP> class Vect<FPP,Gpu>
    {
    public :
        Vect() : dim(0), device(FAUST_DEFAULT_CUDA_DEVICE), data(NULL) {}
        Vect(const FPP  *data_, const faust_unsigned_int dim_, bool dataFromGPU, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
        Vect(const Vect<FPP,Gpu>& cu_v, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
        Vect(const Vect<FPP,Cpu>& v, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
        Vect(const faust_unsigned_int dim_, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );

        void resize(const faust_unsigned_int dim_, int device_);
        void resize(const faust_unsigned_int dim_){resize(dim_, device);}

        void copyFromHost(const FPP *data_, const faust_unsigned_int dim_, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
        void copyFromDevice(const FPP *data_, const faust_unsigned_int dim_, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
        void copyToHost(FPP *data_, const faust_unsigned_int dim_, cudaStream_t stream=0)const;
        void copyToDevice(FPP *data_, const faust_unsigned_int dim_, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0)const;
        void init(const Vect<FPP,Gpu>& cu_v, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
        void moveToDevice(int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);

        void operator=(const Vect<FPP,Gpu>& cu_v);
        void operator=(const Vect<FPP,Cpu>& v);

        FPP max() const;
        FPP min() const;
        FPP sum()const;
        FPP norm() const;
        FPP mean()const{return sum()/dim;}
        FPP dot(const Vect<FPP,Gpu>& cu_v, Faust::BlasHandle<Gpu> blasHandle)const;

        bool operator==(const Vect<FPP,Gpu>& cu_v)const;
        bool operator!=(const Vect<FPP,Gpu>& cu_v)const{return !((*this)==cu_v);}

        void operator+=(const Vect<FPP,Gpu>& cu_v);
        void operator-=(const Vect<FPP,Gpu>& cu_v);
        void operator*=(const Vect<FPP,Gpu>& cu_v);
        void operator/=(const Vect<FPP,Gpu>& cu_v);

        void operator+=(const FPP alpha);
        void operator-=(const FPP alpha);
        void operator*=(const FPP alpha);
        void operator/=(const FPP alpha);
        void scalarMultiply(const FPP alpha){this->operator*=(alpha);}
        void normalize(){scalarMultiply(1/norm());}

        void setValues(const FPP);
        void setValues(const FPP, const faust_unsigned_int);

        void setZeros();
        void setZeros(const faust_unsigned_int);

        void setOnes();
        void setOnes(const faust_unsigned_int);


        FPP mean_relative_error(const Vect<FPP,Gpu>& v)const;
        // multiply (*this) =  A * (*this)
        inline void  multiplyLeft(const Faust::MatDense<FPP,Gpu>& cu_A, Faust::BlasHandle<Gpu> blasHandle);
        #ifdef __COMPILE_SPMAT__
        void  multiplyLeft(const Faust::MatSparse<FPP,Gpu>& cu_S,Faust::SpBlasHandle<Gpu> spblasHandle);
        #endif

        ~Vect(){resize(0);}

        faust_unsigned_int size() const {return dim;}

        bool equality(Vect<FPP,Gpu> const &x, FPP precision=(FPP)0.001) const;

        int getDevice()const{return device;}

        FPP* getData(){return data;}
        const FPP* getData()const {return data;}

        void Display()const;
        void print_file(const char* filename)const;


    private:
        void _create(const faust_unsigned_int dim_, int device_);
        void _clear();


    private:
        faust_unsigned_int dim;
        static const char * m_className;
        FPP* data;
        int device;

        #ifdef __COMPILE_TIMERS__
        public:
            Faust::Timer t_local_multiplyLeft;
        #endif

    };
}


template<typename FPP>
inline void Faust::Vect<FPP,Gpu>::multiplyLeft(const Faust::MatDense<FPP,Gpu>& cu_A, Faust::BlasHandle<Gpu> blasHandle)
{
    Faust::Vect<FPP,Gpu> cu_y_copy;
    gemv(cu_A, *this, cu_y_copy, 1.0, 0.0, 'N', blasHandle);
    this->operator=(cu_y_copy);
}

#include "faust_Vect_gpu.hpp"

#endif
