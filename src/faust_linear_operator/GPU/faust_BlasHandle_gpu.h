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
#ifndef BLASHANDLE_GPU_H
#define BLASHANDLE_GPU_H
#include "faust_constant_gpu.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"


namespace Faust
{

    template <FDevice DEVICE> class BlasHandle;


    template <>
    class BlasHandle<Gpu>
    {

    public :
        BlasHandle(cublasHandle_t cublasHandle):m_cublasHandle(cublasHandle){}
        BlasHandle(){}
        BlasHandle(const BlasHandle<Gpu> & blashandleGPU): m_cublasHandle(blashandleGPU.m_cublasHandle){}
        void operator=(BlasHandle<Gpu> const& blashandleGPU){m_cublasHandle=blashandleGPU.m_cublasHandle;}
        void Init(){cublasStatus_t cublasStat=cublasCreate(&m_cublasHandle);}
        void Destroy(){cublasDestroy(m_cublasHandle);}
        const cublasHandle_t & GetHandle() const{return m_cublasHandle;}

        private :
        cublasHandle_t m_cublasHandle;


    };

}

#endif
