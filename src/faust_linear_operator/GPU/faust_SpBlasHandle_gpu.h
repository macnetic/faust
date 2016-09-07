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
#ifndef SPBLASHANDLE_GPU_H
#define SPBLASHANDLE_GPU_H

#include "faust_constant_gpu.h"

#include "cuda.h"
#include "cusparse.h"


namespace Faust
{

    template <Device DEVICE> class SpBlasHandle;

    template <>
    class SpBlasHandle<Gpu>
    {

    public :

        SpBlasHandle(){}
        SpBlasHandle(cusparseHandle_t cusparseHandle):m_cusparseHandle(cusparseHandle){}
        SpBlasHandle(const SpBlasHandle<Gpu> & spblashandleGPU): m_cusparseHandle(spblashandleGPU.m_cusparseHandle){}
        void operator=(SpBlasHandle<Gpu> const& spblashandleGPU){}//{m_cusparseHandle=spblashandleGPU.m_cusparseHandle;}
        void Init(){cusparseStatus_t cusparseStat=cusparseCreate(&m_cusparseHandle);}
        void Destroy(){cusparseDestroy(m_cusparseHandle);}
        const cusparseHandle_t & GetHandle() const {return m_cusparseHandle;}
        private :
        cusparseHandle_t m_cusparseHandle;


    };

}

#endif
