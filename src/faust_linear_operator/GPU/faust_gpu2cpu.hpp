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
#ifndef __FAUST_CU2FAUST_HPP__
#define __FAUST_CU2FAUST_HPP__

#include "faust_Vect.h"
#include "faust_MatDense.h"
#include "faust_Vect_gpu.h"
#include "faust_MatDense_gpu.h"

#ifdef __COMPILE_SPMAT__
   #include "faust_MatSparse.h"
   #include "faust_MatSparse_gpu.h"
#endif
#include <typeinfo>




template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::Vect<FPP,Cpu>& v, const Faust::Vect<FPP1,Gpu>& cu_v, cudaStream_t stream/*=0*/ )
{
   if(cu_v.getData() != NULL)
   {
      v.resize(cu_v.size());
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(cu_v.getDevice());

      if(typeid(FPP)!=typeid(FPP1))
      {
         FPP1* data_tmp = new FPP1[cu_v.size()];
         faust_cudaMemcpyAsync(data_tmp, cu_v.getData(), cu_v.size()*sizeof(FPP1), cudaMemcpyDeviceToHost, stream);
	      for (int i=0 ; i<v.size() ;i++)
            v[i] = (FPP) data_tmp[i];
	      delete [] data_tmp ; data_tmp=NULL;
      }
      else
      {
         faust_cudaMemcpyAsync(v.getData(), cu_v.getData(), cu_v.size()*sizeof(FPP), cudaMemcpyDeviceToHost, stream);
      }



      faust_cudaSetDevice(currentGPU);
   }
   else if (cu_v.size()==0)
      v.resize(0);
   else
      handleError("faust_gpu2cpu","NULL pointer for non-empty vector");

}

template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::MatDense<FPP,Cpu>& M, const Faust::MatDense<FPP1,Gpu>& cu_M, cudaStream_t stream/*=0*/)
{
   const int dim1 = cu_M.getNbRow();
   const int dim2 = cu_M.getNbCol();

   M.resize(dim1, dim2);

   if(cu_M.estNulle())
      M.setZeros();
   else if(cu_M.estIdentite())
      M.setEyes();
   else if(cu_M.getData() != NULL)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(cu_M.getDevice());


      if(typeid(FPP)!=typeid(FPP1))
      {
   	   FPP1* data_tmp = new FPP1[dim1*dim2];
         faust_cudaMemcpyAsync(data_tmp, cu_M.getData(), dim1*dim2*sizeof(FPP1), cudaMemcpyDeviceToHost, stream);
	      for (int i=0 ; i<dim1*dim2 ;i++)
            M[i] = (FPP) data_tmp[i];
	      delete [] data_tmp ; data_tmp=NULL;
      }
      else
         faust_cudaMemcpyAsync(M.getData(), cu_M.getData(), dim1*dim2*sizeof(FPP), cudaMemcpyDeviceToHost, stream);

      faust_cudaSetDevice(currentGPU);
   }
   else if(dim1*dim2 != 0)
      handleError("faust_gpu2cpu","NULL pointer for non-empty matrix");
}

#ifdef __COMPILE_SPMAT__
/*!
	* \brief copy a faust_cu_spmat cu_S into  a Faust::MatSparse S
	*\warning give a corrupted Faust::MatSparse S (i.e S.nnz not equal to S.mat.nonZeros())
	*/

///WARNING GIVE CORRUPTED Faust::MatSparse (i.e S.nnz != S.mat.nonZeros()
template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::MatSparse<FPP,Cpu>& S, const Faust::MatSparse<FPP1,Gpu>& cu_S, cudaStream_t stream/*=0*/)
{
      const int nnz = cu_S.getNonZeros();
      const int dim1 = cu_S.getNbRow();
      const int dim2 = cu_S.getNbCol();


       if(dim1*dim2>0 && nnz==0)
	{
		S.resize(dim1,dim2);
         	S.setZeros();
	}else
	{
      		if(dim1*dim2>0 && nnz>0)
      		{
         		if(cu_S.getColInd()==NULL || cu_S.getRowPtr()==NULL || cu_S.getValues()==NULL)
            			handleError("faust_gpu2cpu","NULL pointer for non-empty matrix");
		}
			int currentGPU;
            		faust_cudaGetDevice(&currentGPU);
            		faust_cudaSetDevice(cu_S.getDevice());

			FPP1* data_tmp = new FPP1[nnz];
			int* row_ptr_tmp = new int[dim1+1];
			int* id_col_tmp = new int[nnz];

			faust_cudaMemcpyAsync(data_tmp, cu_S.getValues(), nnz*sizeof(FPP1), cudaMemcpyDeviceToHost, stream);
			faust_cudaMemcpyAsync(row_ptr_tmp, cu_S.getRowPtr(), (dim1+1)*sizeof(int), cudaMemcpyDeviceToHost, stream);
            		faust_cudaMemcpyAsync(id_col_tmp, cu_S.getColInd(), nnz*sizeof(int), cudaMemcpyDeviceToHost, stream);




			Faust::MatSparse<FPP,Cpu> spmat_tmp(nnz,dim1,dim2,data_tmp,row_ptr_tmp,id_col_tmp);


			S=spmat_tmp;
			delete [] data_tmp;
			delete [] row_ptr_tmp;
			delete [] id_col_tmp;

			data_tmp=NULL;
			row_ptr_tmp=NULL;
			id_col_tmp=NULL;

      			faust_cudaSetDevice(currentGPU);

	}




}



#endif

#endif
