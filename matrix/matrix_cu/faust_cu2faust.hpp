#ifndef __FAUST_CU2FAUST_HPP__
#define __FAUST_CU2FAUST_HPP__

#include "faust_vec.h"
#include "faust_mat.h"
#include "faust_cu_vec.h"
#include "faust_cu_mat.h"

#ifdef __COMPILE_SPMAT__
   #include "faust_spmat.h"
   #include "faust_cu_spmat.h"
#endif


const char * class_name="faust_cu2faust";

template<typename T, typename U>
void faust_cu2faust(faust_vec<T>& v, const faust_cu_vec<U>& cu_v, cudaStream_t stream/*=0*/ )
{
   if(cu_v.getData() != NULL)
   {
      v.resize(cu_v.size());
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(cu_v.getDevice());

      if(typeid(T)!=typeid(U))
      {
         U* data_tmp = new U[cu_v.size()];
         faust_cudaMemcpyAsync(data_tmp, cu_v.getData(), cu_v.size()*sizeof(U), cudaMemcpyDeviceToHost, stream);
	      for (int i=0 ; i<v.size() ;i++)
            *(v.getData(i)) = (T) data_tmp[i];
	      delete [] data_tmp ; data_tmp=NULL;	
      }
      else
      {
         faust_cudaMemcpyAsync(v.getData(), cu_v.getData(), cu_v.size()*sizeof(T), cudaMemcpyDeviceToHost, stream);
      }



      faust_cudaSetDevice(currentGPU);
   }
   else if (cu_v.size()==0)
      v.resize(0);   
   else
      handleError(class_name,"NULL pointer for non-empty vector");	
      
}

template<typename T, typename U>
void faust_cu2faust(faust_mat<T>& M, const faust_cu_mat<U>& cu_M, cudaStream_t stream/*=0*/)
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


      if(typeid(T)!=typeid(U))
      {
   	   U* data_tmp = new U[dim1*dim2];
         faust_cudaMemcpyAsync(data_tmp, cu_M.getData(), dim1*dim2*sizeof(U), cudaMemcpyDeviceToHost, stream);      
	      for (int i=0 ; i<dim1*dim2 ;i++)
            *(M.getData(i)) = (T) data_tmp[i];
	      delete [] data_tmp ; data_tmp=NULL;
      }
      else
         faust_cudaMemcpyAsync(M.getData(), cu_M.getData(), dim1*dim2*sizeof(T), cudaMemcpyDeviceToHost, stream);      

      faust_cudaSetDevice(currentGPU);
   }
   else if(dim1*dim2 != 0)
      handleError(class_name,"NULL pointer for non-empty matrix");	
}

#ifdef __COMPILE_SPMAT__
template<typename T, typename U>
void faust_cu2faust(faust_spmat<T>& S, const faust_cu_spmat<U>& cu_S, cudaStream_t stream/*=0*/)
{
      const int nnz = cu_S.getNonZeros();
      const int dim1 = cu_S.getNbRow();
      const int dim2 = cu_S.getNbCol();
      S.resize(nnz,dim1,dim2);           
      if(dim1*dim2>0 && nnz==0)
         S.setZeros();

      if(dim1*dim2>0 && nnz>0)
      {
         if(cu_S.getColInd()==NULL || cu_S.getRowPtr()==NULL || cu_S.getValues()==NULL
            || S.getColInd()==NULL || S.getRowPtr()==NULL || S.getValuePtr()==NULL)
            
            handleError(class_name,"NULL pointer for non-empty matrix");

         else
         {
            int currentGPU;
            faust_cudaGetDevice(&currentGPU);
            faust_cudaSetDevice(cu_S.getDevice());
            
            //if eigen sparse matrix S.mat is not in rowMajor
               //handleError(class_name,"Eigen sparse matrix must be in RowMajor");


            faust_cudaMemcpyAsync(S.getRowPtr(), cu_S.getRowPtr(), (dim1+1)*sizeof(int), cudaMemcpyDeviceToHost, stream);      
            faust_cudaMemcpyAsync(S.getColInd(), cu_S.getColInd(), nnz*sizeof(int), cudaMemcpyDeviceToHost, stream);      


	         if (typeid(T) != typeid(U))
	         {
               U* data_tmp = new U[nnz];
               faust_cudaMemcpyAsync(data_tmp, cu_S.getValues(), nnz*sizeof(U), cudaMemcpyDeviceToHost, stream);
               for (int i=0 ; i<nnz ;i++)
                  S.getValuePtr()[i] = (T) data_tmp[i];
               delete [] data_tmp ; data_tmp=NULL;
	         }
	         else
               faust_cudaMemcpyAsync(S.getValuePtr(), cu_S.getValues(), nnz*sizeof(T), cudaMemcpyDeviceToHost, stream);      

            faust_cudaSetDevice(currentGPU);
         }                  
      }  
}
#endif

#endif
