#ifndef __FAUST_CUDA_HPP_
#define __FAUST_CUDA_HPP_

inline void faust_cu_csr2coo(cusparseHandle_t handle, const int *csrRowPtr,
   int nnz, int m, int *cooRowInd, cusparseIndexBase_t idxBase)
{faust_cusparseXcsr2coo(handle, csrRowPtr, nnz, m, cooRowInd, idxBase);}

#endif
