#include "kernels.h"
#include "faust_cuda.h"

//prototypes
template<typename T> __global__ void Add_inria(T *A, const T *B, int numElements);
template<typename T> __global__ void Sub_inria(T *A, const T *B, int numElements);
template<typename T> __global__ void Mult_inria(T *A, const T *B, int numElements);
template<typename T> __global__ void Div_inria(T *A, const T *B, int numElements);

template<typename T> __global__ void AddConst_inria(T *A, T val, int numElements);
template<typename T> __global__ void SubConst_inria(T *A, T val, int numElements);
template<typename T> __global__ void MultConst_inria(T *A, T val, int numElements);
template<typename T> __global__ void DivConst_inria(T *A, T val, int numElements);

template<typename T> __global__ void Square_inria(T *A, int numElements);
template<typename T> __global__ void Square_inria(T *d_cu_dst, const T* d_cu_src, int numElements);
template<typename T> __global__ void Sqrt_inria(T *A, int numElements);
template<typename T> __global__ void Inv_inria(T *A, int numElements);
template<typename T> __global__ void Abs_inria(T *A, int numElements);

template<typename T, typename U> __global__ void Memcpy_inria(T* d_cu_dst, const U* d_cu_src, int length);


template<typename T> __global__ void Memset_inria(T* dev_dst, T valeur, int numElements);

template<typename T> __global__ void Sparse2full_inria(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1);
template<typename T> __global__ void AddSparse2full_inria(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1);
template<typename T> __global__ void SubSparse2full_inria(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1);
template<typename T> __global__ void GetDiag_inria(T* dst_diag, const T* src_M,  int dim1);
template<typename T> __global__ void AddDiagConst_inria(T* dev_dst, T val, int dim1);
template<typename T> __global__ void CopyDiag_inria(T* dev_dst, T* dev_src, int dim1);

template<typename T> __global__ void SetSubmatrix_inria(T* mat_dst, T* mat_src, int src_dim1, int r1, int c1, int nb_rows, int nb_elements);


template<typename T> __global__ void RelativeError_inria(T* data_dst, const T* data_src_th, const T* data_src_mes, const int length);






// Kernel wrappers
template<typename T> void kernel_add(T* d_cu1, const T* d_cu2, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Add_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, d_cu2, length);
   faust_kernelSafe();
}
template<typename T> void kernel_sub(T* d_cu1, const T* d_cu2, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Sub_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, d_cu2, length);
   faust_kernelSafe();
}
template<typename T> void kernel_mult(T* d_cu1, const T* d_cu2, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Mult_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, d_cu2, length);
   faust_kernelSafe();
}
template<typename T>  void kernel_div(T* d_cu1, const T* d_cu2, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Div_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, d_cu2, length);
   faust_kernelSafe();
}
template<typename T> void kernel_add_const(T* d_cu1, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	AddConst_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, valeur, length);
   faust_kernelSafe();
}
template<typename T> void kernel_sub_const(T* d_cu1, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	SubConst_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, valeur, length);
   faust_kernelSafe();
}
template<typename T> void kernel_mult_const(T* d_cu1, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	MultConst_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, valeur, length);
   faust_kernelSafe();
}
template<typename T>  void kernel_div_const(T* d_cu1, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	DivConst_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, valeur, length);
   faust_kernelSafe();
}
template<typename T> void kernel_square(T* d_cu1, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Square_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, length);
   faust_kernelSafe();
}
template<typename T> void kernel_square(T* d_cu_dst, const T* d_cu_src, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Square_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu_dst, d_cu_src, length);
   faust_kernelSafe();
}
template<typename T> void kernel_sqrt(T* d_cu1, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Sqrt_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, length);
   faust_kernelSafe();
}

template<typename T> void kernel_inv(T* d_cu1, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Inv_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, length);
   faust_kernelSafe();
}

template<typename T> void kernel_abs(T* d_cu1, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Abs_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, length);
   faust_kernelSafe();
}



template<typename T, typename U> void kernel_memcpy(T* d_cu_dst, const U* d_cu_src, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Memcpy_inria<T,U><<<blocksPerGrid,threadsPerBlock>>>(d_cu_dst, d_cu_src, length);
   faust_kernelSafe();
}


template<typename T> void kernel_memset(T* d_cu_dst, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Memset_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu_dst, valeur, length);
   faust_kernelSafe();
}	
template<typename T> void kernel_sparse2full(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1, const int src_dim2)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nnz + threadsPerBlock - 1) / threadsPerBlock;
	T zero = (T) 0.0;
	kernel_memset<T>(dev_dst, zero, src_dim1*src_dim2);
	Sparse2full_inria<T><<<blocksPerGrid,threadsPerBlock>>>(dev_dst, dev_src_rowind, dev_src_colind, dev_src_values, nnz, src_dim1);
   faust_kernelSafe();
}
template<typename T> void kernel_add_sparse2full(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nnz + threadsPerBlock - 1) / threadsPerBlock;
	AddSparse2full_inria<T><<<blocksPerGrid,threadsPerBlock>>>(dev_dst, dev_src_rowind, dev_src_colind, dev_src_values, nnz, src_dim1);
   faust_kernelSafe();
}
template<typename T> void kernel_sub_sparse2full(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nnz + threadsPerBlock - 1) / threadsPerBlock;
	SubSparse2full_inria<T><<<blocksPerGrid,threadsPerBlock>>>(dev_dst, dev_src_rowind, dev_src_colind, dev_src_values, nnz, src_dim1);
   faust_kernelSafe();
}
template<typename T> void kernel_get_diag(T* dst_diag, const T* src_M, int dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (dim1 + threadsPerBlock - 1) / threadsPerBlock;
	GetDiag_inria<T><<<blocksPerGrid,threadsPerBlock>>>(dst_diag, src_M, dim1);
   faust_kernelSafe();
}
template<typename T> void kernel_add_diag_const(T* d_cu1, T val, int dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (dim1 + threadsPerBlock - 1) / threadsPerBlock;
	AddDiagConst_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, val, dim1);
   faust_kernelSafe();
}
template<typename T> void kernel_copy_diag(T* d_cu_dst, T* d_cu_src, int dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (dim1 + threadsPerBlock - 1) / threadsPerBlock;
	CopyDiag_inria<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu_dst, d_cu_src, dim1);
   faust_kernelSafe();
}
template<typename T> void kernel_set_submatrix(T* mat_dst, T* mat_src, int src_dim1, int r1, int c1, int nb_rows, int nb_col)
{
	int threadsPerBlock = 256;
	int nb_elements = nb_rows*nb_col;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	SetSubmatrix_inria<T><<<blocksPerGrid,threadsPerBlock>>>(mat_dst, mat_src, src_dim1, r1, c1, nb_rows, nb_elements);
   faust_kernelSafe();
}

template<typename T> void kernel_relative_error(T* data_dst, const T* data_src_th, const T* data_src_mes, const int nb_elements)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	RelativeError_inria<T><<<blocksPerGrid,threadsPerBlock>>>(data_dst, data_src_th, data_src_mes, nb_elements);
   faust_kernelSafe();
}






// Kernels

template<typename T> __global__ void Add_inria(T* A, const T* B, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = A[i] + B[i];
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void Sub_inria(T* A, const T* B, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    
    while (i < numElements)
    {
        A[i] = A[i] - B[i];
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void Mult_inria(T *A, const T *B, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = A[i] * B[i];
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void Div_inria(T *A, const T *B, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = A[i] / B[i];
	i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void AddConst_inria(T *A, T val, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = A[i] + val;
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void SubConst_inria(T *A, T val, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = A[i] - val;
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void MultConst_inria(T *A, T val, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = A[i] * val;
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void DivConst_inria(T* A, T val, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = A[i] / val;
	i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void Square_inria(T *A, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = A[i]*A[i]; 
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void Square_inria(T* dst, const T* src, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        dst[i] = src[i]*src[i]; 
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void Sqrt_inria(T *A, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = sqrt(A[i]); 
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void Inv_inria(T *A, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = 1/A[i] ;
	i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void Abs_inria(T *A, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        A[i] = fabs(A[i]) ;
	i += gridDim.x * blockDim.x;
    }
}

template<typename T, typename U> __global__ void Memcpy_inria(T* dev_dst, const U* dev_src, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        dev_dst[i] = (T)dev_src[i] ;
	i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void Memset_inria(T* dev_dst, T valeur, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < numElements)
    {
        dev_dst[i] = valeur ;
	i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void Sparse2full_inria(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < nnz)
    {
        int j = (dev_src_colind[i]) * src_dim1 + (dev_src_rowind[i]);
        dev_dst[j] =  dev_src_values[i];
	i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void AddSparse2full_inria(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < nnz)
    {
        int j = (dev_src_colind[i]) * src_dim1 + (dev_src_rowind[i]);
        dev_dst[j] +=  dev_src_values[i];
        i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void SubSparse2full_inria(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < nnz)
    {
        int j = (dev_src_colind[i]) * src_dim1 + (dev_src_rowind[i]);
        dev_dst[j] -=  dev_src_values[i];
        i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void GetDiag_inria(T* dst_diag, const T* src_M, int dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < dim1)
    {
        int j = i*(dim1+1);
        dst_diag[i] = src_M[j];
	i += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void AddDiagConst_inria(T *A, T val, int dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < dim1)
    {
        int j = i*(dim1+1);
        A[j] = A[j] + val;
	i += gridDim.x * blockDim.x;
    }
}


template<typename T> __global__ void CopyDiag_inria(T *dst, T* src, int dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    while (i < dim1)
    {
        int j = i*(dim1+1);
        dst[j] = src[j];
	i += gridDim.x * blockDim.x;
    }
}
template<typename T> __global__ void SetSubmatrix_inria(T* mat_dst, T* mat_src, int src_dim1, int r1, int c1, int nb_rows, int nb_elements)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    while(ind<nb_elements)
    {
	    int col = (int)(ind/nb_rows) ;
            int ligne = ind - nb_rows*col;
	    mat_dst[ind] = mat_src[(col+c1)*src_dim1+(ligne+r1)];        
	    ind += gridDim.x * blockDim.x;
    }
}

template<typename T> __global__ void RelativeError_inria(T* data_dst, const T* data_src_th, const T* data_src_mes, const int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    while(i<numElements)
    {
            data_dst[i] = fabs((data_src_th[i]-data_src_mes[i])/data_src_th[i]);
	    i += gridDim.x * blockDim.x;
    }
}


/*template<typename T> __global__ void SetLinfAR1(T* mat_Linf, T* mat_Hinf, T* atur, int nb_n, int nb_elements)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    while(ind<nb_elements)
    {
        int col = (int)(ind/nb_n) ;
        int ligne = ind - nb_n*col;
	int nb_az = nb/2;
	if(ligne < nb_az)
	{
             mat_Linf[ind] = atur[ligne] * mat_Hinf[ind];
	}
	else
	{
	     mat_Linf[ind] = mat_Hinf[col*nb_n + ligne-nb_az];
	}
	ind += gridDim.x * blockDim.x;
    }
}*/


#include "kernel_def.hu"
