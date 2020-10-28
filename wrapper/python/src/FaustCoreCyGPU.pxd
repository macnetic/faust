cdef extern from "FaustCoreCpp.h":
    cdef cppclass FaustCoreCppGPU[FPP](FaustCoreCpp[FPP]):
#        void push_back(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows,
#                       int ncols, bool optimizedCopy)
        void push_back(FPP* valueMat,unsigned int nbrow,unsigned int nbcol, bool optimizedCopy)
        void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,
                      int nbrow_x, int nbcol_x);#,bool isTranspose*/);
        FaustCoreCppGPU[FPP]* normalize_gpu(int ord) const

        @staticmethod
        FaustCoreCppGPU[FPP]* randFaustGPU(unsigned int t,
                                        unsigned int min_num_factors,
                                        unsigned int max_num_factors,
                                        unsigned int min_dim_size,
                                        unsigned int max_dim_size,
                                        float density, bool per_row)
        @staticmethod
        FaustCoreCppGPU[FPP]* hadamardFaustGPU(unsigned int n, const bool norma)
        @staticmethod
        FaustCoreCppGPU[FPP]* fourierFaustGPU(unsigned int n, const bool norma)
        @staticmethod
        FaustCoreCppGPU[FPP]* eyeFaustGPU(unsigned int n, unsigned int m)

        FaustCoreCppGPU[FPP]* mul_faust_gpu(FaustCoreCppGPU[FPP]*)
        FaustCoreCppGPU[FPP]* mul_scal_gpu(const FPP &scal)
