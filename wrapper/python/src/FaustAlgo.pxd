cdef extern from "FaustCoreCpp.h":
    void polyCoeffs[FPP](int d, int K, int n, const FPP* basisX, const FPP*
                         coeffs, FPP* out, bool on_gpu)
    void polyGroupCoeffs_[FPP](int d, int K, int n, const FPP* basisX, const FPP* coeffs, FPP** out, int n_out, bool on_gpu);


cdef extern from "FaustFact.h":
    cdef void prox_blockdiag[FPP](FPP* mat_data, unsigned long mat_nrows,
                                  unsigned long mat_ncols, unsigned long
                                  *m_ptr, unsigned long* n_ptr, unsigned long
                                  size, const bool normalized, const bool pos,
                                  FPP* mat_out);

    cdef cppclass PyxMHTPParams[FPP]:
        bool used
        PyxStoppingCriterion[FPP] stop_crit
        bool constant_step_size
        FPP step_size
        int palm4msa_period
        bool updating_lambda

    cdef cppclass PyxConstraintGeneric:
        int name
        unsigned long num_rows
        unsigned long num_cols
        bool normalizing
        bool pos
        bool is_int_constraint()
        bool is_real_constraint()
        bool is_mat_constraint()

    cdef cppclass PyxConstraintInt(PyxConstraintGeneric):
        unsigned long parameter

    int prox_int[FPP](unsigned int cons_type, unsigned long cons_param, FPP* mat_in, unsigned long num_rows,
                  unsigned long num_cols, FPP* mat_out, const bool normalized, const bool pos)

    int prox_real[FPP,FPP2](unsigned int cons_type, FPP2 cons_param, FPP* mat_in, unsigned long num_rows,
                  unsigned long num_cols, FPP* mat_out, const bool normalized, const bool pos)

    int prox_mat[FPP](unsigned int cons_type, FPP* cons_param, unsigned long cons_param_sz, FPP* mat_in, unsigned long num_rows,
                  unsigned long num_cols, FPP* mat_out, const bool normalized, const bool pos)


    cdef cppclass PyxConstraintScalar[FPP](PyxConstraintGeneric):
        FPP parameter

    cdef cppclass PyxConstraintMat[FPP](PyxConstraintGeneric):
        FPP* parameter # shape = num_rows, num_cols (except if BLOCKDIAG constraint)
        unsigned long parameter_sz

    cdef cppclass PyxStoppingCriterion[FPP]:
        bool is_criterion_error
        int num_its
        FPP error_threshold
        unsigned long max_num_its
        FPP erreps

    cdef cppclass PyxParamsFact[FPP,FPP2]:
        int num_facts
        bool is_update_way_R2L
        FPP2 init_lambda
        FPP2 step_size
        PyxConstraintGeneric** constraints # (num_facts-1)*2 elts
        unsigned int num_constraints
        bool is_verbose
        bool constant_step_size
        unsigned int grad_calc_opt_mode
        int norm2_max_iter
        double norm2_threshold

    cdef cppclass PyxParamsFactPalm4MSA[FPP,FPP2](PyxParamsFact[FPP,FPP2]):
        FPP** init_facts # num_facts elts
        unsigned long* init_fact_sizes
        PyxStoppingCriterion[FPP2] stop_crit

    cdef cppclass \
    PyxParamsFactPalm4MSAFFT[FPP,FPP2](PyxParamsFactPalm4MSA[FPP,FPP2]):
        FPP* init_D

    cdef cppclass PyxParamsHierarchicalFact[FPP,FPP2](PyxParamsFact[FPP,FPP2]):
        unsigned int num_rows
        unsigned int num_cols
        PyxStoppingCriterion[FPP2]* stop_crits #must be of size 2
        bool is_fact_side_left

    cdef cppclass PyxParamsHierarchicalFactFFT[FPP,FPP2](PyxParamsHierarchicalFact[FPP,FPP2]):
        FPP* init_D


    cdef FaustCoreCppCPU[FPP]* fact_palm4MSAFFT[FPP,FPP2](FPP*,unsigned int, unsigned int,
                                             PyxParamsFactPalm4MSAFFT[FPP,FPP2]*,
                                                   FPP*)

    cdef FaustCoreCppCPU[FPP]* fact_palm4MSA[FPP,FPP2](FPP*,unsigned int, unsigned int,
                                             PyxParamsFactPalm4MSA[FPP,FPP2]*,
                                                   FPP*)

    cdef FaustCoreCppCPU[FPP]* fact_hierarchical[FPP,FPP2](FPP*,unsigned int, unsigned int,
                                             PyxParamsHierarchicalFact[FPP,FPP2]*,
                                                       FPP*)

    cdef FaustCoreCppCPU[FPP]* fact_hierarchical_fft[FPP,FPP2](FPP*, FPP*, unsigned int, unsigned int,
                                             PyxParamsHierarchicalFactFFT[FPP,FPP2]*,
                                                       FPP*)

    cdef FaustCoreCppCPU[FPP]* hierarchical2020[FPP, FPP2](FPP* mat, unsigned int
                                                  num_rows, unsigned int
                                                  num_cols,
                                                  #unsigned int nites,
                                                  PyxStoppingCriterion[FPP2]*,
                                                  PyxConstraintGeneric**
                                                  constraints, unsigned int
                                                  num_cons, unsigned int
                                                  num_facts, FPP2*
                                                  inout_lambda, bool
                                                  is_update_way_R2L, bool
                                                  is_fact_side_left, int
                                                  factor_format, bool packing_RL,
                                                  bool no_normalization,
                                                  bool no_lambda,
                                                  PyxMHTPParams[FPP2] mhtpp,
                                                  unsigned int norm2_max_iter,
                                                  FPP2 norm2_threshold, bool
                                                  is_verbose, bool
                                                  constant_step_size,
                                                  FPP2 step_size,
                                                  bool full_gpu)

    cdef FaustCoreCppCPU[FPP]* palm4msa2020[FPP, FPP2](FPP* mat,
                                              unsigned int num_rows,
                                              unsigned int num_cols,
                                              PyxConstraintGeneric** constraints,
                                              unsigned int num_cons,
                                              FPP2* inout_lambda,
                                              PyxStoppingCriterion[FPP2] sc,
                                              bool is_update_way_R2L,
                                              int factor_format,
                                              bool packing_RL,
                                              bool no_normalization,
                                              bool no_lambda,
                                              PyxMHTPParams[FPP2] mhtpp,
                                              unsigned int norm2_max_iter,
                                              FPP2 norm2_threshold,
                                              bool is_verbose,
                                              bool constant_step_size,
                                              FPP2 step_size,
                                              bool full_gpu,
                                              FaustCoreCppCPU[FPP]* cth)

    cdef FaustCoreCppCPU[FPP]* butterfly_hierarchical[FPP](FPP* mat, unsigned int,
                                                        unsigned int, int dir)
    cdef FaustCoreCppCPU[FPP]* butterfly_hierarchical[FPP](FPP* mat, unsigned int,
                                                        unsigned int, int dir,
                                                          int* perm, bool
                                                           mul_perm)


cdef extern from "FaustFactGivensFGFT.h":

    cdef FaustCoreCppCPU[FPP]* fact_givens_fgft[FPP,FPP2](const FPP* Lap, unsigned int num_rows,
                                                       unsigned int num_cols, unsigned int J,
                                                       unsigned int t, FPP* D, unsigned int verbosity,
                                                       const FPP2 stoppingError,
                                                       const bool errIsRel,
                                                       const int order,
                                                       const bool enable_large_Faust)

    cdef FaustCoreCppCPU[FPP]* fact_givens_fgft_sparse[FPP,FPP2](const FPP* data, int* row_ptr,
                                                              int* id_col, int nnz, unsigned int num_rows,
                                                              unsigned int num_cols, unsigned int J,
                                                              unsigned int t, FPP* D,
                                                              unsigned int verbosity,
                                                              const FPP2 stoppingError,
                                                              const bool errIsRel,
                                                              const int order,
                                                              const bool enable_large_Faust)

    cdef FaustCoreCppCPU[FPP]* fact_givens_fgft_cplx[FPP,FPP2](const
                                                            FPP* Lap, unsigned int num_rows,
                                                            unsigned int num_cols, unsigned int J,
                                                            unsigned int t,
                                                            FPP2* D, unsigned int verbosity,
                                                            const FPP2 stoppingError,
                                                            const bool errIsRel,
                                                            const int order,
                                                            const bool enable_large_Faust)

    cdef FaustCoreCppCPU[FPP]* fact_givens_fgft_sparse_cplx[FPP, FPP2](const FPP* data, int* row_ptr,
                                                                    int* id_col, int nnz, unsigned int num_rows,
                                                                    unsigned int num_cols, unsigned int J,
                                                                    unsigned int t,
                                                                    FPP2* D,
                                                                    unsigned int verbosity,
                                                                    const FPP2 stoppingError,
                                                                    const bool errIsRel,
                                                                    const int order,
                                                                    const bool enable_large_Faust)

    cdef void svdtj[FPP, FPP2](FaustCoreCppCPU[FPP]** U, FaustCoreCppCPU[FPP] **V, FPP* S,
                               const FPP* M_data,
                               unsigned int num_rows,
                               unsigned int num_cols,
                               unsigned int J1, unsigned int J2,
                               unsigned int t1, unsigned int t2,
                               unsigned int verbosity,
                               const FPP2 stoppingError,
                               const bool errIsRel,
                               const bool enable_large_Faust,
                               const int err_period)

    cdef void svdtj_sparse[FPP, FPP2](FaustCoreCppCPU[FPP]** U, FaustCoreCppCPU[FPP] **V, FPP* S,
                                      const FPP* data,
                                      int* row_ptr,
                                      int* id_col,
                                      int nnz,
                                      int nrows,
                                      int ncols,
                                      unsigned int J1, unsigned int J2,
                                      unsigned int t1, unsigned int t2,
                                      unsigned int verbosity,
                                      const FPP2 stoppingError,
                                      const bool errIsRel,
                                      const bool enable_large_Faust,
                                      const int err_period)

    cdef void svdtj_cplx[FPP, FPP2](FaustCoreCppCPU[FPP]** U, FaustCoreCppCPU[FPP] **V, FPP* S,
                                    const FPP* M_data,
                                    unsigned int num_rows,
                                    unsigned int num_cols,
                                    unsigned int J1, unsigned int J2,
                                    unsigned int t1, unsigned int t2,
                                    unsigned int verbosity,
                                    const FPP2 stoppingError,
                                    const bool errIsRel,
                                    const bool enable_large_Faust,
                                    const int err_period)

    cdef void svdtj_sparse_cplx[FPP, FPP2](FaustCoreCppCPU[FPP]** U, FaustCoreCppCPU[FPP] **V, FPP* S,
                                           const FPP* data,
                                           int* row_ptr,
                                           int* id_col,
                                           int nnz,
                                           int nrows,
                                           int ncols,
                                           unsigned int J1, unsigned int J2,
                                           unsigned int t1, unsigned int t2,
                                           unsigned int verbosity,
                                           const FPP2 stoppingError,
                                           const bool errIsRel,
                                           const bool enable_large_Faust,
                                           const int err_period)
