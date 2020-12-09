##############################################################################
##                              Description:                                ##
##                                                                          ##
##        Cython Header making the links between C++ class FaustCpp         ##
##        and Cython environment, works like C/C++ Header                   ##
##                                                                          ## 
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.inria.fr>                         ##
##                                                                          ##
##                              License:                                    ##
##  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
##                      Luc Le Magoarou, Remi Gribonval                     ##
##                      INRIA Rennes, FRANCE                                ##
##                      http://www.inria.fr/                                ##
##                                                                          ##
##  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
##  General Public License.                                                 ##
##  This program is free software: you can redistribute it and/or modify    ##
##  it under the terms of the GNU Affero General Public License as          ##
##  published by the Free Software Foundation.                              ##
##                                                                          ##
##  This program is distributed in the hope that it will be useful, but     ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
##  See the GNU Affero General Public License for more details.             ##
##                                                                          ##
##  You should have received a copy of the GNU Affero General Public        ##
##  License along with this program.                                        ##
##  If not, see <http://www.gnu.org/licenses/>.                             ##
##                                                                          ##
##                             Contacts:                                    ##
##      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
##      Adrien Leman    : adrien.leman@inria.fr                             ##
##      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##                                                                          ##
##############################################################################

from libcpp cimport bool, complex

cdef extern from "FaustCoreCpp.h":
    cdef cppclass FaustCoreCpp[FPP]:
        FaustCoreCpp()
        void Display() const
        const char* to_string() const
        void push_back(FPP* valueMat,unsigned int nbrow,unsigned int nbcol, bool optimizedCopy)
        void push_back(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows,
                       int ncols, bool optimizedCopy)
        void get_product(FPP* data, int nbrow, int nbcol);
        void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,
                      int nbrow_x, int nbcol_x);#,bool isTranspose*/);
        void set_FM_mul_mode(const int mode);
        void set_Fv_mul_mode(const int mode);
        # Faust-by-csr product -> dense mat
        void multiply(FPP* y_data, int y_nrows, int y_ncols, FPP* x_data, int* x_row_ptr, int* x_id_col, int x_nnz, int x_nrows, int x_ncols);
        unsigned int getNbRow() const
        unsigned int getNbCol() const
        unsigned int getNBytes() const
#        void setOp(const bool isTransposed,unsigned int& nbRowOp, unsigned int& nbColOp)const;
        unsigned long long nnz() const
        double norm(int ord, double threshold, int max_num_its) const
        double normFro() const
        double normInf() const
        FaustCoreCpp[FPP]* normalize(int ord) const
        double get_nb_factors() const
        unsigned int get_fact_nb_rows(unsigned int& i) const
        unsigned int get_fact_nb_cols(unsigned int& i) const
        double get_fact(unsigned int& i, FPP* fact_ptr) const
        void get_fact_sparse(unsigned int& i,
                      int* rowptr,
                      int* col_ids,
                      FPP* elts,
                      const bool transpose) const
        void get_fact_dense(unsigned int& i,
                      FPP* elts,
                      unsigned int*,
                      unsigned int*,
                      const bool transpose) const
        unsigned int get_fact_nnz(const unsigned int i) const
        bool is_fact_sparse(const unsigned int i) const
        FaustCoreCpp[FPP]* right(const unsigned int) const
        FaustCoreCpp[FPP]* left(const unsigned int) const
        FaustCoreCpp[FPP]* get_slice(unsigned int, unsigned int, unsigned int,
                                    unsigned int) const
        FaustCoreCpp[FPP]* slice(unsigned int, unsigned int, unsigned int,
                                    unsigned int) const
        FaustCoreCpp[FPP]* fancy_idx(unsigned long int* row_ids, unsigned long int
                                  num_rows, unsigned long int* col_ids,
                                  unsigned long int num_cols) const
        bool save_mat_file(const char* filepath) const
        FaustCoreCpp[FPP]* swap_cols(const unsigned int id1, const unsigned int id2,
                                const bool permutation, const bool inplace);
        FaustCoreCpp[FPP]* swap_rows(const unsigned int id1, const unsigned int id2,
                                const bool permutation, const bool inplace);

        FaustCoreCpp[FPP]* optimize_time(const bool transp, const bool inplace,
                                       const int nsamples)
        FaustCoreCpp[FPP]* optimize(const bool transp)
        FaustCoreCpp[FPP]* optimize_storage(const bool time)
        const bool isTransposed()
        FaustCoreCpp[FPP]* transpose()
        FaustCoreCpp[FPP]* conjugate()
        FaustCoreCpp[FPP]* adjoint()
        FaustCoreCpp[FPP]* zpruneout(const int nnz_tres, const int npasses,
                                    const bool only_forward)
        FaustCoreCpp[FPP]* vertcat(FaustCoreCpp[FPP]*) const
        FaustCoreCpp[FPP]* vertcatn(FaustCoreCpp[FPP]**, size_t n) const
        FaustCoreCpp[FPP]* horzcat(FaustCoreCpp[FPP]*) const
        FaustCoreCpp[FPP]* horzcatn(FaustCoreCpp[FPP]**, size_t n) const
        FaustCoreCpp[FPP]* mul_faust(FaustCoreCpp[FPP]*)
        FaustCoreCpp[FPP]* mul_scal(FPP scal)
        void device(char*) const
        @staticmethod
        FaustCoreCpp[FPP]* randFaust(int faust_nrows, int faust_ncols,
                                     unsigned int t,
                                     unsigned int min_num_factors, unsigned int max_num_factors,
                                     unsigned int min_dim_size,
                                     unsigned int max_dim_size,
                                     float density, bool per_row)
        @staticmethod
        FaustCoreCpp[FPP]* randFaust(unsigned int t,
                                             unsigned int min_num_factors, unsigned int max_num_factors,
                                             unsigned int min_dim_size,
                                             unsigned int max_dim_size,
                                             float density, bool per_row)
        @staticmethod
        FaustCoreCpp[FPP]* hadamardFaust(unsigned int n, const bool norma)
        @staticmethod
        FaustCoreCpp[FPP]* fourierFaust(unsigned int n, const bool norma)
        @staticmethod
        FaustCoreCpp[FPP]* eyeFaust(unsigned int n, unsigned int m)

# TODO: all the headers below should be in their own pxd file FaustFact.pxd
cdef extern from "FaustFact.h":
    cdef void prox_blockdiag[FPP](FPP* mat_data, unsigned long mat_nrows,
                                  unsigned long mat_ncols, unsigned long
                                  *m_ptr, unsigned long* n_ptr, unsigned long
                                  size, const bool normalized, const bool pos,
                                  FPP* mat_out);

    cdef cppclass PyxConstraintGeneric:
        int name
        unsigned long num_rows
        unsigned long num_cols
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


    cdef FaustCoreCpp[FPP]* fact_palm4MSAFFT[FPP,FPP2](FPP*,unsigned int, unsigned int,
                                             PyxParamsFactPalm4MSAFFT[FPP,FPP2]*,
                                                   FPP*)

    cdef FaustCoreCpp[FPP]* fact_palm4MSA[FPP,FPP2](FPP*,unsigned int, unsigned int,
                                             PyxParamsFactPalm4MSA[FPP,FPP2]*,
                                                   FPP*)

    cdef FaustCoreCpp[FPP]* fact_hierarchical[FPP,FPP2](FPP*,unsigned int, unsigned int,
                                             PyxParamsHierarchicalFact[FPP,FPP2]*,
                                                       FPP*)

    cdef FaustCoreCpp[FPP]* fact_hierarchical_fft[FPP,FPP2](FPP*, FPP*, unsigned int, unsigned int,
                                             PyxParamsHierarchicalFactFFT[FPP,FPP2]*,
                                                       FPP*)

    cdef FaustCoreCpp[FPP]* hierarchical2020[FPP](FPP* mat, unsigned int
                                                  num_rows, unsigned int
                                                  num_cols,
                                                  #unsigned int nites,
                                                  PyxStoppingCriterion*,
                                                  PyxConstraintGeneric**
                                                  constraints, unsigned int
                                                  num_cons, unsigned int
                                                  num_facts, double*
                                                  inout_lambda, bool
                                                  is_update_way_R2L, bool
                                                  is_fact_side_left, bool
                                                  use_csr, bool packing_RL,
                                                  unsigned int norm2_max_iter,
                                                  double norm2_threshold, bool
                                                  is_verbose, bool
                                                  constant_step_size,
                                                  double step_size,
                                                  const bool on_gpu,
                                                  const bool full_gpu)

    cdef FaustCoreCpp[FPP]* palm4msa2020[FPP](FPP* mat,
                                              unsigned int num_rows,
                                              unsigned int num_cols,
                                              PyxConstraintGeneric** constraints,
                                              unsigned int num_cons,
                                              double* inout_lambda,
                                              PyxStoppingCriterion sc,
                                              bool is_update_way_R2L,
                                              bool use_csr,
                                              bool packing_RL,
                                              unsigned int norm2_max_iter,
                                              double norm2_threshold,
                                              bool is_verbose,
                                              bool constant_step_size,
                                              double step_size,
                                              const bool on_gpu,
                                              const bool full_gpu)



cdef extern from "FaustFactGivensFGFT.h":

    cdef FaustCoreCpp[FPP]* fact_givens_fgft[FPP,FPP2](const FPP* Lap, unsigned int num_rows,
                                                       unsigned int num_cols, unsigned int J,
                                                       unsigned int t, FPP* D, unsigned int verbosity,
                                                       const double stoppingError,
                                                       const bool errIsRel,
                                                       const int order,
                                                       const bool enable_large_Faust)

    cdef FaustCoreCpp[FPP]* fact_givens_fgft_sparse[FPP,FPP2](const FPP* data, int* row_ptr,
                                                              int* id_col, int nnz, unsigned int num_rows,
                                                              unsigned int num_cols, unsigned int J,
                                                              unsigned int t, FPP* D,
                                                              unsigned int verbosity,
                                                              const double stoppingError,
                                                              const bool errIsRel,
                                                              const int order,
                                                              const bool enable_large_Faust)

    cdef FaustCoreCpp[FPP]* fact_givens_fgft_cplx[FPP,FPP2](const
                                                            FPP* Lap, unsigned int num_rows,
                                                            unsigned int num_cols, unsigned int J,
                                                            unsigned int t,
                                                            FPP2* D, unsigned int verbosity,
                                                            const double stoppingError,
                                                            const bool errIsRel,
                                                            const int order,
                                                            const bool enable_large_Faust)

    cdef FaustCoreCpp[FPP]* fact_givens_fgft_sparse_cplx[FPP, FPP2](const FPP* data, int* row_ptr,
                                                                    int* id_col, int nnz, unsigned int num_rows,
                                                                    unsigned int num_cols, unsigned int J,
                                                                    unsigned int t,
                                                                    FPP2* D,
                                                                    unsigned int verbosity,
                                                                    const double stoppingError,
                                                                    const bool errIsRel,
                                                                    const int order,
                                                                    const bool enable_large_Faust)

    cdef void svdtj[FPP, FPP2](FaustCoreCpp[FPP]** U, FaustCoreCpp[FPP] **V, FPP* S,
                         const FPP* M_data, unsigned int num_rows, unsigned int
                         num_cols, unsigned int J, unsigned int t, unsigned int
                         verbosity, const double stoppingError, const bool errIsRel, const bool enable_large_Faust)

    cdef void svdtj_sparse[FPP, FPP2](FaustCoreCpp[FPP]** U, FaustCoreCpp[FPP] **V, FPP* S,
                                      const FPP* data, int* row_ptr, int* id_col, int
                                      nnz, int nrows, int ncols, unsigned int J,
                                      unsigned int t, unsigned int verbosity,
                                      const double stoppingError, const bool errIsRel, const bool enable_large_Faust)

    cdef void svdtj_cplx[FPP, FPP2](FaustCoreCpp[FPP]** U, FaustCoreCpp[FPP] **V, FPP* S,
                         const FPP* M_data, unsigned int num_rows, unsigned int
                         num_cols, unsigned int J, unsigned int t, unsigned int
                         verbosity, const double stoppingError, const bool errIsRel, const bool enable_large_Faust)

    cdef void svdtj_sparse_cplx[FPP, FPP2](FaustCoreCpp[FPP]** U, FaustCoreCpp[FPP] **V, FPP* S,
                                      const FPP* data, int* row_ptr, int* id_col, int
                                      nnz, int nrows, int ncols, unsigned int J,
                                      unsigned int t, unsigned int verbosity,
                                      const double stoppingError, const bool errIsRel, const bool enable_large_Faust)

    cdef void* _enable_gpu_mod(const char* libpath, const bool silent)

