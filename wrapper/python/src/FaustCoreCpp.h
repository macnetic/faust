/****************************************************************************/
/*                              Description:                                */
/*  C++ Header file wrapping the  Class Faust::Transform<FPP,Cpu>,          */
/*  the methods are quite the same but this works with pointer              */
/*  which are easier to handle from Python                                  */
/*                                                                          */     
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
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

#ifndef FAUSTCORECPP_H
#define FAUSTCORECPP_H

#include "faust_MatDense.h"
#ifdef USE_GPU_MOD
#include "faust_TransformHelper_gpu.h"
#endif
#include "faust_TransformHelper.h"
#include "faust_TransformHelperPoly.h"
#include <cstring>
#include <complex>

template<typename FPP, FDevice DEV> class FaustCoreCpp;

#if USE_GPU_MOD
template<typename FPP>
FaustCoreCpp<FPP, Cpu>* clone_gpu2cpu(FaustCoreCpp<FPP, GPU2>* gpu_t);
template<typename FPP>
FaustCoreCpp<FPP, GPU2>* clone_cpu2gpu(FaustCoreCpp<FPP, Cpu>* cpu_t);
#endif

template<typename FPP, FDevice DEV=Cpu>
class FaustCoreCpp
{


    public :


    FaustCoreCpp(): transform(nullptr) {}
    FaustCoreCpp(Faust::TransformHelper<FPP,DEV> *th);
    /**
     * Copies the underlying transform into the argument.
     *
     * \return true to indicate the two objects are different.
     */
#ifdef USE_GPU_MOD
    bool make_transform(Faust::TransformHelper<FPP, GPU2>**) const;
    /**
     * Copies the pointer of the underlying transform into the argument.
     *
     * \return false to indicate the two objects are the same.
     */
#endif
    bool make_transform(Faust::TransformHelper<FPP, Cpu>**) const;
    void Display() const { transform->display();}
    const char* to_string() const;
    void push_back(FPP* valueMat,unsigned int nbrow,unsigned int nbcol, bool optimizedCopy=false);
    void push_back(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, bool optimizedCopy=false);
    void push_back(FPP* bdata, int* brow_ptr, int* bcol_inds, int nrows, int ncols, int bnnz, int bnrows, int bncols, bool optimizedCopy=false);
    unsigned int getNBytes() const;
    unsigned int getNbRow() const;
    unsigned int getNbCol() const;
    void get_product(FPP* y_data, int y_nrows, int y_ncols);
    void multiply(FPP* y_data, int y_nrows, int y_ncols, FPP* x_data, int* x_row_ptr, int* x_id_col, int x_nnz, int x_nrows, int x_ncols);
    void set_FM_mul_mode(const int mode, const bool silent=true);
    void multiply(FPP* value_y,int nbrow_y,int nbcol_y,const FPP* value_x,int nbrow_x,int nbcol_x/*,bool isTranspose*/)const;
    FaustCoreCpp<FPP,DEV>* mul_faust(FaustCoreCpp<FPP,DEV>* right);
    FaustCoreCpp<FPP,DEV>* vertcat(FaustCoreCpp<FPP,DEV>* right) const;
    FaustCoreCpp<FPP,DEV>* vertcatn(FaustCoreCpp<FPP,DEV>** rights, size_t n) const;
    FaustCoreCpp<FPP,DEV>* horzcat(FaustCoreCpp<FPP,DEV>* right) const;
    FaustCoreCpp<FPP,DEV>* horzcatn(FaustCoreCpp<FPP,DEV>** rights, size_t n) const;
    FaustCoreCpp<FPP,DEV>* mul_scal(FPP scal);
    FaustCoreCpp<FPP,DEV>* normalize(int ord) const;
    unsigned long long nnz()const;
    void power_iteration(FPP* out, double threshold=.001, int max_num_its=100) const;
    double norm2(double threshold=.001, int max_num_its=100) const;
    double normFro(const bool full_array, const int batch_size) const;
    double normInf(const bool full_array, const int batch_size) const;
    double norm1(const bool full_array, const int batch_size) const;
    double get_nb_factors() const;
    unsigned int get_fact_nb_rows(unsigned int& i) const;
    unsigned int get_fact_nb_cols(unsigned int& i) const;
    bool is_all_sparse(bool csr, bool bsr) const;
    bool is_all_dense() const;
    void get_fact(const unsigned int& i, FPP* fact_ptr) const;
    void get_fact_sparse(const unsigned int& i,
            int* rowptr,
            int* col_ids,
            FPP* elts,
            const bool transpose = false) const;
    void get_fact_dense(const unsigned int& i, FPP* elts,
            unsigned int* num_rows, unsigned int* num_cols,
            const bool transpose) const;
    void get_fact_bsr_info(const faust_unsigned_int id,
            size_t& bdata_sz,
            size_t& browptr_sz,
            size_t& bcolinds_sz,
            size_t& bnnz,
            size_t& bnrows,
            size_t& bncols ) const;
    void get_fact_bsr(const faust_unsigned_int id,
            FPP* bdata,
            int* brow_ptr,
            int* bcol_inds) const;
    FaustCoreCpp<FPP,DEV>* left(const faust_unsigned_int) const;
    FaustCoreCpp<FPP,DEV>* right(const faust_unsigned_int) const;
    faust_unsigned_int get_fact_nnz(const faust_unsigned_int) const;
    bool is_fact_sparse(const faust_unsigned_int id) const;
    int get_fact_type(const faust_unsigned_int id) const;
    FaustCoreCpp<FPP,DEV>* slice(unsigned int, unsigned int, unsigned int, unsigned int) const;
    FaustCoreCpp<FPP,DEV>* fancy_idx(unsigned long int* row_ids, unsigned long int
                                  num_rows, unsigned long int* col_ids,
                                  unsigned long int num_cols) const;
    bool save_mat_file(const char* filepath) const;
    static FaustCoreCpp<FPP, DEV>* read_from_mat_file(const char* filepath);
    static int get_mat_file_type(const char* filepath);
    FaustCoreCpp<FPP,DEV>* swap_cols(const unsigned int id1, const unsigned int id2,
            const bool permutation, const bool inplace);
    FaustCoreCpp<FPP,DEV>* swap_rows(const unsigned int id1, const unsigned int id2,
            const bool permutation, const bool inplace);

    FaustCoreCpp<FPP,DEV>* optimize_storage(const bool time=false);
    FaustCoreCpp<FPP,DEV>* optimize(const bool transp=false);
    FaustCoreCpp<FPP,DEV>* optimize_time(const bool transp=false, const bool inplace=false, const int nsamples=1);
    FaustCoreCpp<FPP,DEV>* optimize_time(const FPP* value_x,int nbrow_x,int nbcol_x, const bool transp=false, const bool inplace=false, const int nsamples=1);
    FaustCoreCpp<FPP,DEV>* optimize_time(const FPP* x_data, int* x_row_ptr, int* x_id_col, int x_nnz, int x_nrows, int x_ncols, const bool transp=false, const bool inplace=false, const int nsamples=1);

    const bool isTransposed();
    FaustCoreCpp<FPP,DEV>* transpose()const;
    FaustCoreCpp<FPP,DEV>* conjugate()const;
    FaustCoreCpp<FPP,DEV>* adjoint()const;
    FaustCoreCpp<FPP,DEV>* zpruneout(const int nnz_tres, const int npasses, const bool only_forward);
    FaustCoreCpp<FPP,DEV>* clone(int dev=-1) const; //-2 for CPU
    void polyCoeffs(int d, int n, const FPP* basisX, const FPP* coeffs, FPP* out, bool on_gpu) const;
    FaustCoreCpp<FPP,DEV>* polyCoeffs(const FPP* coeffs);
    void mulPolyCoeffs(const FPP* X, int n, FPP* Y, const FPP* coeffs);
    FaustCoreCpp<FPP,DEV>* polyNext() const;
	void device(char* dev) const;
    FaustCoreCpp<Real<FPP>,DEV>* real();
    FPP get_item(unsigned long int i, unsigned long int j);
    void colSliceMultiply(unsigned long int j1, unsigned long int j2, const FPP* data, unsigned long int data_ncols, FPP* out) const;
    void indexMultiply(unsigned long int* d0_ids, size_t d0_ids_len, unsigned long int* d1_ids, size_t d1_ids_len, const FPP* mat_data, int mat_ncols, FPP* mat_out) const;
    ~FaustCoreCpp();
    static FaustCoreCpp<FPP,DEV>* randFaust(unsigned int t,
            unsigned int min_num_factors, unsigned int max_num_factors,
            unsigned int min_dim_size,
            unsigned int max_dim_size, float density, bool per_row);
    static FaustCoreCpp<FPP,DEV>* randFaust(int faust_nrows, int faust_ncols,
            unsigned int t,
            unsigned int min_num_factors, unsigned int max_num_factors,
            unsigned int min_dim_size,
            unsigned int max_dim_size, float density, bool per_row);
    static FaustCoreCpp<FPP,DEV>* randBSRFaust(unsigned int faust_nrows, unsigned int faust_ncols,
            unsigned int min_num_factors, unsigned int max_num_factors,
            unsigned int bnrows, unsigned int bncols, float density);

    static FaustCoreCpp<FPP,DEV>* hadamardFaust(unsigned int n, const bool norma);
    static FaustCoreCpp<FPP,DEV>* fourierFaust(unsigned int n, const bool norma);
    static FaustCoreCpp<FPP,DEV>* eyeFaust(unsigned int n, unsigned int m);
    static FaustCoreCpp<FPP, DEV>* polyBasis(
            unsigned int L_nrows,
            unsigned int L_ncols,
            int* L_rowptr,
            int* L_colind,
            FPP* L_vals,
            unsigned int L_nnz,
            unsigned int K,
            bool on_gpu);
    static FaustCoreCpp<FPP,DEV>* polyBasis_ext(
            unsigned int L_nrows, unsigned int L_ncols,
            int* L_rowptr,
            int* L_colind,
            FPP* L_vals,
            unsigned int L_nnz,
            unsigned int K,
            int* T0_rowptr,
            int* T0_colind,
            FPP* T0_vals,
            unsigned int T0_nnz,
            unsigned int T0_ncols,
            bool on_gpu);


    protected :
    Faust::TransformHelper<FPP,DEV> *transform;

#ifdef USE_GPU_MOD
friend FaustCoreCpp<FPP, Cpu>* clone_gpu2cpu<>(FaustCoreCpp<FPP, GPU2>* gpu_t);
friend FaustCoreCpp<FPP, GPU2>* clone_cpu2gpu<>(FaustCoreCpp<FPP, Cpu>* cpu_t);

#endif

};

template<typename FPP>
void polyCoeffs(int d, int K, int n, const FPP* basisX, const FPP* coeffs, FPP* out, bool on_gpu);
template<typename FPP>
void polyGroupCoeffs_(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP** out, int n_out, bool on_gpu);

template<typename FPP>
using FaustCoreCppCPU = FaustCoreCpp<FPP, Cpu>;

#ifdef USE_GPU_MOD
template<typename FPP>
using FaustCoreCppGPU = FaustCoreCpp<FPP, GPU2>;
#endif

void* _enable_gpu_mod(const char* libpath, const bool silent);
bool _is_gpu_mod_enabled();

template<typename FPP, typename FPP2, FDevice DEV=Cpu>
struct FaustCoreCpp2
{
    static FaustCoreCpp<FPP2,DEV>* real(FaustCoreCpp<FPP,DEV>* core);
};

template<typename FPP, typename FPP2>
using FaustCoreCppCPU2 = FaustCoreCpp2<FPP, FPP2, Cpu>;

#ifdef USE_GPU_MOD
template<typename FPP, typename FPP2>
using FaustCoreCppGPU2 = FaustCoreCpp2<FPP, FPP2, GPU2>;
#endif

#include "FaustCoreCpp.hpp"
#include "FaustCoreCppGPU.hpp"
#endif
