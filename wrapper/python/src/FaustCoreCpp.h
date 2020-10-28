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
#include "faust_TransformHelper.h"
#ifdef USE_GPU_MOD
#include "faust_TransformHelper_gpu.h"
#endif

template<typename FPP, FDevice DEV=Cpu>
class FaustCoreCpp
{

    public :


    FaustCoreCpp(): transform(nullptr) {}
    FaustCoreCpp(Faust::TransformHelper<FPP,DEV> *th);
    void Display() const { transform->display();}
    const char* to_string() const;
    void push_back(FPP* valueMat,unsigned int nbrow,unsigned int nbcol, bool optimizedCopy=false);
    void push_back(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, bool optimizedCopy=false);
    unsigned int getNBytes() const;
    unsigned int getNbRow() const;
    unsigned int getNbCol() const;
    void get_product(FPP* y_data, int y_nrows, int y_ncols);
    void multiply(FPP* y_data, int y_nrows, int y_ncols, FPP* x_data, int* x_row_ptr, int* x_id_col, int x_nnz, int x_nrows, int x_ncols);
    void set_FM_mul_mode(const int mode);
    void set_Fv_mul_mode(const int mode);
    void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x/*,bool isTranspose*/)const;
    FaustCoreCpp<FPP,DEV>* mul_faust(FaustCoreCpp<FPP,DEV>* right);
    FaustCoreCpp<FPP,DEV>* vertcat(FaustCoreCpp<FPP,DEV>* right);
    FaustCoreCpp<FPP,DEV>* horzcat(FaustCoreCpp<FPP,DEV>* right);
    FaustCoreCpp<FPP,DEV>* mul_scal(FPP scal);
    FaustCoreCpp<FPP,DEV>* normalize(int ord) const;
    unsigned long long nnz()const;
    double norm(int ord, double threshold=.001, int max_num_its=100) const;
    double normFro() const;
    double normInf() const;
    double get_nb_factors() const;
    unsigned int get_fact_nb_rows(unsigned int& i) const;
    unsigned int get_fact_nb_cols(unsigned int& i) const;
    void get_fact(const unsigned int& i, FPP* fact_ptr) const;
    void get_fact_sparse(const unsigned int& i,
            int* rowptr,
            int* col_ids,
            FPP* elts,
            const bool transpose = false) const;
    void get_fact_dense(const unsigned int& i, FPP* elts,
            unsigned int* num_rows, unsigned int* num_cols,
            const bool transpose) const;
    FaustCoreCpp<FPP,DEV>* left(const faust_unsigned_int) const;
    FaustCoreCpp<FPP,DEV>* right(const faust_unsigned_int) const;
    faust_unsigned_int get_fact_nnz(const faust_unsigned_int) const;
    bool is_fact_sparse(const faust_unsigned_int id) const;
    FaustCoreCpp<FPP,DEV>* slice(unsigned int, unsigned int, unsigned int, unsigned int);
    FaustCoreCpp<FPP,DEV>* fancy_idx(unsigned long int* row_ids, unsigned long int
                                  num_rows, unsigned long int* col_ids,
                                  unsigned long int num_cols);
    bool save_mat_file(const char* filepath) const;
    FaustCoreCpp<FPP,DEV>* optimize_storage(const bool time=false);
    FaustCoreCpp<FPP,DEV>* optimize(const bool transp=false);
    FaustCoreCpp<FPP,DEV>* optimize_time(const bool transp=false, const bool inplace=false, const int nsamples=1);
    const bool isTransposed();
    FaustCoreCpp<FPP,DEV>* transpose();
    FaustCoreCpp<FPP,DEV>* conjugate();
    FaustCoreCpp<FPP,DEV>* adjoint();
    FaustCoreCpp<FPP,DEV>* zpruneout(const int nnz_tres, const int npasses, const bool only_forward);
    ~FaustCoreCpp();
    static FaustCoreCpp<FPP,DEV>* randFaust(unsigned int t,
            unsigned int min_num_factors, unsigned int max_num_factors,
            unsigned int min_dim_size,
            unsigned int max_dim_size, float density, bool per_row);
    static FaustCoreCpp<FPP,DEV>* hadamardFaust(unsigned int n, const bool norma);
    static FaustCoreCpp<FPP,DEV>* fourierFaust(unsigned int n, const bool norma);
    static FaustCoreCpp<FPP,DEV>* eyeFaust(unsigned int n, unsigned int m);

    protected :
    Faust::TransformHelper<FPP,DEV> *transform;
};

template<typename FPP>
class FaustCoreCppGPU: public FaustCoreCpp<FPP, GPU2>
{
	public:
		FaustCoreCppGPU() {}
		FaustCoreCppGPU(Faust::TransformHelper<FPP,GPU2> *th);
		void get_product(FPP* y_data, int y_nrows, int y_ncols);
		FaustCoreCppGPU<FPP>* mul_faust_gpu(FaustCoreCppGPU<FPP>* right);
		static FaustCoreCppGPU<FPP>* randFaustGPU(unsigned int t,
				unsigned int min_num_factors, unsigned int max_num_factors,
				unsigned int min_dim_size,
				unsigned int max_dim_size, float density, bool per_row);

		static FaustCoreCppGPU<FPP>* hadamardFaustGPU(unsigned int n, const bool norma);
		static FaustCoreCppGPU<FPP>* fourierFaustGPU(unsigned int n, const bool norma);
		static FaustCoreCppGPU<FPP>* eyeFaustGPU(unsigned int n, unsigned int m);
};

void* _enable_gpu_mod(const char* libpath);

#include "FaustCoreCpp.hpp"

#endif
