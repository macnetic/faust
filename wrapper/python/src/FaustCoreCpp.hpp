/****************************************************************************/
/*                              Description:                                */
/*  C++ Hpp file wrapping the  Class Faust::Transform,                      */
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

#include "faust_Transform.h"
#include <iostream>
#include <exception>


template<typename FPP>
FaustCoreCpp<FPP>::FaustCoreCpp(Faust::TransformHelper<FPP,Cpu> *th)
{
    this->transform = th;
}

    template<typename FPP>
void FaustCoreCpp<FPP>::push_back(FPP* valueMat, unsigned int nbrow,unsigned int nbcol, bool optimizedCopy /* false by deft */)
{
    Faust::MatDense<FPP,Cpu> dense_mat(valueMat,nbrow,nbcol);
    // Faust::MatSparse<FPP,Cpu> sparse_mat(dense_mat);
    // sparse_mat.Display();
    if(transform == nullptr) transform = new Faust::TransformHelper<FPP,Cpu>();
    this->transform->push_back(&dense_mat, optimizedCopy);




}

template<typename FPP>
void FaustCoreCpp<FPP>::push_back(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, bool optimizedCopy /* false by deft */)
{
    //    Faust::MatSparse<FPP,Cpu>::MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP1* value, const int* row_ptr, const int* id_col) :
      Faust::MatSparse<FPP, Cpu> sparse_mat(nnz, nrows, ncols, data, row_ptr, id_col);
//      std::cout << "FaustCoreCpp::push_back()" << nrows << " " << sparse_mat.getNbRow() <<  " " << ncols << " " << sparse_mat.getNbCol() << std::endl;
      if(transform == nullptr) transform = new Faust::TransformHelper<FPP,Cpu>();
      this->transform->push_back(&sparse_mat, optimizedCopy);
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::mul_faust(FaustCoreCpp<FPP>* right)
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->multiply(right->transform);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::vertcat(FaustCoreCpp<FPP>* right)
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->vertcat(right->transform);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::horzcat(FaustCoreCpp<FPP>* right)
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->horzcat(right->transform);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::mul_scal(FPP scal)
 {
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->multiply(scal);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}
template<typename FPP>
void FaustCoreCpp<FPP>::get_product(FPP* y_data, int y_nrows, int y_ncols)
{
    Faust::MatDense<FPP, Cpu> Y = this->transform->get_product();
    memcpy(y_data, Y.getData(), sizeof(FPP)*y_ncols*y_nrows);
}

template<typename FPP>
void FaustCoreCpp<FPP>::multiply(FPP* y_data, int y_nrows, int y_ncols, FPP* x_data, int* x_row_ptr, int* x_id_col, int x_nnz, int x_nrows, int x_ncols)
{
    Faust::MatSparse<FPP, Cpu> X(x_nnz, x_nrows, x_ncols, x_data, x_row_ptr, x_id_col);
    Faust::MatDense<FPP, Cpu> Y;
    Y = this->transform->multiply(X);
    memcpy(y_data, Y.getData(),sizeof(FPP)*y_ncols*y_nrows);
}

template<typename FPP>
void FaustCoreCpp<FPP>::multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x)const
{



    unsigned int nbRowThis,nbColThis;


    nbRowThis = this->getNbRow();
    nbColThis = this->getNbCol();

    if ( (nbrow_y != nbRowThis) | (nbrow_x != nbColThis) | (nbcol_y != nbcol_x) )
    {
        std::cout<<"nbRowThis "<<nbRowThis<<" must be equal to nb row y  "<<nbrow_y<<std::endl;
        std::cout<<"nbColThis "<<nbColThis<<" must be equal to nb row x  "<<nbrow_x<<std::endl;
        std::cout<<"nbcol_y "<<nbcol_y<<" must be equal to nbcol_x  "<<nbcol_x<<std::endl;
        handleError("FaustCpp"," multiply : invalid dimension");
    }
    if (nbcol_x == 1)
    {
        Faust::Vect<FPP,Cpu> X(nbrow_x,value_x);
        Faust::Vect<FPP,Cpu> Y;


        Y = this->transform->multiply(X);

        memcpy(value_y,Y.getData(),sizeof(FPP)*nbrow_y);
    }else
    {
        Faust::MatDense<FPP,Cpu> X(value_x,nbrow_x,nbcol_x);
        Faust::MatDense<FPP,Cpu> Y;

        Y = this->transform->multiply(X);

        memcpy(value_y,Y.getData(),sizeof(FPP)*nbrow_y*nbcol_y);
    }
}

template<typename FPP>
void FaustCoreCpp<FPP>::set_mul_order_opt_mode(const int mode)
{
    this->transform->set_mul_order_opt_mode(mode);
}

template<typename FPP>
unsigned long long FaustCoreCpp<FPP>::nnz() const
{
    return this->transform->get_total_nnz();
}


template<typename FPP>
double FaustCoreCpp<FPP>::norm(int ord, double threshold, int max_num_its) const
{
    int flag; //not used yet
    switch(ord) {
        case 2:
            return this->transform->spectralNorm(max_num_its, threshold, flag);
        case 1:
            return this->transform->normL1();
        default:
            handleError("FaustCoreCpp", "norm(int ord) unvalid norm order.");
    }
    return -1;
}

template<typename FPP>
double FaustCoreCpp<FPP>::normFro() const
{
    return this->transform->normFro();
}

template<typename FPP>
double FaustCoreCpp<FPP>::normInf() const
{
    return this->transform->normInf();
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::normalize(int ord) const
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->normalize(ord);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
double FaustCoreCpp<FPP>::get_nb_factors() const
{
    double nb_fact = this->transform->size();
    return nb_fact;
}

template<typename FPP>
const char* FaustCoreCpp<FPP>::to_string() const
{
    std::string str = this->transform->to_string();
    char * c_str = (char*) malloc(str.size()+1);
    memcpy(c_str, str.c_str(), str.size()+1);
    return (const char*) c_str;
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::get_fact_nb_rows(unsigned int& i) const
{
    return this->transform->get_fact_nb_rows(i);
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::get_fact_nb_cols(unsigned int& i) const
{
    return this->transform->get_fact_nb_cols(i);
}

template<typename FPP>
faust_unsigned_int FaustCoreCpp<FPP>::get_fact_nnz(const faust_unsigned_int i) const
{
    return transform->get_fact_nnz(i);
}

template<typename FPP>
void FaustCoreCpp<FPP>::get_fact(const unsigned int& i, FPP* fact_ptr) const
{
    Faust::MatDense<FPP,Cpu> dense_factor = this->transform->get_fact(i);
    // not optimized here (we have three copies from C++ object to Py, the first in MatDense::Clone()
    // (called from Faust::Transform::get_fact()) the second when converting to
    // MatDense (even if not a MatSparse, the copy is made)
    // and finally a third copy here)
    memcpy(fact_ptr, dense_factor.getData(),
            sizeof(FPP)*dense_factor.getNbCol()*dense_factor.getNbRow());
    //optimized versions are get_fact_dense(), get_fact_sparse()
}

template<typename FPP>
void FaustCoreCpp<FPP>::get_fact_dense(const unsigned int& i, FPP* fact_ptr,
        unsigned int* num_rows, unsigned int* num_cols, const bool transpose) const
{
    faust_unsigned_int nrows, ncols;
    // only one copy, directly in numpy buffer (optimization)
    this->transform->get_fact(i, fact_ptr, &nrows, &ncols, transpose);
    if(num_rows != NULL)
        *num_rows = transpose?nrows:ncols;
    if(num_cols != NULL)
        *num_cols = transpose?ncols:nrows;
}

template<typename FPP>
void FaustCoreCpp<FPP>::get_fact_sparse(const unsigned int& id,
        int* rowptr,
        int* col_ids,
        FPP* elts,
        const bool transpose) const
{
    faust_unsigned_int size, nrows, ncols; //TODO: delete nrows,ncols when NULL arg's ok
    // only one copy per buffer (optimization) directly from underlying
    // MatSparse buffers to scipy buffers
    transform->get_fact(id, rowptr, col_ids, elts, &size, &nrows, &ncols, transpose);
}

template<typename FPP>
bool FaustCoreCpp<FPP>::is_fact_sparse(const faust_unsigned_int id) const
{
    return transform->is_fact_sparse(id);
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::slice(unsigned int start_row_id, unsigned int end_row_id, unsigned int start_col_id, unsigned int end_col_id)
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->slice(start_row_id, end_row_id, start_col_id, end_col_id);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::right(const faust_unsigned_int id) const
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->right(id);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::left(const faust_unsigned_int id) const
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->left(id);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::fancy_idx(unsigned long int* row_ids, unsigned long int
        num_rows, unsigned long int* col_ids,
        unsigned long int num_cols)
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->fancy_index(row_ids, num_rows, col_ids, num_cols);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
bool FaustCoreCpp<FPP>::save_mat_file(const char* filepath) const
{
    //    std::cout << "FaustCoreCpp::save_mat_file()" << std::endl;
    try
    {
        this->transform->save_mat_file(filepath);
        return true;
    }
    catch(exception& e) {
        //cerr << e.what() << endl;
        return false;
    }
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::getNbRow() const
{
    return this->transform->getNbRow();
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::getNBytes() const
{
    return this->transform->getNBytes();
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::getNbCol() const
{
    return this->transform->getNbCol();
}

    template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::transpose()
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->transpose();
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::transpose() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP>
void FaustCoreCpp<FPP>::optimize_mul(const bool transp /* deft to false*/)
{
    this->transform->optimize_multiply(transp);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::optimize_mul() th=" << th << "core=" << core << std::endl;
#endif
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::optimize(const bool transp /* deft to false*/)
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->optimize(transp);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::optimize() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

    template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::optimize_storage(const bool time)
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->optimize_storage(time);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::optimize_storage() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP>
const bool FaustCoreCpp<FPP>::isTransposed()
{
    return transform->isTransposed();
}

    template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::conjugate()
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->conjugate();
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::conjugate() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP>
  FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::randFaust(unsigned int t,
            unsigned int min_num_factors, unsigned int max_num_factors,
            unsigned int min_dim_size, unsigned int max_dim_size, float density, bool per_row) {
      Faust::TransformHelper<FPP,Cpu>* th = Faust::TransformHelper<FPP,Cpu>::randFaust(Faust::RandFaustType(t), min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
      if(!th) return NULL;
      FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
      return core;
  }

//TODO: hadamardFaust and fourierFaust shouldn't be here... or at least it
//should be reconsidered
template<typename FPP>
  FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::hadamardFaust(unsigned int n, const bool norma)
{
      Faust::TransformHelper<FPP,Cpu>* th = Faust::TransformHelper<FPP,Cpu>::hadamardFaust(n, norma);
      if(!th) return NULL;
      FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
      return core;
}

template<typename FPP>
  FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::fourierFaust(unsigned int n, const bool norma)
{
      Faust::TransformHelper<FPP,Cpu>* th = Faust::TransformHelper<FPP,Cpu>::fourierFaust(n, norma);
      if(!th) return NULL;
      FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
      return core;
}

template<typename FPP>
  FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::eyeFaust(unsigned int n, unsigned int m)
{
      Faust::TransformHelper<FPP,Cpu>* th = Faust::TransformHelper<FPP,Cpu>::eyeFaust(n, m);
      if(!th) return NULL;
      FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
      return core;
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::adjoint()
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->adjoint();
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::adjoint() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::zpruneout(const int nnz_tres, const int npasses, const bool only_forward)
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->pruneout(nnz_tres, npasses, only_forward);
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP>::~FaustCoreCpp()
{
    if(transform) delete transform;
}

