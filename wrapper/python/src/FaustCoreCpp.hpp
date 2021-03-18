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
#include <vector>


template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP, DEV>::FaustCoreCpp(Faust::TransformHelper<FPP,DEV> *th)
{
    this->transform = th;
}

    template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::push_back(FPP* valueMat, unsigned int nbrow,unsigned int nbcol, bool optimizedCopy /* false by deft */)
{
    Faust::MatDense<FPP,DEV> dense_mat(nbrow, nbcol, valueMat);
    // Faust::MatSparse<FPP,DEV> sparse_mat(dense_mat);
    // sparse_mat.Display();
    if(transform == nullptr) transform = new Faust::TransformHelper<FPP,DEV>();
    this->transform->push_back(&dense_mat, optimizedCopy);
}

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::push_back(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, bool optimizedCopy /* false by deft */)
{
    //    Faust::MatSparse<FPP,DEV>::MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP1* value, const int* row_ptr, const int* id_col) :
//      Faust::MatSparse<FPP, DEV> sparse_mat(nnz, nrows, ncols, data, row_ptr, id_col);
//      std::cout << "FaustCoreCpp::push_back()" << nrows << " " << sparse_mat.getNbRow() <<  " " << ncols << " " << sparse_mat.getNbCol() << std::endl;
      if(transform == nullptr) transform = new Faust::TransformHelper<FPP,DEV>();
//      this->transform->push_back(&sparse_mat, optimizedCopy);
      this->transform->push_back(data, row_ptr, id_col, nnz, nrows, ncols, optimizedCopy);
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::mul_faust(FaustCoreCpp<FPP,DEV>* right)
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->multiply(right->transform);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::vertcat(FaustCoreCpp<FPP,DEV>* right) const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->vertcat(right->transform);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::vertcatn(FaustCoreCpp<FPP,DEV>** rights, size_t n) const
{
//    Faust::TransformHelper<FPP,DEV>* th = this->transform->vertcat(right->transform);
    Faust::TransformHelper<FPP,DEV>* th = nullptr;
    std::vector<Faust::TransformHelper<FPP,DEV>*> ths;
    ths.push_back(this->transform);
    for(int i =0; i < n; i++)
        ths.push_back(rights[i]->transform);
    th = Faust::vertcat(ths);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::horzcatn(FaustCoreCpp<FPP,DEV>** rights, size_t n) const
{
//    Faust::TransformHelper<FPP,DEV>* th = this->transform->vertcat(right->transform);
    Faust::TransformHelper<FPP,DEV>* th = nullptr;
    std::vector<Faust::TransformHelper<FPP,DEV>*> ths;
    ths.push_back(this->transform);
    for(int i =0; i < n; i++)
        ths.push_back(rights[i]->transform);
    th = Faust::horzcat(ths);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::horzcat(FaustCoreCpp<FPP,DEV>* right) const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->horzcat(right->transform);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::mul_scal(FPP scal)
 {
    Faust::TransformHelper<FPP,DEV>* th = this->transform->multiply(scal);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::get_product(FPP* y_data, int y_nrows, int y_ncols)
{
    Faust::MatDense<FPP, DEV> Y = this->transform->get_product();
    memcpy(y_data, Y.getData(), sizeof(FPP)*y_ncols*y_nrows);
}

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::multiply(FPP* y_data, int y_nrows, int y_ncols, FPP* x_data, int* x_row_ptr, int* x_id_col, int x_nnz, int x_nrows, int x_ncols)
{
    Faust::MatSparse<FPP, DEV> X(x_nnz, x_nrows, x_ncols, x_data, x_row_ptr, x_id_col);
    Faust::MatDense<FPP, DEV> Y;
    Y = this->transform->multiply(X);
    memcpy(y_data, Y.getData(),sizeof(FPP)*y_ncols*y_nrows);
}

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::multiply(FPP* value_y,int nbrow_y,int nbcol_y, const FPP* value_x,int nbrow_x,int nbcol_x)const
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
//        Faust::Vect<FPP,Cpu> X(nbrow_x,value_x);
//        Faust::Vect<FPP,Cpu> Y;

//        Y = this->transform->multiply(X);
//        Y = this->transform->multiply(value_x);
//        memcpy(value_y,Y.getData(),sizeof(FPP)*nbrow_y);
        // assuming that x and y are pointers to memory allocated to the proper
        // sizes
        this->transform->multiply(value_x, value_y);
    }else
    {
        Faust::MatDense<FPP,Cpu> X(value_x,nbrow_x,nbcol_x);
        Faust::MatDense<FPP,Cpu> Y;

        Y = this->transform->multiply(X);

        memcpy(value_y,Y.getData(),sizeof(FPP)*nbrow_y*nbcol_y);
    }
}
template<typename FPP, FDevice DEV>
 FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::polyNext() const
{
    Faust::TransformHelperPoly<FPP> *transform_poly = dynamic_cast<Faust::TransformHelperPoly<FPP>*>(this->transform);
    if(nullptr)
        throw std::runtime_error("polyCoeffs can only be used on a Poly. specialized Faust.");
    auto th = transform_poly->next();
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::polyCoeffs(int d, int n, const FPP* basisX, const FPP* coeffs, FPP* out) const
{
    Faust::TransformHelperPoly<FPP> *transform_poly = dynamic_cast<Faust::TransformHelperPoly<FPP>*>(this->transform);
    if(nullptr)
        throw std::runtime_error("polyCoeffs can only be used on a Poly. specialized Faust.");
    transform_poly->poly(d, n, basisX, coeffs, out);
}

template<typename FPP>
void polyCoeffs(int d, int K, int n, const FPP* basisX, const FPP* coeffs, FPP* out)
{
    Faust::poly(d, K, n, basisX, coeffs, out);
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::polyCoeffs(const FPP* coeffs)
{
    //verify if transform is a TransformHelperPoly otherwise raise an exception
    Faust::TransformHelperPoly<FPP> *transform_poly = dynamic_cast<Faust::TransformHelperPoly<FPP>*>(this->transform);
    if(nullptr)
        throw std::runtime_error("polyCoeffs can only be used on a Poly. specialized Faust.");
    Faust::TransformHelper<FPP,DEV>* th = transform_poly->polyFaust(coeffs);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}


template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::set_FM_mul_mode(const int mode)
{
    this->transform->set_FM_mul_mode(mode);
}

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::set_Fv_mul_mode(const int mode)
{
    this->transform->set_Fv_mul_mode(mode);
}

template<typename FPP, FDevice DEV>
unsigned long long FaustCoreCpp<FPP,DEV>::nnz() const
{
    return this->transform->get_total_nnz();
}


template<typename FPP, FDevice DEV>
double FaustCoreCpp<FPP,DEV>::norm(int ord, double threshold, int max_num_its) const
{
    int flag; //not used yet
    switch(ord)
	{
        case 2:
            return this->transform->spectralNorm(max_num_its, threshold, flag);
        case 1:
            return this->transform->normL1();
        default:
            handleError("FaustCoreCpp", "norm(int ord) unvalid norm order.");
    }
	//TODO: call normFro and normInf from here (and delete the two functions below)
    return -1;
}

template<typename FPP, FDevice DEV>
double FaustCoreCpp<FPP,DEV>::normFro() const
{
    return this->transform->normFro();
}

template<typename FPP, FDevice DEV>
double FaustCoreCpp<FPP,DEV>::normInf() const
{
    return this->transform->normInf();
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::normalize(int ord) const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->normalize(ord);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
double FaustCoreCpp<FPP,DEV>::get_nb_factors() const
{
    double nb_fact = this->transform->size();
    return nb_fact;
}

template<typename FPP, FDevice DEV>
const char* FaustCoreCpp<FPP,DEV>::to_string() const
{
	//TODO: ideally a pre-allocated str should be passed by the caller (no allocation inside this function)
	//TODO: until that the caller is responsible to free the returned memory buffer
    std::string str = this->transform->to_string();
    char * c_str = (char*) malloc(str.size()+1);
    memcpy(c_str, str.c_str(), str.size()+1);
    return (const char*) c_str;
}

template<typename FPP, FDevice DEV>
unsigned int FaustCoreCpp<FPP,DEV>::get_fact_nb_rows(unsigned int& i) const
{
    return this->transform->get_fact_nb_rows(i);
}

template<typename FPP, FDevice DEV>
unsigned int FaustCoreCpp<FPP,DEV>::get_fact_nb_cols(unsigned int& i) const
{
    return this->transform->get_fact_nb_cols(i);
}

template<typename FPP, FDevice DEV>
faust_unsigned_int FaustCoreCpp<FPP,DEV>::get_fact_nnz(const faust_unsigned_int i) const
{
    return transform->get_fact_nnz(i);
}

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::get_fact(const unsigned int& i, FPP* fact_ptr) const
{
    Faust::MatDense<FPP,DEV> dense_factor = this->transform->get_fact(i);
    // not optimized here (we have three copies from C++ object to Py, the first in MatDense::Clone()
    // (called from Faust::Transform::get_fact()) the second when converting to
    // MatDense (even if not a MatSparse, the copy is made)
    // and finally a third copy here)
    memcpy(fact_ptr, dense_factor.getData(),
            sizeof(FPP)*dense_factor.getNbCol()*dense_factor.getNbRow());
    //optimized versions are get_fact_dense(), get_fact_sparse()
}

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::get_fact_dense(const unsigned int& i, FPP* fact_ptr,
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

template<typename FPP, FDevice DEV>
void FaustCoreCpp<FPP,DEV>::get_fact_sparse(const unsigned int& id,
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

template<typename FPP, FDevice DEV>
bool FaustCoreCpp<FPP,DEV>::is_fact_sparse(const faust_unsigned_int id) const
{
    return transform->is_fact_sparse(id);
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::slice(unsigned int start_row_id, unsigned int end_row_id, unsigned int start_col_id, unsigned int end_col_id) const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->slice(start_row_id, end_row_id, start_col_id, end_col_id);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::right(const faust_unsigned_int id) const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->right(id);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::left(const faust_unsigned_int id) const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->left(id);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::fancy_idx(unsigned long int* row_ids, unsigned long int
        num_rows, unsigned long int* col_ids,
        unsigned long int num_cols) const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->fancy_index(row_ids, num_rows, col_ids, num_cols);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
bool FaustCoreCpp<FPP,DEV>::save_mat_file(const char* filepath) const
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

template<typename FPP, FDevice DEV>
unsigned int FaustCoreCpp<FPP,DEV>::getNbRow() const
{
    return this->transform->getNbRow();
}

template<typename FPP, FDevice DEV>
unsigned int FaustCoreCpp<FPP,DEV>::getNBytes() const
{
    return this->transform->getNBytes();
}

template<typename FPP, FDevice DEV>
unsigned int FaustCoreCpp<FPP,DEV>::getNbCol() const
{
    return this->transform->getNbCol();
}

    template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::transpose() const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->transpose();
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::transpose() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::optimize_time(const bool transp /* deft to false*/, const bool inplace /* default to false */, const int nsamples /* default to 1*/)
{
    if(inplace)
        this->transform->optimize_time(transp, inplace, nsamples);
    else
    {
        auto th = this->transform->optimize_time(transp, inplace, nsamples);
        return new FaustCoreCpp<FPP,DEV>(th);
    }
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::optimize_time() th=" << th << "core=" << core << std::endl;
#endif
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::optimize(const bool transp /* deft to false*/)
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->optimize(transp);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::optimize() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

    template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::optimize_storage(const bool time)
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->optimize_storage(time);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::optimize_storage() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP, FDevice DEV>
const bool FaustCoreCpp<FPP,DEV>::isTransposed()
{
    return transform->isTransposed();
}

    template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::conjugate() const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->conjugate();
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::conjugate() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP, FDevice DEV>
  FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::randFaust(
          int faust_nrows,
          int faust_ncols,
          unsigned int t,
          unsigned int min_num_factors, unsigned int max_num_factors,
          unsigned int min_dim_size, unsigned int max_dim_size, float density, bool per_row) {
      Faust::TransformHelper<FPP,DEV>* th = Faust::TransformHelper<FPP,DEV>::randFaust(faust_nrows, faust_ncols, Faust::RandFaustType(t), min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
      if(!th) return NULL;
      FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
      return core;
  }

template<typename FPP, FDevice DEV>
  FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::randFaust(
          unsigned int t,
          unsigned int min_num_factors, unsigned int max_num_factors,
          unsigned int min_dim_size, unsigned int max_dim_size, float density, bool per_row) {
      Faust::TransformHelper<FPP,DEV>* th = Faust::TransformHelper<FPP,DEV>::randFaust(Faust::RandFaustType(t), min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
      if(!th) return NULL;
      FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
      return core;
  }

//TODO: hadamardFaust and fourierFaust shouldn't be here... or at least it
//should be reconsidered
template<typename FPP, FDevice DEV>
  FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::hadamardFaust(unsigned int n, const bool norma)
{
      Faust::TransformHelper<FPP,DEV>* th = Faust::TransformHelper<FPP,DEV>::hadamardFaust(n, norma);
      if(!th) return NULL;
      FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
      return core;
}

template<typename FPP, FDevice DEV>
  FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::fourierFaust(unsigned int n, const bool norma)
{
      Faust::TransformHelper<FPP,DEV>* th = Faust::TransformHelper<FPP,DEV>::fourierFaust(n, norma);
      if(!th) return NULL;
      FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
      return core;
}

template<typename FPP, FDevice DEV>
  FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::eyeFaust(unsigned int n, unsigned int m)
{
      Faust::TransformHelper<FPP,DEV>* th = Faust::TransformHelper<FPP,DEV>::eyeFaust(n, m);
      if(!th) return NULL;
      FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
      return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::adjoint() const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->adjoint();
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::adjoint() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::zpruneout(const int nnz_tres, const int npasses, const bool only_forward)
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->pruneout(nnz_tres, npasses, only_forward);
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
void  FaustCoreCpp<FPP,DEV>::device(char* dev) const
{
	std::strncpy(dev, "cpu", 4);
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::swap_cols(const unsigned int id1, const unsigned int id2,
            const bool permutation, const bool inplace)
{
    if(inplace)
    {
      this->transform->swap_cols(id1, id2, permutation, inplace);
      return this;
    }
    else
    {
        auto th = this->transform->swap_cols(id1, id2, permutation, inplace);
        return new FaustCoreCpp<FPP,DEV>(th);
    }
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::swap_rows(const unsigned int id1, const unsigned int id2,
            const bool permutation, const bool inplace)
{
    if(inplace)
    {
      this->transform->swap_rows(id1, id2, permutation, inplace);
      return this;
    }
    else
    {
        auto th = this->transform->swap_rows(id1, id2, permutation, inplace);
        return new FaustCoreCpp<FPP,DEV>(th);
    }
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP, DEV>::clone() const
{
    Faust::TransformHelper<FPP,DEV>* th = this->transform->clone();
    auto core = new FaustCoreCpp<FPP, DEV>(th);
    return core;
}

template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>* FaustCoreCpp<FPP,DEV>::polyBasis(
        unsigned int L_nrows, unsigned int L_ncols,
        int* L_rowptr,
        int* L_colind,
        FPP* L_vals,
        unsigned int L_nnz,
        unsigned int K)
{
    Faust::MatSparse<FPP, DEV> L(L_nnz, L_nrows, L_ncols, L_vals, L_rowptr, L_colind);
    Faust::TransformHelper<FPP,DEV>* th = Faust::basisChebyshev(&L, K);
    if(!th) return NULL;
    FaustCoreCpp<FPP,DEV>* core = new FaustCoreCpp<FPP,DEV>(th);
    return core;
}


template<typename FPP, FDevice DEV>
FaustCoreCpp<FPP,DEV>::~FaustCoreCpp()
{
    if(transform) delete transform;
}

void* _enable_gpu_mod(const char* libpath, const bool silent)
{
#ifdef USE_GPU_MOD
    return Faust::enable_gpu_mod(libpath, silent);
#else
    return nullptr;
#endif
}

bool _is_gpu_mod_enabled()
{
#ifdef USE_GPU_MOD
    return Faust::is_gpu_mod_enabled();
#else
    return false;
#endif
}
