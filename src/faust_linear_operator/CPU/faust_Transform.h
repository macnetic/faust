/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
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
#ifndef __FAUST_Transform_H__
#define __FAUST_Transform_H__

#include <vector>
#include "faust_LinearOperator.h"


#include "faust_transform_algebra.h"
#include "faust_RefManager.h"

//#include "faust_Vect.h"
//#include "faust_MatDense.h"
//#include "faust_MatSparse.h"

//! \class Faust::Transform
//! \brief contains the data of the FAUST method (sparses matrix, lambda, number of non-zeros...) */


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

	template<typename FPP,Device DEVICE> class LinearOperator;

	template<typename FPP,Device DEVICE> class Transform;
	template<typename FPP,Device DEVICE> class Vect;
	template<typename FPP,Device DEVICE> class MatDense;
	template<typename FPP,Device DEVICE> class MatSparse;
	template<Device DEVICE> class BlasHandle;
	template<Device DEVICE> class SpBlasHandle;
	template<typename FPP, Device DEVICE> class TransformHelper;


	// forward definition of friend function
	template<typename FPP>
		Faust::Vect<FPP,Cpu> operator*(const Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu>& v);
	template<typename FPP>
		Faust::MatDense<FPP,Cpu> operator*(const Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu>& M);


	template<typename FPP>
		class Transform<FPP,Cpu> : public Faust::LinearOperator<FPP,Cpu>
		{

			public:
				void faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB)const;

				/** \brief Constructor
				 * \param data : Vector including sparse matrix
				 * \param totalNonZeros : Number of nonzeros Value in the data (all factorized matrix) */
				Transform();

				/** \brief Constructor
				 * \tparam facts : std::Vector including MatGeneric pointer representing the factor of the Transform that will be copied in the Transform
				 * \tparam lambda (optional) : the multiplicative scalar (default value 1.0)
				 * \tparam optimizedCopy (optional) : boolean to control which type of copy of facts is made,
				 if True, the copy is optimized, the dynamic type of the factor can changed
				 if False, the dynamic type stay the same
				 (default value false)*/
				Transform(const std::vector<Faust::MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact=true);

				/** \brief Copy constructor. */
				Transform(const Transform<FPP,Cpu> & A);

				/** \brief Move constructor. Factors are not duplicated in memory but T loses its.*/
				Transform(Transform<FPP,Cpu> && T);

				Transform(const Transform<FPP, Cpu>* A, const bool transpose_A, const bool conj_A, const Transform<FPP, Cpu>* B, const bool transpose_B, const bool conj_B);
				/** \brief
				 * check the factors validity of the faust, if the list of factors represents a valid matrix
				 * */
				void check_factors_validity() const;


				/** \brief Constructor
				 * \param facts : Vector including dense matrix*/
				Transform(const std::vector<Faust::MatDense<FPP,Cpu> >&facts, const bool optimizedCopy=false);
				Transform(const std::vector<Faust::MatSparse<FPP,Cpu> >&facts, const bool optimizedCopy=false);

				//void get_facts(std::vector<Faust::MatSparse<FPP,Cpu> >& sparse_facts)const{sparse_facts = data;}
				//void get_facts(std::vector<Faust::MatDense<FPP,Cpu> >& facts)const;

				faust_unsigned_int size()const{return data.size();}
				void size(int size_)const{ data.resize(size_);}


				//template<Device DEVICE> class BlasHandle;
				//template<Device DEVICE> class SpBlasHandle;
				/** \brief Perform the product of all factorized matrix. */
				Faust::MatDense<FPP,Cpu> get_product(const char opThis='N', const bool isConj=false)const;
				Faust::MatDense<FPP,Cpu> get_product(Faust::BlasHandle<Cpu> blas_handle,Faust::SpBlasHandle<Cpu> spblas_handle)const;


				/** \brief return a copy of the factor of index id
				//  \warning dynamic memory allocation is made for the return pointer if cloning_fact == true*/
				Faust::MatGeneric<FPP,Cpu>* get_fact(faust_unsigned_int id, const bool cloning_fact = true) const;
				bool is_fact_sparse(const faust_unsigned_int id) const;
				bool is_fact_dense(const faust_unsigned_int id) const;
				faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
				void get_fact(const faust_unsigned_int id,
						 const int** row_ids,
						 const int** col_ids,
						 const FPP** elts,
						 faust_unsigned_int* size,
						 faust_unsigned_int* num_rows,
						 faust_unsigned_int* num_cols) const;
				void get_fact(const faust_unsigned_int id,
						int* d_outer_count_ptr, int* d_inner_ptr, FPP* d_elts,
						faust_unsigned_int* nnz,
						faust_unsigned_int* num_rows, faust_unsigned_int* num_cols,
						const bool transpose=false) const;
				void get_fact(const faust_unsigned_int id,
						const FPP** elts,
						faust_unsigned_int* num_rows,
						faust_unsigned_int* num_cols) const;
				void get_fact(const faust_unsigned_int id,
						FPP* elts,
						faust_unsigned_int* num_rows,
						faust_unsigned_int* num_cols,
						const bool transpose=false) const;
				faust_unsigned_int getNbRow() const;
				faust_unsigned_int getNbCol() const;
				void print_file(const char* filename) const;
				void print_data_ptrs() const
				{
					for(int i=0;i<data.size();i++)
						cout << data[i] << " ";
					cout << endl;
				}
				void init_from_file(const char* filename);
				/**
				 *	\brief Writes the FAuST into a Matlab file. The product is written as a cell array with the matrix factors as elements.
				 *	\arg \c filename The filepath to the output file (preferably with a .mat suffix).
				 *	\arg \c transpose The boolean to set to true if you want to save the transpose Faust of this instance or false otherwise.
				 *	\throw std::logic_error if any problem occurs.
				 */ 
				void save_mat_file(const char* filename, bool transpose, bool conjugate=false) const;
				long long int get_total_nnz()const{return totalNonZeros;}
				void clear(){data.resize(0);totalNonZeros=0;}

				/** \brief add M to the end of the Faust F, work as for std::vector::push_back
				  F={S_0,S_1,S_2,..,S_n} becomes {S_0,S_1,S_2,..,S_n,copy(M)}
				 * \tparam M : MatGeneric pointer representing the factor that will be copied at the end of the Faust

				 * \param optimizedCopy (optional) : boolean to control which type of copy of the fact is made,
				 if True, the copy is optimized, the dynamic type of the factor can changed
				 if False, the dynamic type stay the same
				 (default value false)*/
				void push_back(const Faust::MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const bool conjugate=false);



				/** \brief Constructor
				 * \tparam M : MatGeneric pointer representing the factor that will be copied
				 * \param optimizedCopy (optional) : boolean to control which type of copy of the fact is made,
				 if True, the copy is optimized, the dynamic type of the factor can changed
				 if False, the dynamic type stay the same
				 (default value false)*/
				void push_first(const Faust::MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false);
				//void pop_back(Faust::MatGeneric<FPP,Cpu>* M);
				//void pop_first(Faust::MatGeneric<FPP,Cpu>* M);
				//void pop_first(Faust::MatGeneric<FPP,Cpu>* M) const;
				void Display(const bool transpose=false, const bool displaying_small_mat_elts=false)const;
				std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false)const;
				void transpose();
				/**
				 * \brief Replaces the Faust by its conjugate.
				 */
				void conjugate();
				void updateNonZeros();
				void setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const;
				/**
				 * \brief Returns the first and last indices of matrices that are not orthogonal and where between them none or not all matrices are orthogonal.
				 */
				void get_nonortho_interior_prod_ids(int &start_id, int &end_id);
				///(*this) = (*this) * A
				void multiply(const Transform<FPP,Cpu> & A);
				///(*this) = A * (*this)
				void multiplyLeft(const Transform<FPP,Cpu> & A);
				void scalarMultiply(const FPP scalar);
				float getRCG() const{return ((float)(getNbRow()*getNbCol()))/((float) get_total_nnz());}
				double spectralNorm(const int nbr_iter_max, double threshold, int &flag) const;
				double normL1(const bool transpose=false) const;
				double normL1(const char opThis, const bool isConj) const;
				double normFro() const;
				double normFro(const char opThis, const bool isConj) const;

				static void delete_fact(void * fact)
				{
#ifdef DEBUG
					cout << "Faust::Transform delete_fact" << endl;
#endif
					delete static_cast<MatGeneric<FPP,Cpu>*>(fact);
				}

				~Transform(){
#ifdef FAUST_VERBOSE
					std::cout << "~Transform()" << std::endl;
#endif
					for (int i=0;i<data.size();i++)
						if(this->dtor_delete_data)
							delete data[i];
						else
							ref_man.release(data[i]);
				}

				/*!
				 * \brief multiplication between with vector x
				 *  x = op((*this)) * x
				 <br>
				 * op((*this)) = (*this) if opThis='N', op((*this)= = transpose((*this)) if opThis='T'<br>
				 *! \tparam  x :  the vector to be multiplied
				 *! \tparam opThis : character
				 */
				Faust::Vect<FPP,Cpu> multiply(const Faust::Vect<FPP,Cpu> x,const char opThis) const;

				Faust::Vect<FPP,Cpu> multiply(const Faust::Vect<FPP,Cpu>& x) const
				{
					return this->multiply(x,'N');
				}


				/*!
				 * \brief multiplication between with vector x
				 *  x = op((*this)) * x
				 <br>
				 * op((*this)) = (*this) if opThis='N', op((*this)= = transpose((*this)) if opThis='T'<br>
				 *! \tparam  x :  the vector to be multiplied
				 *! \tparam opThis : character
				 */
				Faust::MatDense<FPP,Cpu> multiply(const Faust::MatDense<FPP,Cpu> A,const char opThis='N') const;

				/** \brief Move assign operator overload. Factors are not duplicated in memory but T loses its.*/
				Transform<FPP,Cpu>& operator=(Faust::Transform<FPP,Cpu>&& T);

				void operator=(const Transform<FPP,Cpu>&  f);//{data=f.data;totalNonZeros=f.totalNonZeros;}
				/// add all of the sparse matrices from f.data to this->data
				void operator*=(const FPP  scalar){scalarMultiply(scalar);};
				void operator*=(const Transform<FPP,Cpu>&  f){multiply(f);};


#ifdef __COMPILE_TIMERS__
				void print_timers() const;
#endif

			private:
				void disable_data_deletion() { this->dtor_delete_data = false; }
				void enable_data_deletion() { this->dtor_delete_data = true; }

			private:
				long long int totalNonZeros;
				static const char * m_className;
				std::vector<Faust::MatGeneric<FPP,Cpu>*> data;
				bool dtor_delete_data;
				static RefManager ref_man;

#ifdef __COMPILE_TIMERS__
				mutable std::vector<Faust::Timer> t_multiply_vector;
#endif


				// friend function
				friend Faust::Vect<FPP,Cpu> Faust::operator*<>(const Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu>& v);
				friend Faust::MatDense<FPP,Cpu> Faust::operator*<>(const Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu>& M);
				friend TransformHelper<FPP,Cpu>;
		};

}



#include "faust_Transform.hpp"


#endif
