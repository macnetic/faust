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

#include "faust_FFT.h"
#include "faust_WHT.h"

namespace Faust {

	//TODO: move these functions to Faust::utils ?
	template<typename FPP>
		void conjugate(complex<FPP>* elts, faust_unsigned_int n)
		{
			for(faust_unsigned_int i=0; i< n; i++)
				elts[i] = complex<FPP>(elts[i].real(), - elts[i].imag());
		}

	template<typename FPP>
		void conjugate(FPP* elts, faust_unsigned_int n)
		{
		//nothing to do for real numbers
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts,
				const FPP lambda_, const bool optimizedCopy, const bool cloning_fact) : is_transposed(false),
																						is_conjugate(false),
																						is_sliced(false),
																						is_fancy_indexed(false)
	{
		if(lambda_ != FPP(1.0))
			cerr << "WARNING: the constructor argument for multiplying the Faust by a scalar is DEPRECATED and might not be supported in next versions of FAµST." << endl;
		this->transform = make_shared<Transform<FPP,Cpu>>(facts, lambda_, optimizedCopy, cloning_fact);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper() : is_transposed(false), is_conjugate(false), is_sliced(false), is_fancy_indexed(false)
	{
		this->transform = make_shared<Transform<FPP,Cpu>>();
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(Transform<FPP,Cpu> &t, const bool moving /* default to false */) : is_transposed(false), is_conjugate(false), is_sliced(false), is_fancy_indexed(false)
	{
		if(moving)
			this->transform = make_shared<Transform<FPP,Cpu>>(std::move(t));
		else
			this->transform = make_shared<Transform<FPP,Cpu>>(t);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th_left, TransformHelper<FPP,Cpu>* th_right)
		: is_transposed(false), is_conjugate(false), is_sliced(false), is_fancy_indexed(false)
		{
			this->transform = make_shared<Transform<FPP,Cpu>>(th_left->transform.get(), th_left->is_transposed, th_left->is_conjugate,
					th_right->transform.get(), th_right->is_transposed, th_right->is_conjugate);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate)
		{
			this->transform = th->transform;
			this->is_transposed = transpose?!th->is_transposed:th->is_transposed;
			this->is_sliced = th->is_sliced;
			this->is_fancy_indexed = false;
			if(th->is_sliced)
				copy_slices(th);
			this->is_conjugate = conjugate?!th->is_conjugate:th->is_conjugate;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th)
		{
			this->transform = th->transform;
			this->is_transposed = th->is_transposed;
			this->is_conjugate = th->is_conjugate;
			this->is_sliced = th->is_sliced;
			if(th->is_sliced)
				copy_slices(th);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th, Slice s[2])
		{
			this->transform = th->transform; //do not remove this line, necessary for eval_sliced_transform()
			this->is_transposed = th->is_transposed;
			this->is_conjugate = th->is_conjugate;
			if(! (s[0].belong_to(0, th->getNbRow()) || s[1].belong_to(0, th->getNbCol())))
				handleError("Faust::TransformHelper::TransformHelper(TransformHelper,Slice)", "Slice overflows a Faust dimension.");
			this->slices[0] = s[0];
			this->slices[1] = s[1];
			this->is_sliced = true;
			eval_sliced_Transform();
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols)
		{
			this->transform = th->transform; //do not remove this line, necessary for eval*()
			this->is_transposed = th->is_transposed;
			this->is_conjugate = th->is_conjugate;
			this->is_sliced = false;
			//TODO: check indices
//				handleError("Faust::TransformHelper::TransformHelper(TransformHelper,Slice)", "Fancy indexing overflows a Faust dimension.");
			unsigned int id0=0, id1=1;
			this->fancy_num_cols = num_cols;
			this->fancy_num_rows = num_rows;
			if(this->is_transposed)
			{
				id0 = 1;
				id1 = 0;
				this->fancy_num_cols = num_rows;
				this->fancy_num_rows = num_cols;
			}
			this->fancy_indices[id0] = new faust_unsigned_int[num_rows];
			this->fancy_indices[id1] = new faust_unsigned_int[num_cols];
			this->is_fancy_indexed= true;
			memcpy(this->fancy_indices[id0], row_ids, num_rows*sizeof(faust_unsigned_int));
			memcpy(this->fancy_indices[id1], col_ids, num_cols*sizeof(faust_unsigned_int));
			eval_fancy_idx_Transform();
			delete[] this->fancy_indices[0];
			delete[] this->fancy_indices[1];
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> A) const
		{
			MatDense<FPP,Cpu> M = this->transform->multiply(A, isTransposed2char());
			if(is_conjugate) M.conjugate();
			return M;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> x) const
		{
			Vect<FPP,Cpu> v = this->transform->multiply(x, isTransposed2char());
			if(is_conjugate) v.conjugate();
			return v;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> x, const bool transpose)
		{
			is_transposed ^= transpose;
			Vect<FPP,Cpu> v = this->transform->multiply(x, isTransposed2char());
			is_transposed ^= transpose;
			return v;
		}


	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> A, const bool transpose)
		{
			is_transposed ^= transpose;
			MatDense<FPP,Cpu> M = this->transform->multiply(A, isTransposed2char());
			is_transposed ^= transpose;
			return M;
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(TransformHelper<FPP, Cpu>* th_right)
		{
			return new TransformHelper<FPP,Cpu>(this, th_right);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(FPP& scalar)
		{
			const vector<MatGeneric<FPP,Cpu>*>& vec = transform->data; //TransformHelper is a friend class of Transform // we can access private attribute data
			//the point here is to minimize the number of copies (with direct access)
			// the constructor then will copy the factors from the vector
			Transform<FPP,Cpu>* t = new Transform<FPP,Cpu>(vec, scalar, false, true); //optimizedCopy == false, cloning_fact == true
			TransformHelper<FPP,Cpu>* th  = new TransformHelper<FPP,Cpu>(*t);
			//TODO: refactor the attributes copy ? 
			th->is_transposed = is_transposed;
			th->is_conjugate = is_conjugate;
			th->is_sliced = is_sliced;
			if(is_sliced)
			{
				th->slices[0] = Slice(slices[0]);
				th->slices[1] = Slice(slices[1]);
			}
			return th;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const MatGeneric<FPP,Cpu>* M, bool optimizedCopy /* false by default */)
		{
			//warning: should not be called after initialization of factors (to respect the immutable property)
			//this function is here only for python wrapper (TODO: see how to modify that wrapper in order to delete this function after)
			this->transform->push_back(M, optimizedCopy, is_conjugate); //2nd argument is for opt. (possibly converting dense <-> sparse)
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::getNbRow() const
		{
			if(is_sliced){
				faust_unsigned_int id = is_transposed?1:0;
				return	this->slices[id].end_id-this->slices[id].start_id;
			}
			else
				return is_transposed?transform->getNbCol():transform->getNbRow();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::getNbCol() const
		{
			if(is_sliced)
			{
				faust_unsigned_int id = is_transposed?0:1;
				return this->slices[id].end_id-this->slices[id].start_id;
			}
			else
				return is_transposed?transform->getNbRow():transform->getNbCol();
		}

	template<typename FPP>
		bool TransformHelper<FPP,Cpu>::isReal() const
		{
			return this->transform->isReal();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::get_total_nnz() const
		{
			return this->transform->get_total_nnz();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::size() const
		{
			return this->transform->size();
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::display() const
		{
			this->transform->Display(is_transposed, false);
		}

	template<typename FPP>
		string TransformHelper<FPP,Cpu>::to_string() const
		{
			return this->transform->to_string(is_transposed);
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::get_fact_nnz(const faust_unsigned_int id) const
		{
			if(id == 0 || id == this->size()-1)
				return this->transform->get_fact_nnz(is_transposed?size()-id-1:id);
			else
				return this->transform->get_fact_nnz(is_transposed?size()-id-1:id);
		}

	template<typename FPP>
		unsigned int TransformHelper<FPP,Cpu>::get_fact_nb_rows(const faust_unsigned_int id) const
		{
			return get_fact_dim_size(id,0);
		}

	template<typename FPP>
		unsigned int TransformHelper<FPP,Cpu>::get_fact_nb_cols(const faust_unsigned_int id) const
		{
			return get_fact_dim_size(id,1);
		}

	template<typename FPP>
		unsigned int TransformHelper<FPP,Cpu>::get_fact_dim_size(const faust_unsigned_int id, unsigned short dim) const
		{
			//dim == 0 to get num of cols, >0 otherwise
			faust_unsigned_int rid; //real id
			if(is_transposed) {
				rid = size()-id-1;
				dim = !dim;
			}
			else {
				rid = id;
				//dim = dim;
			}
			Faust::MatGeneric<FPP,Cpu>* mat;
			if(rid == 0 || rid == this->size()-1)
				mat = this->transform->get_fact(rid, false);
			else
				mat = this->transform->get_fact(rid, false);
			if(dim)
				return mat->getNbCol();
			else
				return mat->getNbRow();
		}

	//private
	template<typename FPP>
		const MatGeneric<FPP,Cpu>* TransformHelper<FPP,Cpu>::get_gen_fact(const faust_unsigned_int id) const
		{
			return this->transform->data[is_transposed?size()-id-1:id];
		}

	template<typename FPP>
		unsigned long long TransformHelper<FPP,Cpu>::get_fact_addr(const faust_unsigned_int id) const
		{
			return (unsigned long long) this->transform->data[is_transposed?size()-id-1:id];
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
				const int** rowptr,
				const int** col_ids,
				const FPP** elts,
				faust_unsigned_int* nnz,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols) const
		{
			this->transform->get_fact(is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols);

		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
				int* rowptr,
				int* col_ids,
				FPP* elts,
				faust_unsigned_int* nnz,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose /* = false*/) const
		{
			this->transform->get_fact(is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols, is_transposed ^ transpose);
			if(is_conjugate)
				Faust::conjugate(elts, *nnz);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
				const FPP** elts,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols) const
		{
			this->transform->get_fact(is_transposed?this->size()-id-1:id, elts, num_rows, num_cols);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
				FPP* elts,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose /* default to false */) const
		{
			this->transform->get_fact(is_transposed?this->size()-id-1:id, elts, num_rows, num_cols, is_transposed ^ transpose);
			if(is_conjugate)
				Faust::conjugate(elts,*num_cols*(*num_rows));
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_fact(faust_unsigned_int id) const
		{
			MatDense<FPP,Cpu> dense_factor;
			MatGeneric<FPP,Cpu>* factor_generic;
			factor_generic = this->transform->get_fact(is_transposed?size()-id-1:id);

			switch (factor_generic->getType())
			{
				case MatType::Dense :
					{
						MatDense<FPP,Cpu>* factor_dense_ptr = dynamic_cast<MatDense<FPP,Cpu>* > (factor_generic);
						dense_factor = (*factor_dense_ptr); //copy
					}
					break;

				case MatType::Sparse :
					{
						MatSparse<FPP,Cpu>* factor_sparse_ptr = dynamic_cast<MatSparse<FPP,Cpu>* > (factor_generic);
						dense_factor = (*factor_sparse_ptr); //copy
					}
					break;

				default:
					handleError("Faust::TransformHelper", "get_fact : unknown type of the factor matrix");
			}
			delete factor_generic;
			if(is_transposed) dense_factor.transpose();
			if(is_conjugate) dense_factor.conjugate();
			return dense_factor;
		}

	template<typename FPP>
		bool Faust::TransformHelper<FPP,Cpu>::is_fact_sparse(const faust_unsigned_int id) const
		{
			return this->transform->is_fact_sparse(is_transposed?size()-id-1:id);
		}


	template<typename FPP>
		bool Faust::TransformHelper<FPP,Cpu>::is_fact_dense(const faust_unsigned_int id) const
		{
			return this->transform->is_fact_dense(is_transposed?size()-id-1:id);
		}

	template<typename FPP>
	void TransformHelper<FPP, Cpu>::copy_slices(TransformHelper<FPP, Cpu> *th, const bool transpose /* default to false */)
	{
		this->slices[0].copy(th->slices[0]);
		this->slices[1].copy(th->slices[1]);
	}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP, Cpu>::slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
				faust_unsigned_int start_col_id, faust_unsigned_int end_col_id)
		{
			Slice sr(start_row_id, end_row_id);
			Slice sc(start_col_id, end_col_id);
			Slice s[2];
			if(is_transposed)
			{
				s[0] = sc;
				s[1] = sr;
			}
			else
			{
				s[0] = sr;
				s[1] = sc;
			}
			return new TransformHelper<FPP, Cpu>(this, s);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::fancy_index(faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols)
		{
			return new TransformHelper<FPP,Cpu>(this, row_ids, num_rows, col_ids, num_cols);
		}

	template<typename FPP>
	void TransformHelper<FPP, Cpu>::eval_fancy_idx_Transform()
	{
		bool cloning_fact = false;
		faust_unsigned_int size = this->size();
		std::vector<MatGeneric<FPP,Cpu>*> factors((size_t) size);
		MatGeneric<FPP,Cpu>* gen_fac, *first_sub_fac, *last_sub_fac;
		gen_fac = this->transform->get_fact(0, cloning_fact);
		//				first_sub_fac = gen_fac->get_rows(slices[0].start_id, slices[0].end_id-slices[0].start_id);
		//		first_sub_fac->Display();
		first_sub_fac = gen_fac->get_rows(this->fancy_indices[0], this->fancy_num_rows);
		if(cloning_fact)
			delete gen_fac;
		if(size > 1) {
			gen_fac = this->transform->get_fact(size-1, cloning_fact);
			//					last_sub_fac = gen_fac->get_cols(slices[1].start_id, slices[1].end_id-slices[1].start_id);
			last_sub_fac = gen_fac->get_cols(this->fancy_indices[1], this->fancy_num_cols);	//		std::cout << "---" << std::endl;
			//		last_sub_fac->Display();
			if(cloning_fact)
				delete gen_fac;
			factors.reserve(size);
			factors.insert(factors.begin(), first_sub_fac);
			if(size > 2)
			{
				auto it = factors.begin();
				for(faust_unsigned_int i = 1; i < size-1; i++)
				{
					gen_fac = this->transform->get_fact(i, cloning_fact);
					factors[i] = gen_fac;
				}
			}
			factors.insert(factors.begin()+(size-1), last_sub_fac);
			factors.resize(size);
		}
		else { //only one factor
			last_sub_fac = first_sub_fac->get_cols(this->fancy_indices[1], this->fancy_num_cols);
			delete first_sub_fac;
			factors[0] = last_sub_fac;
			factors.resize(1);
		}
		this->transform = make_shared<Transform<FPP, Cpu>>(factors, 1.0, false, cloning_fact);
		if(cloning_fact) {
			for(faust_unsigned_int i = 0; i < size; i++)
				delete factors[i];
		}
	}

	template<typename FPP>
	void TransformHelper<FPP, Cpu>::eval_sliced_Transform()
	{
		bool cloning_fact = true;
		std::vector<MatGeneric<FPP,Cpu>*> factors((size_t) this->size());
		faust_unsigned_int size = this->size();
		MatGeneric<FPP,Cpu>* gen_fac, *first_sub_fac, *last_sub_fac;
		gen_fac = this->transform->get_fact(0, cloning_fact);
		first_sub_fac = gen_fac->get_rows(slices[0].start_id, slices[0].end_id-slices[0].start_id);
		//		first_sub_fac->Display();
		//
		if(cloning_fact)
			delete gen_fac;
		if(size > 1) {
			gen_fac = this->transform->get_fact(size-1, cloning_fact);
			last_sub_fac = gen_fac->get_cols(slices[1].start_id, slices[1].end_id-slices[1].start_id);
			//		std::cout << "---" << std::endl;
			//		last_sub_fac->Display();
			if(cloning_fact)
				delete gen_fac;
			factors.reserve(size);
			factors.insert(factors.begin(), first_sub_fac);
			if(size > 2)
			{
				for(faust_unsigned_int i = 1; i < size-1; i++)
				{
					gen_fac = this->transform->get_fact(i, cloning_fact);
					factors[i] = gen_fac;
				}

			}
			factors.insert(factors.begin()+(size-1), last_sub_fac);
			factors.resize(size);
		}
		else { //only one factor
			last_sub_fac = first_sub_fac->get_cols(slices[1].start_id, slices[1].end_id-slices[1].start_id);
			delete first_sub_fac;
			factors[0] = last_sub_fac;
			factors.resize(1);
		}
		this->transform = make_shared<Transform<FPP,Cpu>>(factors, 1.0, false, cloning_fact);
		if(cloning_fact) {
			for(faust_unsigned_int i = 0; i < size; i++)
				delete factors[i];
		}
	}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_product() const {
			return this->transform->get_product(isTransposed2char(), is_conjugate);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::save_mat_file(const char* filename) const
		{
			this->transform->save_mat_file(filename, is_transposed, is_conjugate);
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::spectralNorm(const int nbr_iter_max, double threshold, int &flag) const
		{
			vector <MatGeneric<FPP, Cpu>*>& orig_facts = this->transform->data;
			int start_id, end_id;
			this->transform->get_nonortho_interior_prod_ids(start_id, end_id);
//			cout << "start_id=" << start_id << "end_id=" << end_id << endl;
			if(start_id < 0)
				return 1.0;
			else if(start_id == 0)
				return this->transform->spectralNorm(nbr_iter_max, threshold, flag);
//			cout << "optimized norm2" << endl;
			vector<MatGeneric<FPP,Cpu>*> non_ortho_start(orig_facts.begin()+start_id, orig_facts.end());
			TransformHelper<FPP, Cpu> t(non_ortho_start, 1.0, false, false);
			return t.transform->spectralNorm(nbr_iter_max, threshold, flag);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::transpose()
		{
			return new TransformHelper<FPP,Cpu>(this, true, false);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::vertcat(const TransformHelper<FPP,Cpu>* G)
		{
			// NOTE: this function is in this class and not in Faust::Transform
			// because it's simpler to handle the transpose and conjugate cases here
			// the functions get_gen_fact() and get_fact_nb_cols()/rows are
			// really helpful for that.
			//TODO: too big function, refactor
			//TODO: clean the code
			TransformHelper<FPP,Cpu>* T = nullptr;
			// take the greater number of factors from this and G as the concatened Faust number
			// of factors
			std::vector<MatGeneric<FPP,Cpu>*> facts(max(G->size(),this->size())+1);
			const MatSparse<FPP,Cpu> * F_fac, *G_fac, *tmp_fac;
			const MatGeneric<FPP,Cpu> *tmp_fac_gen;
			MatSparse<FPP,Cpu>* T_fac;
			faust_unsigned_int T_fac_nnz;
			faust_unsigned_int T_fac_nb_rows, T_fac_nb_cols;
			bool F_inter_fac_allocated, G_inter_fac_allocated;
			bool F_last_fac=false, G_last_fac=false;
			int * colind, *rowptr;
			FPP* values;

			if(this->getNbCol() != G->getNbCol()) handleError("TransformHelper::vertcat()","The dimensions must agree.");
			for(faust_unsigned_int i=0; i < facts.size(); i++){
//				cout << "factor i=" << i <<  " facts.size()=" << facts.size() << endl;
				// F (*this) factor
				if(i < this->size())
				{
					tmp_fac_gen = get_gen_fact(i); //this->transform->data[i];
					if(F_inter_fac_allocated = (tmp_fac_gen->getType() == MatType::Dense))
					{
						F_fac = new MatSparse<FPP,Cpu>(*dynamic_cast<const MatDense<FPP,Cpu>*>(tmp_fac_gen));
						if(is_transposed) const_cast<MatSparse<FPP,Cpu>*>(F_fac)->transpose();
						if(is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(F_fac)->conjugate();
					}
					else
					{
						F_fac = dynamic_cast<const MatSparse<FPP,Cpu>*>(tmp_fac_gen);
						if(is_transposed || is_conjugate)
						{
							F_fac = new MatSparse<FPP,Cpu>(F_fac->nnz, F_fac->getNbRow(), F_fac->getNbCol(), F_fac->getValuePtr(), F_fac->getRowPtr(), F_fac->getColInd(), is_transposed);
							if(is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(F_fac)->conjugate();

							F_inter_fac_allocated = true;
						}
					}
					T_fac_nnz = F_fac->getNonZeros(); //+G_fac->getNonZeros()
					T_fac_nb_cols = this->get_fact_nb_cols(i);//F_fac->getNbCol(); //+G_fac part
					T_fac_nb_rows = this->get_fact_nb_rows(i);//F_fac->getNbRow(); //+G_fac part
				}
				else
				{
					//fill remaining factors by identity
					tmp_fac_gen = get_gen_fact(this->size()-1);//this->transform->data[this->size()-1];
					// number of ones
					T_fac_nnz = this->get_fact_nb_cols(this->size()-1);//tmp_fac_gen->getNbCol(); //+G_fac->getNonZeros()
					T_fac_nb_rows = T_fac_nb_cols = T_fac_nnz; //+G_fac part
					// opt. create the id factor only once
					if(!F_last_fac)
					{
						F_last_fac = true;
						F_fac = MatSparse<FPP,Cpu>::eye(T_fac_nnz, T_fac_nnz);
					}
				}
				// G factor
				if(i < G->size())
				{
					tmp_fac_gen = G->get_gen_fact(i);//G->transform->data[i];
					if((G_inter_fac_allocated = (tmp_fac_gen->getType() == MatType::Dense)))
					{
						G_fac = new MatSparse<FPP,Cpu>(*dynamic_cast<const MatDense<FPP,Cpu>*>(tmp_fac_gen));
						if(G->is_transposed) const_cast<MatSparse<FPP,Cpu>*>(G_fac)->transpose();
						if (G->is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(G_fac)->conjugate();

					}
					else
					{
						G_fac = dynamic_cast<const MatSparse<FPP,Cpu>*>(tmp_fac_gen);
						if(G->is_transposed || G->is_conjugate)
						{
							G_fac = new MatSparse<FPP,Cpu>(G_fac->nnz, G_fac->getNbRow(), G_fac->getNbCol(), G_fac->getValuePtr(), G_fac->getRowPtr(), G_fac->getColInd(), G->is_transposed);
							if (G->is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(G_fac)->conjugate();
							G_inter_fac_allocated = true;
						}

					}
					T_fac_nnz += G_fac->getNonZeros();
					T_fac_nb_cols += G->get_fact_nb_cols(i);//G_fac->getNbCol();
					T_fac_nb_rows += G->get_fact_nb_rows(i);//G_fac->getNbRow();
				}
				else
				{ //G has less or equal factors than F

					//fill remaining factors by identity
					tmp_fac_gen = G->get_gen_fact(G->size()-1);//G->transform->data[G->size()-1];
					// number of ones
					T_fac_nnz += G->get_fact_nb_cols(G->size()-1);//tmp_fac_gen->getNbCol();
					if(i < facts.size()-1 || !F_last_fac) {
						T_fac_nb_cols = F_fac->getNbCol()+G->get_fact_nb_cols(G->size()-1);//F_fac->getNbCol()+tmp_fac_gen->getNbCol();
					}
					//G_fac will be square id. => nb cols == nb rows
					T_fac_nb_rows = F_fac->getNbRow()+G->get_fact_nb_cols(G->size()-1);//tmp_fac_gen->getNbCol();
					// opt. create the id factor only once
					if(!G_last_fac)
					{
						G_last_fac = true;
						G_fac = MatSparse<FPP,Cpu>::eye(G->get_fact_nb_cols(G->size()-1),G->get_fact_nb_cols(G->size()-1));//(tmp_fac_gen->getNbCol(), tmp_fac_gen->getNbCol());
					}
				}
				values = new FPP[T_fac_nnz];
				colind = new int[T_fac_nnz];
				rowptr = new int[T_fac_nb_rows+1];
				// group the data from F and G
				memcpy(values, F_fac->getValuePtr(), sizeof(FPP)*F_fac->getNonZeros());
				memcpy(values+F_fac->getNonZeros(), G_fac->getValuePtr(),
						sizeof(FPP)*G_fac->getNonZeros());
				assert(T_fac_nnz == F_fac->getNonZeros()+G_fac->getNonZeros());
				assert(T_fac_nb_rows == F_fac->getNbRow()+G_fac->getNbRow());
				/****** col indices ******/
				// F indices don't change
				memcpy(colind, F_fac->getColInd(), sizeof(int)*F_fac->getNonZeros());
				// G indices are shifted by F->getNbCol()
				memcpy(colind+F_fac->getNonZeros(), G_fac->getColInd(), sizeof(int)*G_fac->getNonZeros());
				//shift G indices
				if(i < facts.size()-1)
					for(faust_unsigned_int j=F_fac->getNonZeros();j<F_fac->getNonZeros()+G_fac->getNonZeros();j++)
						colind[j] += F_fac->getNbCol();
				/***** row indices *****/
				memcpy(rowptr, F_fac->getRowPtr(), sizeof(int)*(F_fac->getNbRow()+1));
				//ignore first ele == 0 of G_fac->getRowPtr()
				memcpy(rowptr+F_fac->getNbRow()+1, G_fac->getRowPtr()+1, sizeof(int)*G_fac->getNbRow());
				// shift number of elements to take account of already added F_fac elts
				for(faust_unsigned_int j=F_fac->getNbRow()+1;j<F_fac->getNbRow()+G_fac->getNbRow()+1;j++)
					rowptr[j] += F_fac->getRowPtr()[F_fac->getNbRow()];
				// concatened Faust factor
//				cout << "T_fac_nb_rows:" << T_fac_nb_rows << endl;
//				cout << "T_fac_nb_cols:" << T_fac_nb_cols << endl;
//				cout << "nnz:" << T_fac_nnz << endl;
//				cout << "rowptr=";
//				for(faust_unsigned_int j=0;j<T_fac_nb_rows+1;j++)
//					cout << rowptr[j] << ",";
//				cout << endl;
//				cout << "colind=";
//				for(faust_unsigned_int j=0;j<T_fac_nnz;j++)
//					cout << colind[j] << ",";
//				cout << endl;
//				cout << "values=";
//				for(faust_unsigned_int j=0;j<T_fac_nnz;j++)
//					cout << values[j] << ",";
//				cout << endl;
				T_fac = new MatSparse<FPP,Cpu>(T_fac_nnz, T_fac_nb_rows, T_fac_nb_cols, values, rowptr, colind);
//				cout << "stage 4: ok" << endl;
//				cout << "T_fac:"<< endl;
//				T_fac->Display();
				facts[i] = T_fac;
				if(F_inter_fac_allocated)
				{
					delete F_fac;
					F_inter_fac_allocated = ! F_inter_fac_allocated;
				}
				if(G_inter_fac_allocated)
				{
					delete G_fac;
					G_inter_fac_allocated = ! G_inter_fac_allocated;
				}
				delete values;
				delete colind;
				delete rowptr;
			}
			T = new TransformHelper<FPP,Cpu>(facts, 1.0, false, false);
//			cout << "final stage ok" << endl;
//			T->display();
			//delete last factors (identity for each Faust)
			if(!F_inter_fac_allocated) delete F_fac;
			if(!F_inter_fac_allocated) delete G_fac;
			return T;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::horzcat(const TransformHelper<FPP,Cpu>* G)
		{
			TransformHelper<FPP,Cpu> *Ft, *Gt, *C, *Ct;
			Ft = this->transpose();
			Gt = const_cast<TransformHelper<FPP,Cpu> *>(G)->transpose(); //no harm
			C = Ft->vertcat(Gt);
			Ct = C->transpose();
			delete Ft;
			delete Gt;
			delete C;
			return Ct;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::conjugate()
		{
			return new TransformHelper<FPP,Cpu>(this, false, true);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::adjoint()
		{
			return new TransformHelper<FPP,Cpu>(this, true, true);
		}


	template<typename FPP>
		bool TransformHelper<FPP,Cpu>::isTransposed() const
		{
			return is_transposed;
		}

	template<typename FPP>
		bool TransformHelper<FPP,Cpu>::isConjugate() const
		{
			return is_conjugate;
		}

	template<typename FPP>
		const char TransformHelper<FPP,Cpu>::isTransposed2char() const
		{
			return is_transposed?'T':'N';
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::normL1() const {
			return this->transform->normL1(is_transposed);
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::normInf() const {
			return this->transform->normL1(!is_transposed);
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::normFro() const {
			vector <MatGeneric<FPP, Cpu>*>& orig_facts = this->transform->data;
			int start_id, end_id;
			this->transform->get_nonortho_interior_prod_ids(start_id, end_id);
			if(start_id < 0)
				return Faust::fabs(MatDense<FPP,Cpu>::eye(this->getNbCol(), this->getNbCol()).norm());
			else if(start_id == 0)
				return this->transform->normFro();
			vector<MatGeneric<FPP,Cpu>*> non_ortho_start(orig_facts.begin()+start_id, orig_facts.end());
			TransformHelper<FPP, Cpu> t(non_ortho_start, 1.0, false, false);
			return t.transform->normFro();
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::normalize(const int meth /* 1 for 1-norm, 2 for 2-norm (2-norm), MAX for inf-norm */) const
		{
			const int MAX = (1u << 31) - 1;// numeric_limits<int>::max();
			const MatGeneric<FPP,Cpu>* last_fact = get_gen_fact(this->size()-1);
			unsigned int ncols = this->getNbCol();
			unsigned int nrows = this->getNbRow();
			//2-norm parameters
			double precision =  0.001;
			faust_unsigned_int nbr_iter_max=100;
			int flag;

			vector<FPP> norm_invs(ncols);
			FPP norm;
			vector<int> coords(ncols);
			TransformHelper<FPP,Cpu>* normalizedTh = nullptr;
			//			last_fact->Display();
			for(faust_unsigned_int i=0;i<ncols;i++)
			{
				//TODO: we shouldn't have to const_cast, slice() must be const
				const TransformHelper<FPP,Cpu> * col = const_cast<TransformHelper<FPP,Cpu> *>(this)->slice(0,nrows, i, i+1);
//				cout << "TransformHelper normalize, meth=" << meth << endl;
				switch(meth)
				{
					case 1: //1-norm
//						cout << "normalize with 1-norm" << endl;
						norm = col->normL1();
						break;
					case 2: //2-norm
//						cout << "normalize with 2-norm" << endl;
						norm = col->spectralNorm(nbr_iter_max, precision, flag);
						break;
					case -2: // fro norm
//						cout << "normalize with fro-norm" << endl;
						norm = col->normFro();
						break;
					case -1: //inf norm
//						cout << "normalize with inf-norm" << endl;
						norm = col->normInf();
						break;
					default:
						handleError("Faust::TransformHelper::normalize()", "order for the norm to use is not valid");
				}
				if(norm != FPP(0))
					norm_invs[i] = 1./norm;
				else
					norm_invs[i] = 1;
				coords[i] = i;
				delete col;
			}
			MatSparse<FPP,Cpu>* norm_diag = new MatSparse<FPP,Cpu>(coords, coords,
					norm_invs, ncols, ncols);
			//			norm_diag.Display();
			vector<MatGeneric<FPP,Cpu>*> factors;
			for(int i =0; i < size(); i++)
				//NOTE: about const cast: no problem, we know we won't write it
				//NOTE: don't use get_gen_fact() to avoid transpose auto-handling
				factors.push_back(const_cast<Faust::MatGeneric<FPP, (Device)0u>*>(this->transform->data[i]));
			if(this->is_transposed)
				factors.insert(factors.begin(), norm_diag);
			else
				factors.push_back(norm_diag);
			normalizedTh = new TransformHelper<FPP,Cpu>(factors, FPP(1.0), false, false);
			normalizedTh->is_transposed = this->is_transposed;
			//			normalizedTh->display();
			return normalizedTh;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::~TransformHelper() {
			// transform is deleted auto. when no TransformHelper uses it (no more weak refs left)
#ifdef FAUST_VERBOSE
			cout << "Destroying Faust::TransformHelper object." << endl;
#endif
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density /* 1.f */, bool per_row /* true */) {
			if(!TransformHelper<FPP,Cpu>::seed_init) {
				srand(time(NULL)); //seed init needed for MatDense rand generation
				TransformHelper<FPP,Cpu>::seed_init = true;

			}
			// pick randomly the number of factors into {min_num_factors, ..., max_num_factors}
			std::uniform_int_distribution<int> num_fac_distr(min_num_factors, max_num_factors);
			std::uniform_int_distribution<int> dim_distr(min_dim_size, max_dim_size);
			std::uniform_int_distribution<int> bin_distr(0,1);
			unsigned int num_factors = num_fac_distr(generator);
			// create factors randomly respecting the RandFaustType asked and dims interval
			std::vector<MatGeneric<FPP,Cpu>*> factors((size_t) num_factors);
			unsigned int num_rows, num_cols = dim_distr(generator);
			float fact_density;
			for(unsigned int i=0;i<num_factors;i++) {
				num_rows = num_cols;
				num_cols = dim_distr(generator);
#ifdef FAUST_VERBOSE
				cout << "TransformHelper<FPP,Cpu>::randFaust() per_row: " <<  per_row << endl;
#endif
				if(density == -1.)
					fact_density = per_row?5.1/num_cols:5.1/num_rows;
				else
					fact_density = density;
				switch(t){
					case DENSE:
						factors[i] = MatDense<FPP,Cpu>::randMat(num_rows, num_cols, fact_density, per_row);
						break;
					case SPARSE:
						factors[i] = MatSparse<FPP,Cpu>::randMat(num_rows, num_cols, fact_density, per_row);
						break;
					case MIXED:
						if(bin_distr(generator))
							factors[i] = MatDense<FPP,Cpu>::randMat(num_rows, num_cols, fact_density, per_row);
						else
							factors[i] = MatSparse<FPP,Cpu>::randMat(num_rows, num_cols, fact_density, per_row);
						break;
					default:
						handleError("Faust::TransformHelper", "randFaust(): Unknown RandFaustType");
						break;

				}
				if(factors[i] == NULL) return NULL;
			}
			TransformHelper<FPP,Cpu>* randFaust = new TransformHelper<FPP, Cpu>(factors,1.0,false,false);
			return randFaust;
		}


	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::hadamardFaust(unsigned int n)
		{
			TransformHelper<FPP,Cpu>* hadamardFaust = nullptr;
			vector<MatGeneric<FPP,Cpu>*> factors;
			bool cloning_fact = false; // big opt. allowed only because of the RefManager used in Transform class
			//this opt. avoids to duplicate the same factor
			try {
				wht_factors(n, factors, cloning_fact);
				hadamardFaust = new TransformHelper<FPP, Cpu>(factors, 1.0, false, false);
			}
			catch(std::bad_alloc)
			{
				// nothing to do, out of memory return nullptr
			}
			return hadamardFaust;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::fourierFaust(unsigned int n)
		{

			vector<MatGeneric<FPP,Cpu>*> factors(n+1);
			TransformHelper<FPP,Cpu>* fourierFaust = nullptr;
			try
			{
				fft_factors(n, factors);
				//			for(int i=0;i<n+1;i++)
				//				factors[i]->Display();
				fourierFaust = new TransformHelper<FPP, Cpu>(factors, 1.0, false, false);
			}
			catch(std::bad_alloc e)
			{
				//nothing to do, out of memory, return nullptr
			}
			return fourierFaust;
		}

	template<typename FPP> bool TransformHelper<FPP,Cpu>::seed_init = false;
	template<typename FPP> std::default_random_engine TransformHelper<FPP,Cpu>::generator(time(NULL));
}
