namespace scipy
{
	template<class I, class J>
		void csr_column_index1(const J n_idx,
				const J col_idxs[],
				const J n_row,
				const J n_col,
				const I Ap[],
				const I Aj[],
				I col_offsets[],
				I Bp[])
		{
			// bincount(col_idxs)
			for(I jj = 0; jj < n_idx; jj++){
				const I j = col_idxs[jj];
				col_offsets[j]++;
			}

			// Compute new indptr
			I new_nnz = 0;
			Bp[0] = 0;
			for(I i = 0; i < n_row; i++){
				for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
					new_nnz += col_offsets[Aj[jj]];
				}
				Bp[i+1] = new_nnz;
			}

			// cumsum in-place
			for(I j = 1; j < n_col; j++){
				col_offsets[j] += col_offsets[j - 1];
			}
		}

	template<class I, class J, class T>
		void csr_column_index2(const J col_order[],
				const I col_offsets[],
				const J nnz,
				const I Aj[],
				const T Ax[],
				I Bj[],
				T Bx[])
		{
			I n = 0;
			for(I jj = 0; jj < nnz; jj++){
				const I j = Aj[jj];
				const I offset = col_offsets[j];
				const I prev_offset = j == 0 ? 0 : col_offsets[j-1];
				if (offset != prev_offset) {
					const T v = Ax[jj];
					for(I k = prev_offset; k < offset; k++){
						Bj[n] = col_order[k];
						Bx[n] = v;
						n++;
					}
				}
			}
		}

	template<class I, class J, class T>
		void csr_row_index(const J n_row_idx,
				const J rows[],
				const I Ap[],
				const I Aj[],
				const T Ax[],
				I Bj[],
				T Bx[])
		{
			for(J i = 0; i < n_row_idx; i++){
				const J row = rows[i];
				const I row_start = Ap[row];
				const I row_end   = Ap[row+1];
				Bj = std::copy(Aj + row_start, Aj + row_end, Bj);
				Bx = std::copy(Ax + row_start, Ax + row_end, Bx);
			}
		}
}
