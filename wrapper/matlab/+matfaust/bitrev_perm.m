%==========================================================================================
%> @brief Bitreversal permutation.
%>
%> @param n: the size of the permutation, it must be a power of two. P dimensions will be n x n.
%>
%> @retval P a sparse matrix defining the bit-reversal permutation.
%>
%> @b See also matfaust.dft
%==========================================================================================
function P = bitrev_perm(n)
    index = 1:n;
    new_index = BitReversalPermutation(index);
    P = sparse(index, new_index, ones(1, n), n, n);
end
