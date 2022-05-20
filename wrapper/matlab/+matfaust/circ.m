%==========================================================================================
%> @brief Returns a circulant Faust C defined by the vector c (which is the first column of full(C)).
%>
%> @b Example:
%>
%> @code
%> >> import matfaust.circ
%> >> c = rand(1, 8);
%> >> C = circ(c)
%> @endcode
%>
%> C =
%>
%> Faust size 8x8, density 1.75, nnz_sum 112, 8 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 1 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 2 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 3 (complex) SPARSE, size 8x8, density 0.125, nnz 8
%> - FACTOR 4 (complex) SPARSE, size 8x8, density 0.125, nnz 8
%> - FACTOR 5 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 6 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 7 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%>
%> @code
%> >> full_C = full(C);
%> >> all(full_C(:, 1).' - c < 1e-15)
%>
%> ans =
%>
%>   logical
%>
%>    1
%>
%> >> c
%>
%> c =
%>
%>     0.2630    0.6541    0.6892    0.7482    0.4505    0.0838    0.2290
%>     0.9133
%>
%> >> real(full_C)
%>
%> ans =
%>
%>     0.2630    0.9133    0.2290    0.0838    0.4505    0.7482    0.6892    0.6541
%>     0.6541    0.2630    0.9133    0.2290    0.0838    0.4505    0.7482    0.6892
%>     0.6892    0.6541    0.2630    0.9133    0.2290    0.0838    0.4505    0.7482
%>     0.7482    0.6892    0.6541    0.2630    0.9133    0.2290    0.0838    0.4505
%>     0.4505    0.7482    0.6892    0.6541    0.2630    0.9133    0.2290    0.0838
%>     0.0838    0.4505    0.7482    0.6892    0.6541    0.2630    0.9133    0.2290
%>     0.2290    0.0838    0.4505    0.7482    0.6892    0.6541    0.2630    0.9133
%>     0.9133    0.2290    0.0838    0.4505    0.7482    0.6892    0.6541    0.2630
%> @endcode
%>
%> @b See also matfaust.anticirc, matfaust.toeplitz
%==========================================================================================
function C = circ(c)
    log2c = log2(numel(c));
    if(log2c ~= floor(log2c))
        error('Only power of two length vectors are supported')
    end
    if ~ ismatrix(c) || ~ isnumeric(c)
        error('c must be numeric vector')
    end
    if size(c, 2) == 1
        c = c.';
    elseif size(c, 1) ~= 1
        error('c must be a vector')
    end
    n = numel(c);
    F = matfaust.dft(n, 'normed', false);
    FH = F';
    if (size(c, 1) < size(c, 2))
        c = c.';
    end
    S = sparse(diag(FH*(c/n)));
    C = F * matfaust.Faust(S*factors(FH, 1)) * right(FH, 2);
end
