%=========================================
%> @brief Constructs a toeplitz Faust whose first column is c and first row r.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b T = toeplitz(c), T is a symmetric Toeplitz Faust whose the first column is c. <br/>
%> &nbsp;&nbsp;&nbsp; @b T = toeplitz(c, r), T is a Toeplitz Faust whose the first column is c and the first row is [c(1), r(2:)] <br/>
%> @param c: the first column of the toeplitz Faust.
%> @param r: (2nd argument) the first row of the toeplitz Faust. Defaulty r = conj(c).
%> r(1) is ignored, the first row is always [c(1), r(2:)].
%>
%>
%> @retval T the toeplitz Faust.
%>
%> @b Example
%>
%> @code
%> % in a matlab terminal
%> >> import matfaust.toeplitz
%> >> c = 1:10;
%> >> T = toeplitz(c)
%>
%> T =
%> @endcode
%>
%> Faust size 10x10, density 5.52, nnz_sum 552, 10 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 10x32, density 0.0625, nnz 20, addr: 0x7f278c6aec20
%> - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278ded1440
%> - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278defd560
%> - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278f687dc0
%> - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f27899a94d0
%> - FACTOR 5 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f2724d04e30
%> - FACTOR 6 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278deb8fb0
%> - FACTOR 7 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278f69e8f0
%> - FACTOR 8 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278ded1840
%> - FACTOR 9 (complex) SPARSE, size 32x10, density 0.0625, nnz 20, addr: 0x7f278c6b2070
%>
%> @code
%> >> full_T = full(T);
%> >> all(full_T(:,1).' - c < 1e-15)
%>
%> ans =
%>
%>   logical
%>
%>      1
%>
%> >> all(full_T(1, :) - c < 1e-15)
%>
%> ans =
%>
%>   logical
%>
%>      1
%> @endcode
%>
%> @code
%> >> r = 11:20;
%> >> T2 = toeplitz(c, r)
%> @endcode
%>
%> T2 =
%>
%> Faust size 10x10, density 5.52, nnz_sum 552, 10 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 10x32, density 0.0625, nnz 20, addr: 0x7f278dee2640
%> - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278c694ff0
%> - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278debd710
%> - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278deb5e70
%> - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278c6dc0d0
%> - FACTOR 5 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278decd230
%> - FACTOR 6 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278deb44a0
%> - FACTOR 7 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278dee9c00
%> - FACTOR 8 (complex) SPARSE, size 32x32, density 0.0625, nnz 64, addr: 0x7f278deb6090
%> - FACTOR 9 (complex) SPARSE, size 32x10, density 0.0625, nnz 20, addr: 0x7f278c682e60
%>
%> @code
%> >> full_T2 = full(T2);
%> >> all(full_T2(1, :) - [c(1) , r(2:end)] < 1e-15)
%>
%> ans =
%>
%>   logical
%>
%>    1
%>
%> >> all(full_T2(:, 1).' - c < 1e-14)
%>
%> ans =
%>
%>   logical
%>
%>    1
%> >> % Look at the density of a larger toeplitz Faust
%> >> % it indicates a speedup of the Faust-matrix/vector product
%> >> density(toeplitz(rand(1, 1024)))
%> 0.0820
%> @endcode
%>
%>
%> @b See also matfaust.circ, matfaust.anticirc
%=========================================
function T = toeplitz(c, varargin)
    if (length(varargin) > 0)
        r = varargin{1};
        if(~ ismatrix(r) || size(r, 1) ~= 1 && size(r, 2) ~= 1)
            error('The second argument must be a vector')
        end
    else
        r = conj(c); % default r
    end
    if ~ ismatrix(r) || ~ ismatrix(c) || ~ isnumeric(c) || ~ isnumeric(r)
        error('r and c must be numeric vectors')
    end
    if size(c, 2) == 1
        c = c.';
    elseif size(c, 1) ~= 1
        error('c must be a vector')
    end
    if size(r, 2) == 1
        r = r.';
    elseif size(r, 1) ~= 1
        error('r must be a vector')
    end
    m = numel(c);
    n = numel(r);
    N = 2 ^ ceil(log2(max(m, n)));
    c_ = [c, zeros(1, N-m+1+N-n), r(end:-1:2)];
    C = matfaust.circ(c_);
    T = C(1:m, 1:n);
end
