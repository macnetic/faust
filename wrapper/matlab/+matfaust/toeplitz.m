%=========================================
%> @brief Returns  a toeplitz Faust whose the first column is c and the first row r.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b T = toeplitz(c), T is a symmetric Toeplitz Faust whose the first column is c. <br/>
%> &nbsp;&nbsp;&nbsp; @b T = toeplitz(c, r), T is a Toeplitz Faust whose the first column is c and the first row is [c(1), r(2:)] <br/>

%> @param c: the first column of the toeplitz Faust.
%> @param r: (2nd argument) the first row of the toeplitz Faust. Defaulty r = c.
%> r(1) is ignored, the first row is always [c(1),
%>         r(2:)].
%>
%>
%> @b Example
%>
%> @code
%> % in a matlab terminal
%> >> c = rand(1, 10);
%> >> T = toeplitz(c)
%>
%> T =
%> @endcode
%>
%> Faust size 10x10, density 6.16, nnz_sum 616, 12 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 10x32, density 0.0625, nnz 20
%> - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 5 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
%> - FACTOR 6 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
%> - FACTOR 7 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 8 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 9 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 10 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 11 (complex) SPARSE, size 32x10, density 0.0625, nnz 20
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
%> >> r = rand(1, 10);
%> >> T2 = toeplitz(c, r)
%> @endcode
%>
%> T2 =
%>
%> Faust size 10x10, density 6.16, nnz_sum 616, 12 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 10x32, density 0.0625, nnz 20
%> - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 5 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
%> - FACTOR 6 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
%> - FACTOR 7 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 8 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 9 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 10 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 11 (complex) SPARSE, size 32x10, density 0.0625, nnz 20
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
%> >> all(full_T2(:, 1).' - c < 1e-15)
%>
%> ans =
%>
%>   logical
%>
%>    1
%> @endcode
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
        r = c; % default r
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
