%==========================================================================================
%> @brief Returns a circulant Faust C defined by the vector c (which is the first column of full(C)).
%>
%> @param c: the vector to define the circulant Faust. Its length must be a power of two.
%> @param 'dev', str: 'gpu' or 'cpu' to create the Faust on CPU or GPU ('cpu' is the default).
%>
%> @b Example:
%>
%> @code
%> >> import matfaust.circ
%> >> c = 1:8;
%> >> C = circ(c)
%> @endcode
%>
%> C =
%>
%> Faust size 8x8, density 1.5, nnz_sum 96, 6 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 1 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 2 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 3 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 4 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 5 (complex) SPARSE, size 8x8, density 0.25, nnz 16
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
%> >> real(full_C)
%>
%> ans =
%>
%>     1.0000    8.0000    7.0000    6.0000    5.0000    4.0000    3.0000    2.0000
%>     2.0000    1.0000    8.0000    7.0000    6.0000    5.0000    4.0000    3.0000
%>     3.0000    2.0000    1.0000    8.0000    7.0000    6.0000    5.0000    4.0000
%>     4.0000    3.0000    2.0000    1.0000    8.0000    7.0000    6.0000    5.0000
%>     5.0000    4.0000    3.0000    2.0000    1.0000    8.0000    7.0000    6.0000
%>     6.0000    5.0000    4.0000    3.0000    2.0000    1.0000    8.0000    7.0000
%>     7.0000    6.0000    5.0000    4.0000    3.0000    2.0000    1.0000    8.0000
%>     8.0000    7.0000    6.0000    5.0000    4.0000    3.0000    2.0000    1.0000
%>
%> >> % Look at the density of a larger circulant Faust
%> >> % it indicates a speedup of the Faust-matrix/vector product
%> >> density(circ(rand(1, 1024)))
%> 0.0391
%> @endcode
%>
%>
%> @b See also matfaust.anticirc, matfaust.toeplitz
%==========================================================================================
function C = circ(c, varargin)
	import matfaust.Faust
	dev = 'cpu';
	argc = length(varargin);
	if(argc > 0)
		for i=1:2:argc
			if(argc > i)
				% next arg (value corresponding to the key varargin{i})
				tmparg = varargin{i+1};
			end
			switch(varargin{i})
				case 'dev'
					if(argc == i || ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error('dev keyword argument is not followed by a valid value: cpu, gpu*.')
					else
						dev = tmparg;
					end
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end
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
%    C = F * matfaust.Faust(S*factors(FH, 1)) * right(FH, 2);
	nf = numfactors(F);
	if(nf > 3)
		C = left(F, nf-1) * Faust(factors(F, nf) * S * factors(FH, 1) * factors(FH, 2)) * right(FH, 3);
	elseif(nf > 2)
		C = left(F, nf-1) * Faust(factors(F, nf) * S * factors(FH, 1) * factors(FH, 2)) * Faust(right(FH, 3));
	else
		C = Faust(left(F, nf-1)) * Faust(factors(F, nf) * S * factors(FH, 1) * factors(FH, 2));
	end
	if startsWith(dev, 'gpu')
		C = clone(C, 'dev', 'gpu');
	end
end
