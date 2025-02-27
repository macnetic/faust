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
%> @param 'dev', str: 'gpu' or 'cpu' to create the Faust on CPU or GPU ('cpu' is the default).
%> @param 'diag_opt', logical: cf. matfaust.circ.
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
    import matfaust.Faust
	dev = 'cpu';
    diag_opt = false;
    argc = length(varargin);
    if (argc > 0)
        i = 1;
        if ismatrix(varargin{1}) && isnumeric(varargin{1})
            r = varargin{1};
            if(size(r, 1) ~= 1 && size(r, 2) ~= 1)
                error('The second argument must be a vector')
            end
            i = 2;
        else
            r = conj(c); % default r
        end
		if(argc > 1)
            while i <= argc
                %for i=2:2:argc
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
                    case 'diag_opt'

                        if(argc == i || ~ islogical(tmparg))
                            error('diag_opt keyword argument is not followed by a logical')
                        else
                            diag_opt = tmparg;
                        end
                    otherwise
                        if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
                            error([ tmparg ' unrecognized argument'])
                        end
                end
                i = i + 2;
            end
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
    C = matfaust.circ(c_, 'diag_opt', diag_opt);
    if diag_opt
        % see issue #335
        T = Faust(speye(m, size(C, 1))) * C * Faust(speye(size(C, 2), n))
    else
        T = C(1:m, 1:n);
    end
    if startsWith(dev, 'gpu')
	    T = clone(T, 'dev', 'gpu');
    end
end
