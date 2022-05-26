%==========================================================================
%> @brief Factorizes the matrix M according to a butterfly support.

%>        @param 'type', str: the type of factorization 'right'ward, 'left'ward or 'bbtree'.
%>        More precisely: if 'left' (resp. 'right') is used then at each stage of the
%>        factorization the most left factor (resp. the most right factor) is split in two.
%>        If 'bbtree' is used then the matrix is factorized according to a Balanced
%>        Binary Tree (which is faster as it allows parallelization).
%>
%> @retval F the Faust which is an approximate of M according to a butterfly support.
%>
%> @b Example:
%> @code
%> >> import matfaust.wht
%> >> import matfaust.dft
%> >> import matfaust.fact.butterfly
%> >> M = full(wht(32)); % it works with dft too!
%> >> F = butterfly(M, 'type', 'bbtree');
%> >> err = norm(full(F)-M)/norm(M)
%> err =
%>
%>    1.4311e-15
%> @endcode
%>
%> <b>Reference:</b> Leon Zheng, Elisa Riccietti, and Remi Gribonval, <a href="https://arxiv.org/pdf/2110.01230.pdf">Hierarchical Identifiability in Multi-layer Sparse Matrix Factorization</a>
%==========================================================================
function F = butterfly(M, varargin)
        import matfaust.Faust
        nargin = length(varargin);
        type = 'right';
		perm = [];
        if(nargin > 0)
                for i=1:2:nargin
                        switch(varargin{i})
                                case 'type'
                                        if(nargin < i+1 || ~ any(strcmp(varargin{i+1}, {'right', 'left', 'bbtree'})))
                                                error('keyword argument ''type'' must be followed by ''left'' or ''right'' or ''bbtree''')
                                        else
                                                type = varargin{i+1};
                                        end
								case 'perm'
									if(nargin < i+1 ||  ~ is_array_of_indices(varargin{i+1}, M) && ~ is_cell_arrays_of_indices(varargin{i+1}, M) && ~ strcmp(varargin{i+1}, 'default_8'))
                                                error('keyword argument ''perm'' must be followed by ''default_8'', an array of permutation indices or a cell array of arrays of permutation indices')
                                        else
                                                perm = varargin{i+1};
                                        end
                        end
                end
        end
		if iscell(perm) % perm is a cell of arrays, each one defining a permutation to test
			% evaluate butterfly factorisation using the permutations and
			% keep the best Faust
			min_err = inf;
			nM = norm(M);
			for i=1:length(perm)
				F = matfaust.fact.butterfly(M, 'type', type, 'perm', perm{i});
				err = norm(full(F)-M)/nM;
				if err < min_err
					min_err = err;
					best_F = F;
				end
			end
			F = best_F;
			return;
		end
        if(strcmp(type, 'right'))
                type = 1;
        elseif(strcmp(type, 'left'))
                type = 0;
        elseif(strcmp(type, 'bbtree'))
                type = 2;
        end
		if(strcmp(class(M), 'single'))
				core_obj = mexButterflyRealFloat(M, type, perm);
				F = Faust(core_obj, isreal(M), 'cpu', 'float');
		else
			if(isreal(M))
				core_obj = mexButterflyReal(M, type, perm);
			else
				core_obj = mexButterflyCplx(M, type, perm);
			end
			F = Faust(core_obj, isreal(M));
		end
end

function b = is_array_of_indices(value, M)
	s = size(M,2);
	b = ismatrix(value) && isnumeric(value) && isreal(value);
	b = b && all(value <= s);
	b = b && all(value > 0);
	b = b && length(value) == s;
end

function b = is_cell_arrays_of_indices(value, M)
	b = iscell(value);
	for i=1:length(value)
		b = b && is_array_of_indices(value{i}, M);
		if ~ b
			break
		end
	end
end
