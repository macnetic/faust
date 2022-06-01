%==========================================================================
%> @brief Factorizes the matrix M according to a butterfly support.
%>        @param 'type', str: the type of factorization 'right'ward, 'left'ward or 'bbtree'.
%>        More precisely: if 'left' (resp. 'right') is used then at each stage of the
%>        factorization the most left factor (resp. the most right factor) is split in two.
%>        If 'bbtree' is used then the matrix is factorized according to a Balanced
%>        Binary Tree (which is faster as it allows parallelization).
%>
%>
%> @param 'perm', value	four kind of values are possible for this argument (Note that this argument is made only for the bbtree type of factorization).
%>
%> 1. perm is an array of column indices of the permutation matrix P which is such that the returned Faust is F = G * P.' where G is the Faust butterfly approximation of M*P.
%> 2. perm is a cell array of arrays of permutation column indices as defined in 1. In that case, all permutations passed to the function are used as explained in 1, each one producing a Faust, the best one (that is the best approximation of M) is kept and returned by butterfly.
%> 3. perm is 'default_8', this is a particular case of 2. Eight default permutations are used. For the definition of those permutations please refer to [1].
%> 4. By default this argument is empty, no permutation is used.
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
%> Use butterfly with a permutation factor defined by J:
%> @code
%> >> J = 32:-1:1;
%> >> F = butterfly(H, 'type', 'bbtree', 'perm', J);
%> F =
%>
%> Faust size 32x32, density 0.34375, nnz_sum 352, 6 factor(s):
%> - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 5 (double) SPARSE, size 32x32, density 0.03125, nnz 32
%> @endcode
%>
%> Use butterfly with successive permutations J1 and J2
%> and keep the best approximation:
%>
%> @code
%> >> J1 = J;
%> >> J2 = randperm(32);
%> >> F = butterfly(H, 'type', 'bbtree', 'perm', {J1, J2})
%>
%> F =
%>
%> Faust size 32x32, density 0.34375, nnz_sum 352, 6 factor(s):
%> - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 5 (double) SPARSE, size 32x32, density 0.03125, nnz 32
%> @endcode
%>
%>
%>
%> <b>Reference: [1]</b> Leon Zheng, Elisa Riccietti, and Remi Gribonval, <a href="https://arxiv.org/pdf/2110.01230.pdf">Hierarchical Identifiability in Multi-layer Sparse Matrix Factorization</a>
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
        if strcmp(perm, 'default_8')
            permutations = cell(1, 8);
            pchoices = {'000', '001', '010', '011', '100', '101', '110', '111'};
            for i=1:8
                 P = get_permutation_matrix(floor(log2(size(M, 1))), pchoices{i});
                 [permutations{i}, ~, ~] = find(P.'); % don't get the column indices directly because it would always 1 to size(P, 1) (indeed P is in CSC format), rather get them in the proper order (row 0 to size(P, 1)) by getting the row indices of P transpose
                 permutations{i} = permutations{i}.'; % just for readibility in case of printing
            end
            F = matfaust.fact.butterfly(M, 'type', type, 'perm', permutations);
            return;
        elseif iscell(perm) % perm is a cell of arrays, each one defining a permutation to test
            % evaluate butterfly factorisation using the permutations and
            % keep the best Faust
            min_err = inf;
            nM = norm(M, 'fro');
            for i=1:length(perm)
				%                perm{i}
				m = numel(perm{i});
				F = matfaust.fact.butterfly(M, 'type', type, 'perm', perm{i});
				err = norm(F-M, 'fro')/nM;
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

function perm = perm_type(i, type)
	%% Type 0 is c in paper. Type 1 is b in paper. Type 2 is a in paper.
	size = 2^i;
	size_o2 = size / 2;
    switch(type)
        case 0
            row_inds = [1:size_o2, size_o2+1:size];
            col_inds = [1:size_o2, size:-1:size-size/2+1];
        case 1
            row_inds = [size_o2:-1:1, size_o2 + (1:size_o2)];
            col_inds = [1:size_o2, size_o2 + (1:size_o2)];
        case 2
            row_inds = [1:size_o2, size_o2 + (1:size_o2)];
            col_inds = [1:2:size-1, 2:2:size];
        otherwise
            error('perm_type received an invalid value for type argument')
    end
	m = numel(row_inds);
	n = m;
	nnz = m;
	perm = sparse(row_inds, col_inds, ones(1, m), m, n, nnz);
end

function permutations = shared_logits_permutation(num_factors, choices)
	% choices: array of three logical-s
	permutations = {};
	for i=2:num_factors
		block = speye(2^i);
		for j=1:3
			if choices(j)
				block = block * perm_type(i, j-1);
			end
		end
		permutations = [permutations,  {kron(speye(2 ^(num_factors - i)), block)}];
	end
end

function p = get_permutation_matrix(num_factors, perm_name)
	% perm_name: str 000, 001, ..., 111
	choices = zeros(1, 3);
	for i=1:length(perm_name)
		choices(i) = str2double(perm_name(i));
	end
	permutations = shared_logits_permutation(num_factors, choices);
	p = sparse(full(matfaust.Faust(permutations)));
end
