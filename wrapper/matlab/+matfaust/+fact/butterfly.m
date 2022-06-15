%==========================================================================
%> @brief Factorizes the matrix M according to a butterfly support and optionally a permutation using the algorithms described in [1].
%>
%> The result is a Faust F of the form BP where B has a butterfly structure and P is a permutation matrix determined by the optional parameter ‘perm'.
%>
%> @param M: the matrix to factorize. It can be real (single or double, the class might have a large impact on performance) or complex. M must be square and its dimension must be a power of two.
%>@param 'type', str: the type of factorization 'right'ward, 'left'ward or 'bbtree'.
%>        More precisely: if 'left' (resp. 'right') is used then at each stage of the
%>        factorization the most left factor (resp. the most right factor) is split in two.
%>        If 'bbtree' is used then the matrix is factorized according to a Balanced
%>        Binary Tree (which is faster as it allows parallelization).
%> @param 'perm', value	five kinds of values are possible for this argument.
%>
%> 1. perm is an array of column indices of the permutation matrix P which is such that the returned Faust is F = B * P where B is the Faust butterfly approximation of M*P.'.  If the array of indices is not a valid permutation the behaviour is undefined (however an invalid size or an out of bound index raise an exception).
%> 2. perm is a cell array of arrays of permutation column indices as defined in 1. In that case, all permutations passed to the function are used as explained in 1, each one producing a Faust, the best one (that is the best approximation of M in the Frobenius norm) is kept and returned by butterfly.
%> 3. perm is 'default_8', this is a particular case of 2. Eight default permutations are used. For the definition of those permutations please refer to [2].
%> 4. perm is 'bitrev': in that case the permutation is the bit-reversal permutation (cf. matfaust.tools.bitrev_perm).
%> 5. By default this argument is empty, no permutation is used (this is equivalent to using the identity permutation matrix in 1).
%>
%> @note Below is an example of how to create a permutation matrix from a permutation list
%> of indices (as defined by the perm argument) and conversely how to convert
%> a permutation matrix to a list of indices of permutation.
%> @code
%> >> I = randperm(8) % random permutation as a list indices
%> I =
%>
%>      6     1     5     3     7     2     8     4
%>
%> >> % convert a permutation as a list of indices to a permutation matrix P as a csr_matrix
%> >> n = numel(I);
%> >> P = sparse(I, 1:n, 1);
%> >> full(P)
%>
%> ans =
%>
%>      0     1     0     0     0     0     0     0
%>      0     0     0     0     0     1     0     0
%>      0     0     0     1     0     0     0     0
%>      0     0     0     0     0     0     0     1
%>      0     0     1     0     0     0     0     0
%>      1     0     0     0     0     0     0     0
%>      0     0     0     0     1     0     0     0
%>      0     0     0     0     0     0     1     0
%>
%> >> [I_, ~, ~] = find(P);
%> >> I_ = I_.'
%> I_ =
%>
%>      6     1     5     3     7     2     8     4
%> >> all(I_ == I)
%>
%> ans =
%>
%>   logical
%>
%>      1
%> @endcode
%>
%> @retval F the Faust which is an approximation of M according to a butterfly support.
%>
%> @b Examples:
%> @code
%> >> import matfaust.wht
%> >> import matfaust.dft
%> >> import matfaust.fact.butterfly
%> >> H = full(wht(32));
%> >> F = butterfly(H, 'type', 'bbtree');
%> >> err = norm(F-H)/norm(H)
%> err =
%>
%>    1.3947e-15
%> >> % it works with dft too!
%> >> % all you need is to specify the bit-reversal permutation
%> >> % since the Discrete Fourier Transform is the product of a butterfly factors with this particular permutation
%> >> DFT = full(dft(32));
%> >> F = butterfly(DFT, 'type', 'bbtree', 'perm', 'bitrev');
%> >> err = norm(full(F)-DFT)/norm(DFT)
%> err =
%> 	  2.9339e-15
%> @endcode
%>
%> Use butterfly with simple permutations:
%> @code
%> >> M = rand(4, 4);
%> >> % without any permutation
%> >> F1 = butterfly(M, 'type', 'bbtree');
%> >> % which is equivalent to using the identity permutation
%> >> p = 1:4;
%> >> F2 = butterfly(M, 'type', 'bbtree', 'perm', p);
%> >> % compute the relative diff
%> >> norm(F2-F1)/norm(F1)
%> ans =
%>         0
%> >> % then try another permutation
%> >> p2 = [2, 1, 4, 3];
%> >> F3 = butterfly(M, 'type', 'bbtree', 'perm', p2);
%> @endcode
%>
%> Use butterfly with a permutation defined by a list of indices J:
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
%>
%> >> % this is equivalent to passing a list containing a single permutation:
%> >> % F = butterfly(H, 'type', 'bbtree', 'perm', {J})
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
%> @endcode
%>
%>
%> Or to to use the 8 default permutations (keeping the best approximation resulting Faust)
%> @code
%> >> F = butterfly(H, 'type', 'bbtree', 'perm', 'default_8')
%> F =
%>
%> Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
%> - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> @endcode
%>
%>
%> <b>References:</b>
%> <br/> <b>[1]</b> Leon Zheng, Elisa Riccietti, and Remi Gribonval, <a href="https://arxiv.org/pdf/2110.01230.pdf">Hierarchical Identifiability in Multi-layer Sparse Matrix Factorization</a> <br/>
%> <b>[2]</b> T. Dao, A. Gu, M. Eichhorn, A. Rudra, and C. Re,
%> “Learning Fast Algorithms for Linear Transforms Using
%> Butterfly Factorizations,” in Proceedings of the 36th
%> International Conference on Machine Learning. June
%> 2019, pp. 1517–1527, PMLR
%>
%> @b See also: matfaust.wht, matfaust.dft, matfaust.rand_butterfly
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
									if(nargin < i+1 ||  ~ is_array_of_indices(varargin{i+1}, M) && ~ is_cell_arrays_of_indices(varargin{i+1}, M) && ~ strcmp(varargin{i+1}, 'default_8') && ~ strcmp(varargin{i+1}, 'bitrev'))
                                                error('keyword argument ''perm'' must be followed by ''default_8'', ''bitrev'', an array of permutation indices or a cell array of arrays of permutation indices')
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
                 [permutations{i}, ~, ~] = find(P);
                 permutations{i} = permutations{i}.'; % just for readibility in case of printing
            end
            F = matfaust.fact.butterfly(M, 'type', type, 'perm', permutations);
            return;
		elseif strcmp(perm, 'bitrev')
			P = bitrev_perm(size(M, 2));
			[perm, ~, ~] = find(P);
			perm = perm.';
			F = matfaust.fact.butterfly(M, 'type', type, 'perm', perm);
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
				err = norm(full(F)-M, 'fro')/nM;
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
