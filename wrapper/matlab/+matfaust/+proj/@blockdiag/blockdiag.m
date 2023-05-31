% ==================================================
%> @brief Functor for the BLOCKDIAG projector.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b blockdiag(shape, block_shapes): returns a BLOCKDIAG projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), block_shapes defines the diagonal blocks shapes.<br/>
%>
%> @param shape: vector of size 2, to define the size of the input matrix.
%> @param block_shapes: a cell array in the defined as follows in order to get n diagonal blocks: {{b1_nrows, b1_ncols}, …, {bn_nrows, bn_ncols}}. The sum of bi_nrows must equal shape(1), likewise the sum of bi_ncols must equal shape(2).
%>
%> @retval blockdiag projector
%>
%> @b Example
%> @code
%> >> import matfaust.proj.blockdiag
%> >> rng(42)
%> >> M = rand(5, 5)
%>
%> M =
%>
%>     0.3745    0.1560    0.0206    0.1834    0.6119
%>     0.9507    0.0581    0.9699    0.3042    0.1395
%>     0.7320    0.8662    0.8324    0.5248    0.2921
%>     0.5987    0.6011    0.2123    0.4319    0.3664
%>     0.1560    0.7081    0.1818    0.2912    0.4561
%>
%> >> p = blockdiag(size(M), {{1, 1}, {3, 3}, {1, 1}}, 'normalized', false);
%> >> p(M)
%>
%> ans =
%>
%>     0.1948         0         0         0         0
%>          0    0.0302    0.5045    0.1582         0
%>          0    0.4505    0.4330    0.2729         0
%>          0    0.3127    0.1104    0.2247         0
%>          0         0         0         0    0.2372
%> >>
%> @endcode
% ==================================================
 classdef blockdiag < matfaust.proj.proj_gen
	properties(Access = private)
		m_vec
		n_vec
		normalized
		pos
	end
	methods

		function proj = blockdiag(shape, mn_cell, varargin)
			M = zeros(shape(1), shape(2));
			m_vec = zeros(1, length(mn_cell));
			n_vec = zeros(1, length(mn_cell));
			mn_mat = zeros(length(mn_cell), 2);
			m_vec(1) = mn_cell{1}{1};
			n_vec(1) = mn_cell{1}{2};
			for i=2:length(mn_cell)
				m_vec(i) = mn_cell{i}{1}+m_vec(i-1);
				n_vec(i) = mn_cell{i}{2}+n_vec(i-1);
				mn_mat(i,1) = m_vec(i);
				mn_mat(i,2) = n_vec(i);
			end
			proj.m_vec = m_vec;
			proj.n_vec = n_vec;
			proj.normalized = false;
			proj.pos = false;
			proj.constraint = matfaust.factparams.ConstraintMat('blockdiag', mn_mat, varargin{:});
			argc = length(varargin);
			if(argc > 0)
				for i=1:argc
					switch(varargin{i})
						case 'normalized'
							if(argc == i || ~ islogical(varargin{i+1}))
								error('normalized keyword arg. is not followed by a boolean')
							end
							proj.normalized = varargin{i+1};
						case 'pos'
							if(argc == i || ~ islogical(varargin{i+1}))
								error('pos keyword arg. is not followed by a boolean')
							end
							proj.pos = varargin{i+1};
						otherwise
							if(isstr(varargin{i}))
								error([ varargin{i} ' unrecognized argument'])
							end
					end
				end
			end

		end

		function pM = subsref(self, S)
			if(~ strcmp(S.type, '()') && ~ strcmp(S.type, '.'))
				error('Invalid use of projector functor object: only () or . are handled')
			end
			if(iscell(S.subs) && length(S.subs) >= 1 && ismatrix(S.subs{1}))
					%error('The projector must be called on a matrix')
				M = S.subs{1};
				% pM = self.constraint.project(M);
				if(isreal(M))
					if(strcmp(class(M), 'single'))
						pM = mexFaustRealFloat('prox_blockdiag', M, self.m_vec, self.n_vec, self.normalized, self.pos);
					else
						pM = mexFaustReal('prox_blockdiag', M, self.m_vec, self.n_vec, self.normalized, self.pos);
					end
				else
					pM = mexFaustCplx('prox_blockdiag', M, self.m_vec, self.n_vec, self.normalized, self.pos);
				end
			elseif(ischar(S.subs) && strcmp('constraint', S.subs))
				pM = self.constraint;
			else
				error('wrong use of projector: must be projector(matrix) or projector.constraint.')
			end
		end

	end
end
