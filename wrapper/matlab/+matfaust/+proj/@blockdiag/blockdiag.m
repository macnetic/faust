%==================================================
%> @brief Functor for the BLOCKDIAG projector.
%>
%> TODO
%==================================================
 classdef blockdiag < matfaust.proj.proj_gen
	properties
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
			for i=1:length(mn_cell)
				m_vec(i) = mn_cell{i}{1};
				n_vec(i) = mn_cell{i}{2};
				mn_mat(i,1) = mn_cell{i}{1};
				mn_mat(i,2) = mn_cell{i}{2};
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
					pM = mexFaustReal('prox_blockdiag', M, self.m_vec, self.n_vec, self.normalized, self.pos);
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
