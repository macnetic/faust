%==================================================
%> Functor for the BLOCKDIAG projector.
%==================================================
 classdef blockdiag %< matfaust.proj.proj_gen
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
			for i=1:length(mn_cell)
				m_vec(i) = mn_cell{i}{1};
				n_vec(i) = mn_cell{i}{2};
			end
			proj.m_vec = m_vec;
			proj.n_vec = n_vec;
			proj.normalized = false;
			proj.pos = false;
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
				error('bad use of projector: must be projector(matrix) or projector.constraint.')
			end
		end

	end
end
