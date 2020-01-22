%> @package matfaust.proj @brief This module provides matrix projectors.

%=========================================================
%> @brief The parent abstract class to represent projectors.
%=========================================================
classdef (Abstract) proj_gen
	properties (SetAccess = protected)
		constraint
	end
	methods
		function pM = subsref(self, S)
			if(~ strcmp(S.type, '()') && ~ strcmp(S.type, '.'))
				error('Invalid use of projector functor object: only () or . are handled')
			end
			if(iscell(S.subs) && length(S.subs) == 1 && ismatrix(S.subs{1}))
					%error('The projector must be called on a matrix')
				M = S.subs{1};
				pM = self.constraint.project(M);
			elseif(ischar(S.subs) && strcmp('constraint', S.subs))
				pM = self.constraint;
			else
				error('bad use of projector: must be projector(matrix) or projector.constraint.')
			end
		end
	end
end
