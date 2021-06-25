classdef ParamsFactFactory
	properties(Constant)
		SIMPLIFIED_PARAMS = {
		{'squaremat', 'hadamard'},
		{'rectmat', 'meg'},
		{'dft'}
		};
		SQRTMAT_ID = 1;
		RECTMAT_ID = 2;
		DFTMAT_ID = 3;
	end
	methods(Static)
		function sp = createParams(M, p)
			import matfaust.factparams.*
			import matfaust.factparams.ParamsFactFactory.*
			param_id = -1;
			if(~ is_a_valid_simplification(p))
				error('Invalid p to represent a simplified parametrization')
			end
			if(ischar(p))
				param_id = lower(p);
			else % p{1} is a cell with a char array at first pos
				param_id = lower(p{1});
			end
			if(strcmp_anycell(param_id, ParamsFactFactory.SIMPLIFIED_PARAMS{ParamsFactFactory.SQRTMAT_ID}))
				sp = ParamsHierarchicalSquareMat.createParams(M, p);
			elseif(strcmp_anycell(param_id, ParamsFactFactory.SIMPLIFIED_PARAMS{ParamsFactFactory.RECTMAT_ID}))
				sp = ParamsHierarchicalRectMat.createParams(M, p);
			elseif(strcmp_anycell(param_id, ParamsFactFactory.SIMPLIFIED_PARAMS{ParamsFactFactory.DFTMAT_ID}))
				sp = ParamsHierarchicalDFT.createParams(M, p);
			else
				error('p is not a known simplified parametrization')
			end
		end
		function bool = is_a_valid_simplification(p)
			bool = ischar(p) || iscell(p) && length(p) > 0 && ischar(p{1});
		end
	end
	methods(Access = private, Static)
		function bool = strcmp_anycell(str, c)
			bool = false;
			if(~ iscell(c))
				bool = false;
			else
				for i=1:length(c)
					if(strcmp(c{i}, str))
						bool = true;
						return
					end
				end
			end
		end
	end
end

