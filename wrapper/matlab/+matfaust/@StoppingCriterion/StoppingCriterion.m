%% class StoppingCriterion
%%
classdef StoppingCriterion
	properties (SetAccess = public)
		is_criterion_error
		num_its
		error_treshold
		max_num_its
	end
	methods
		function stop_crit = StoppingCriterion(varargin)
			%set default values
			is_criterion_error = false;
			error_treshold = 0.3;
			num_its = 500;
			max_num_its = 1000;
			if(nargin < 1)
				error('matfaust.StoppingCriterion needs at least one argument.')
			else
				if(~ isscalar(varargin{1}))
						error('matfaust.StoppingCriterion 1th argument must be a scalar.')
				end
				if(nargin == 1 && (varargin{1}-floor(varargin{1})) == 0) % num_its criterion
					num_its = varargin{1};
					num_its = floor(num_its);
				else
					is_criterion_error = true; %varargin{1};
					%if( ~islogical(is_criterion_error))
					%	error('matfaust.StoppingCriterion 1th argument (is_criterion_error) must be logical.')
					%end
					error_treshold = varargin{1};
					if(~ isscalar(error_treshold) || ~ isreal(error_treshold))
						error('matfaust.StoppingCriterion 2nd argument (error_treshold) must be a real number.')
					end
					if(nargin > 1)
						max_num_its = varargin{2};
						if(~ isscalar(max_num_its) || ~ isreal(max_num_its))
							error('matfaust.StoppingCriterion 3rd argument (max_num_its) must be an integer.')
						end
						max_num_its = floor(max_num_its);
					end
				end
			end
			stop_crit.is_criterion_error = is_criterion_error;
			stop_crit.num_its = num_its;
			stop_crit.max_num_its = max_num_its;
			stop_crit.error_treshold = error_treshold;
		end
	end
end
