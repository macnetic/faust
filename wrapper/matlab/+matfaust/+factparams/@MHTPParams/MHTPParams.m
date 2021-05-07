% experimental block start
% =========================================================
%> @brief This class defines the set of parameters to run the MHTP-PAL4MSA algorithm.
%>
%> See also matfaust.fact.palm4msa_mhtp, matfaust.fact.hierarchical_mhtp.
% =========================================================
classdef MHTPParams
	properties (SetAccess = public, Hidden = false)
		num_its;
		constant_step_size;
		step_size;
		palm4msa_period;
		updating_lambda;
	end
	methods
		% =========================================================
		%> Constructor of the MHTPParams class.
		%>
		%> See also matfaust.fact.palm4msa_mhtp, matfaust.fact.hierarchical_mhtp.
		%>
		%> @param 'num_its', int: (optional) the number of iterations to run the MHTP algorithm.
		%> @param 'constant_step_size', bool: (optional) true to use a constant step for the gradient descent, False otherwise. If false the step size is computed dynamically along the iterations (according to a Lipschitz criterion).
		%> @param 'step_size', real: (optional) The step size used when constant_step_size==true.
		%> @param 'palm4msa_period', int: (optional) The period (in term of iterations) according to the MHTP algorithm is ran (i.e.: 0 <= i < N being the PALM4MSA iteration, MHTP is launched every i = 0 (mod palm4msa_period). Hence the algorithm is ran one time at least – at PALM4MSA iteration 0).
		%> @param 'updating_lambda', bool: (optional) if true then the scale factor of the Faust resulting of the factorization is updated after each iteration of MHTP (otherwise it never changes during the whole MHTP execution).
		% =========================================================
 		function p = MHTPParams(varargin)
			argc = length(varargin);
			% default parameter values
			p.num_its = 50;
			p.constant_step_size = false;
			p.step_size = 1e-3;
			p.palm4msa_period = 1000;
			p.updating_lambda = true;
			if(argc > 0)
				for i=1:2:argc
					if(argc > i)
						% next arg (value corresponding to the key varargin{i})
						tmparg = varargin{i+1};
					end
					switch(varargin{i})
						case 'step_size'
							if(argc == i || ~ isscalar(tmparg) || ~ isreal(tmparg) || tmparg < 0)
								error('step_size argument must be followed by a positive real value')
							else
								p.step_size = tmparg;
							end
						case 'num_its'
							if(argc == i || ~ isscalar(tmparg) || ~isreal(tmparg) || tmparg < 0 || tmparg-floor(tmparg) > 0)
								error('num_its argument must be followed by an integer')
							else
								p.num_its = tmparg;
							end
						case 'palm4msa_period'
							if(argc == i || ~ isscalar(tmparg) || ~isreal(tmparg) || tmparg < 0 || tmparg-floor(tmparg) > 0)
								error('palm4msa_period argument must be followed by an integer')
							else
								p.palm4msa_period = tmparg;
							end
						case 'updating_lambda'
							if(argc == i || ~ islogical(tmparg))
								error('updating_lambda argument must be followed by a logical')
							else
								p.updating_lambda = tmparg;
							end
						case 'constant_step_size'
							if(argc == i || ~ islogical(tmparg))
								error('constant_step_size argument must be followed by a logical')
							else
								p.constant_step_size = tmparg;
							end
						otherwise
							if((isstr(varargin{i}) || ischar(varargin{i})))
								error([ tmparg ' unrecognized argument'])
							end
						end
					end
				end
		end
	end
end
% experimental block end
