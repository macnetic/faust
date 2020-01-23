%==================================================
%> @brief A helper class for constructing a list of consistent matfaust.proj.proj_gen projectors or ConstraintGeneric objects.
%>
%> NOTE: ConstraintGeneric use is not advised (these objects are not well documented). Use rather the projectors functors (from matfaust.proj namespace).
%==================================================
classdef ConstraintList
	properties(SetAccess = public)
		clist
	end
	methods
		function this = ConstraintList(varargin)
			import matfaust.factparams.*
			tuple_len = 4; % name, value, nrows, ncols
			i = 1;
			j = 1; % the number of processed constraints
			nargs = length(varargin);
			this.clist = {};
			while(i <= length(varargin))
				if(isa(varargin{i}, 'ConstraintGeneric'))
					this.clist = [ this.clist, {varargin{i}} ];
					i = i + 1;
					continue
				end
				cname = ConstraintName(varargin{i});
				if(i+1 > nargs)
					% ENOTE: throw() is a function (it needs ())
					throw(MException('ConstraintList:ErroneousVarargin',...
					['No value/parameter given to define the ', int2str(j), '-th constraint.']))
				end
				cval = varargin{i+1};
				if(i+2 > nargs)
					throw MException('ConstraintList:ErroneousVarargin',...
					['No number of rows given to define the ', int2str(j), '-th constraint.'])
				end
				nrows = varargin{i+2};
				if(i+3 > nargs)
					throw(MException('ConstraintList:ErroneousVarargin',...
					['No number of columns given to define the ', int2str(j), '-th constraint.']))
				end
				ncols = varargin{i+3};
				if(cname.is_int_constraint())
					cons = ConstraintInt(cname, nrows, ncols, cval);
				elseif(cname.is_real_constraint())
					cons = ConstraintReal(cname, nrows, ncols, cval);
				elseif(cname.is_mat_constraint())
					cons = ConstraintMat(cname, cval);
				else
					% it shouldn't happen because ConstraintName has verified the input data
					msg = ' Not a valid name for a ConstraintGeneric object.';
					throw(MExecption('ConstraintList:ErroneousVarargin', msg))
				end
				this.clist = [ this.clist, {cons} ];
				i = i + tuple_len;
			end
		end
	end
end
