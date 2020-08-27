%================================
%> @brief This function loads explicitely the gpu_mod library in memory.
%>
%>    Normally it's not required to load the library manually, but it could be
%>    useful to set a non-default path or to diagnose a loading issue.
%>
%>    @param libpath the absolute or relative path where to find the dynamic
%>    library (gm) to load. By default, it's none to auto-find the library
%>    (if possible).
%>    @param backend the GPU backend to use, only cuda is available for now.
%>    @param silent if True nothing or almost will be displayed on loading
%>    (e.g. silent errors), otherwise all messages are visible.
%>
%===============================
function enable_gpu_mod(varargin)
	% default values
	backend = 'cuda';
	silent = true;
	osstr = computer;
	switch(osstr)
		case 'PCWIN64'
			libpath = 'C:\Program Files\Faust\lib\gm.dll';
		case 'GLNXA64'
			libpath = '/opt/local/faust/lib/libgm.so';
		case 'MACI64'
			libpath = '/opt/local/faust/lib/libgm.dylib';
		otherwise
			% not reachable
	end
	if(nargin > 0)
		i = 1;
		while (i <= nargin)
			switch(varargin{i})
				case 'libpath'
					libpath = varargin{i+1};
					if(nargin == i || ~ isstr(libpath))
						error('libpath must be a char array or string')
					end
				case 'backend'
					if(nargin == i || ~ isstr(libpath))
						error('libpath must be a char array or string')
					end
					if(~ strcmp(varargin{i+1}, 'cuda'))
						error('For now the only supported backend is cuda')
					end
				case 'silent'
					if(nargin == i || ~ islogical(varargin{i+1}))
						error('silent keyword argument is not followed by a logical')
					else
						silent = varargin{i+1};
					end
				otherwise
					if(isstr(varargin{i}) && (~ strcmp(varargin{i}, 'libpath') && ~ strcmp(varargin{i}, 'backend') && ~ strcmp(varargin{i}, 'silent')) )
						error([ varargin{i} ' unrecognized argument'])
					end
			end
			i = i + 2;
		end
	end
	mexFaustReal('enable_gpu_mod', libpath, backend, silent); % it should load the library for complex Faust too
end
