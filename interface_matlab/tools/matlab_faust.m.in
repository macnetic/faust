%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef matlab_faust < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = matlab_faust(varargin)
            this.objectHandle = mexFaust('new',varargin{:});
			%this.objectHandle = mexLoadFaust(varargin{:});
		end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            mexFaust('delete', this.objectHandle);
        end

        %% Multiplication faust-vector or faust-matrix
        function varargout = mtimes(this, varargin)
             [varargout{1:nargout}] = mexFaust('multiply', this.objectHandle, varargin{:});
         end
		 
		 %% Evaluate the product of a faust_core
		 function varargout = get_product(this)
				[varargout{1:nargout}]=mexFaust('get_product',this.objectHandle);
	     end
		 
		 function trans=transpose(this)
			if (nargout	 == 0)
				mexFaust('transpose',this.objectHandle);
			else
				trans = matlab_faust({});
				trans.objectHandle = mexFaust('transpose',this.objectHandle);
				
			end
				
		 end
		function res=solve(this,varargin)
			[varargout{1:nargout}] = mexFaust('solve', this.objectHandle, varargin{:});
		end
		 
		function Size=size(this,varargin);
			if (nargin == 1)
				Size=mexFaust('size',this.objectHandle);
			else (nargin == 2)
				Size=mexFaust('size',this.objectHandle,varargin);
			end
		end
        
    end
end