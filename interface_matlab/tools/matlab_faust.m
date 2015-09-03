%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef matlab_faust < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = matlab_faust(varargin)
            this.objectHandle = faust_mex('new',varargin{:});
			%this.objectHandle = mexLoadFaust(varargin{:});
		end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            faust_mex('delete', this.objectHandle);
        end

        %% Multiplication faust-vector or faust-matrix
        function varargout = mtimes(this, varargin)
             [varargout{1:nargout}] = faust_mex('multiply', this.objectHandle, varargin{:});
         end
		 
		 %% Evaluate the product of a faust_core
		 function varargout = get_product(this)
				[varargout{1:nargout}]=faust_mex('get_product',this.objectHandle);
	     end
		 
		 function trans=transpose(this)
			if (nargout	 == 0)
				faust_mex('transpose',this.objectHandle);
			else
				trans = matlab_faust({});
				trans.objectHandle = faust_mex('transpose',this.objectHandle);
				
			end
				
		 end
		 
		function Size=size(this,varargin);
			if (nargin == 1)
				Size=faust_mex('size',this.objectHandle);
			else (nargin == 2)
				Size=faust_mex('size',this.objectHandle,varargin);
			end
		end
        
%         %% Train - an example class method call
%         function varargout = train(this, varargin)
%             [varargout{1:nargout}] = faust_mex('train', this.objectHandle, varargin{:});
%         end
% 
%         %% Test - another example class method call
%         function varargout = test(this, varargin)
%             [varargout{1:nargout}] = faust_mex('test', this.objectHandle, varargin{:});
%         end
    end
end