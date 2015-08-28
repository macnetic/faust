%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef matlab_faust < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = matlab_faust(varargin)
            this.objectHandle = mexLoadFaust(varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            faust_mex('delete', this.objectHandle);
        end

        %% Multiplication faust-vector or faust-matrix
        function varargout = mtimes(this, varargin)
             [varargout{1:nargout}] = faust_mex('multiply', this.objectHandle, varargin{:});
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