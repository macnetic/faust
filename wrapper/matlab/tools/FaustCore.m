%class FAUSTCORE
% this class only stored the Faust C++ object and 
% inherit from matlab handle class which allows different
% object to the same reference
% 
% WARNING : the user must not directly use this class,
% he must use the matlab class Faust
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://faust.gforge.inria.fr>

classdef FaustCore < handle
    properties (SetAccess = public, Hidden = false)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = FaustCore(varargin)
            % Constructor - build a faust from a cell array of matrix and a scalar (optional)
            %                1st input : 1D cell array of matrix (sparse or dense)
            %                2nd input : (optional) multiplicative scalar
	    %              - or from a filename (mat file) where a faust is stored with save_faust
	    if (nargin == 1) && ischar(varargin{1})
		filename=varargin{1};
		load(filename);
		if (~exist('faust_factors','var') || ~exist('transpose_flag','var'))
			error('Faust : invalid file');
		end
		this=FaustCore(faust_factors);
	    else				
		this.objectHandle = mexFaust('new',varargin{:});
	    end
		
	end

        %% Destructor - Destroy the C++ class instance
        function delete(this)
            % destructor delete the faust
            mexFaust('delete', this.objectHandle);
        end
        
        
        
    end
    
end












