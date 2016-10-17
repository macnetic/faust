%% class FAUSTCORE
% this class only stored the Faust C++ object and 
% inherit from matlab handle class which allows different
% object to have the same reference
% 
% WARNING : the user must not directly use this class,
% he must use the matlab class Faust
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot, Adrien Leman, Thomas Gautrais, Luc Le Magoarou, Remi Gribonval
%			INRIA Rennes, FRANCE
%			http://www.inria.fr/
%
% The FAuST Toolbox is distributed under the terms of the GNU Affero 
% General Public License.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published 
% by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% Contacts:
%   Nicolas Bellot	: nicolas.bellot@inria.fr
%   Adrien Leman	: adrien.leman@inria.fr
%   Thomas Gautrais : thomas.gautrais@inria.fr
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%

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
		if (~exist('Faust_factors','var') )
			error('FaustCore : invalid file');
		end
		this=FaustCore(Faust_factors);
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












