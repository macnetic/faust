%% class FAUSTCORE
% this class only stored the Faust C++ object and 
% inherit from matlab handle class which allows different
% object to have the same reference
% 
% WARNING : the user must not directly use this class,
% he must use the matlab class Faust
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.inria.fr>
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
		isRealFlag;
		dev;
		dtype;
	end
	methods
		%% Constructor - Create a new C++ class instance
		function this = FaustCore(varargin)
			if(nargin == 4) && ~iscell(varargin{1}) %&& isa(varargin{1}, 'handle'))
				%if(~ isvalid(varargin{1}))
				%    error('FaustCore: invalid handle to copy passed to the constructor.')
				%end
				this.objectHandle = varargin{1};
				if ((varargin{2} ~= 1) && (varargin{2} ~= 0))
					error('FaustCore: invalid argument 2 (isReal) passed to the constructor, must be equal to 0 or 1');
				end
				this.isRealFlag = varargin{2};
				this.dev = varargin{3};
				this.dtype = varargin{4};
			elseif(nargin >= 1)
				factors = varargin{1};
				scale = varargin{2};
				optCopy = varargin{3};
				isRealFlag = varargin{4};
				this.dev = varargin{5};
				onGPU = startsWith(this.dev, 'gpu');
				this.dtype = 'double';
				if(nargin >= 6)
					this.dtype = varargin{6};
				end
				% varargin{2} is lambda (optional)
				% varargin{3} is optimizedCopy boolean (not documented).
				if(onGPU)
					if (isRealFlag)
						if(strcmp(this.dtype, 'double'))
							this.objectHandle = mexFaustGPUReal('new',varargin{1:3});
						else
							this.objectHandle = mexFaustGPURealFloat('new',varargin{1:3});
						end
					else
						this.objectHandle = mexFaustGPUCplx('new',varargin{1:3});
					end
				else
					if (isRealFlag)
						if(strcmp(this.dtype, 'double'))
							this.objectHandle = mexFaustReal('new',varargin{1:3});
						else % float
							this.objectHandle = mexFaustRealFloat('new',varargin{1:3});
						end
					else
						this.objectHandle = mexFaustCplx('new',varargin{1:3});
					end
				end
				this.isRealFlag = isRealFlag;
			end
		end

		%% Destructor - Destroy the C++ class instance
		function delete(this)
			% destructor delete the faust
			if(isa(this.objectHandle, 'integer'))
				if(startsWith(this.dev, 'gpu'))
					if (this.isRealFlag)
						if(strcmp(this.dtype, 'double'))
							mexFaustGPUReal('delete', this.objectHandle);
						else % float
							mexFaustGPURealFloat('delete', this.objectHandle);
						end
					else
						mexFaustGPUCplx('delete', this.objectHandle);
					end

				else
					if (this.isRealFlag)
						if(strcmp(this.dtype, 'double'))
							mexFaustReal('delete', this.objectHandle);
						else % float
							mexFaustRealFloat('delete', this.objectHandle);
						end
					else
						mexFaustCplx('delete', this.objectHandle);
					end
				end
			end
		end

		function set_FM_mul_mode(this, mode)
			if(isa(this.objectHandle, 'integer'))
				if (this.isRealFlag)
					if(strcmp(this.dtype, 'double'))
						mexFaustReal('set_FM_mul_mode', this.objectHandle, mode);
					else % float
						mexFaustRealFloat('set_FM_mul_mode', this.objectHandle, mode);
					end
				else
					mexFaustCplx('set_FM_mul_mode', this.objectHandle, mode);
				end
			end
		end

	end

end












