%% Description mexFaust
% this mexfunction should not be called directly by the user,
% all the call of this mexfunction is done by the Faust.m  class
%
%  mexfunction handling all the method, function of the Faust class
%  arargout{1:nargout} = mexFaust(CMD,vargin{:}),
% for instance, 1st input CMD  can be equal to :
% 'new'-> create a new faust
%'delete'-> delete a faust
%'multiply'-> multiply a faust by a matrix
%'full'-> evaluate the product of a faust
%'transpose' -> transpose a faust
%'size' -> get the dimension of a faust
%
% for more information about the Faust class, 
%type 'doc Faust'

%% License:
% Copyright (2016):	Nicolas Bellot, Adrien Leman, Luc Le Magoarou, Remi Gribonval
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
%   Nicolas Bellot : nicolas.bellot@inria.fr
%   Adrien Leman   : adrien.leman@inria.fr
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%
