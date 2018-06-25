%% Description normalizeCol
% A_normalize = normalizeCol(A,varargin)
% each column of the ouput matrix is normalized (norm(A_normalize(:,i))=1)
% - A can be a matrix or cell-array representing a matrix
% - an optional scalar multiplicator can added in input 
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot, Adrien Leman, Thomas Gautrais, 
%			Luc Le Magoarou, Remi Gribonval
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
%   Thomas Gautrais	: thomas.gautrais@inria.fr
%   Luc Le Magoarou	: luc.le-magoarou@inria.fr
%   Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%
function A_normalize = normalizeCol(A,varargin)

nb_input=length(varargin);

if (nb_input == 1)
    lambda=varargin{1};
else
    lambda=1;
end

if(iscell(A))
    X =lambda*dvp(A);
else
    X = A;
end

X(:,sum(X.^2,1)==0)=1;% null column are fixed to 1

if (iscell(A))
    A_normalize=A;
    A_normalize{end}=A_normalize{end}./repmat(sqrt(sum(X.^2,1)),size(X,1),1);
else
    A_normalize = A./repmat(sqrt(sum(X.^2,1)),size(X,1),1);
end

end

