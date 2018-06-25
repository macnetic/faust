%% Description prox_spcol 
%  Projection onto the set of matrices with sparse columns and 
%  unit Frobenius norm.
%  Xprox = prox_spcol(X,s) projects the input matrix X onto the set of
%  matrices which have at most s non-zero entries per column and unit 
%  Frobenius norm. Xprox is the projection of X onto this set.
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Luc Le Magoarou, Remi Gribonval
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
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%


function Xprox = prox_spcol(X,s)

Xabs = abs(X);
Xprox = zeros(size(X));

[~,sortIndex] = sort(Xabs,'descend');
maxIndex = sortIndex(1:round(s),:);
incre = 0:size(X,1):size(X,1)*(size(X,2)-1);
incremat = repmat(incre,s,1);
maxIndex = maxIndex + incremat;
Xprox(maxIndex(:)) = X(maxIndex(:));
Xprox = Xprox/norm(Xprox,'fro');
end