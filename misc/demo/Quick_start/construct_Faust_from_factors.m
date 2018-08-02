%% Description construct_Faust_from_factor
%
%  This demo presents the method to build a Faust 
%  from its factors
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
%   Thomas Gautrais	: thomas.gautrais@inria.fr
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%
%%
import matfaust.Faust
% number of row of the Faust
dim1=300;
% number of column of the Faust
dim2=100;
% number of factor of the Faust
nb_factor=4;
% density of each factor
density_factor=0.1;

%cell-array representing its factors
factors = cell(1,nb_factor);

% 1st factor is a rectangular sparse matrix of density equal to density_factor
factors{1} = sprand(dim1,dim2,density_factor);

% all the others factor are square sparse matrix of density equal to density_factor
for i=2:nb_factor
  factors{i} = sprand(dim2,dim2,density_factor); 
end

%% construct the Faust
A_faust=Faust(factors);


% a multiplicative scalar can be taken into account to construct the Faust 
lambda=2;

% B_faust=lambda*A_faust
B_faust=Faust(factors,lambda);
