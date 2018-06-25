%% Description grad_comp
%  Computation of the gradient and Lipschitz modulus
%  [grad, LC] = grad_comp(L,S,R,X,lambda) computes the gradient grad of
%  H(L,S,R,lambda) = || X - lambda*L*S*R || and its Lipschitz modulus LC.
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


function [grad, LC] = grad_comp(L,S,R,X,lambda)

% Compute the gradient
gradtemp = lambda*mult_left(L,S);
gradtemp = mult_right(gradtemp,R);
gradtemp = gradtemp-X;
gradtemp = lambda*mult_right(gradtemp',L);
gradtemp = mult_left(R,gradtemp);
grad = gradtemp';

% Compute the Lipschitz constant
LC = lambda^2 * norm(dvp(R))^2* norm(dvp(L))^2;
%LC=1;
end
