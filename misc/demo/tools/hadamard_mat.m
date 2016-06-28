%% Description hadamard_mat
%  Computation of tha Hadamard matrix and its "native" factorization
%  [H, Fact] = hadamard_mat(n) computes the Hadamard matrix H of size 
%  2^n*2^n and its factorization Fact.
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
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


function [H, Fact] = hadamard_mat(n)

bloc = [1 1;1 -1];
matbase = bloc;
for i = 1:2^(n-1)-1
    matbase = blkdiag(matbase,bloc);
end
matbase = (1/sqrt(2))*matbase;

Perm = zeros(size(matbase));
for l = 1:size(Perm,1)/2
    Perm(2*l-1,l)=1;
    Perm(2*l,l+size(Perm,1)/2)=1;
end

Fact = zeros(2^n,2^n,n);
for i=1:n
    Fact(:,:,i) = matbase*Perm;
end
H = Fact(:,:,1)^n;
end