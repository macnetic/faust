%% Description hadamard_mat
%  Computation of tha Hadamard matrix and its "native" factorization
%  [H, Fact] = hadamard_mat(M) computes the Hadamard matrix H of size 
%  2^M*2^M and its factorization Fact.
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


function [H,Fact] = hadamard_mat(M)

bloc = (1/sqrt(2))*[1 1;1 -1];
matbase = bloc;
matbase=kron(speye(2^(M-1)),matbase);
n=size(matbase,1);

L=(1:n/2);
id_i=[2*L-1,2*L];
id_j=[L,L+n/2];
values=ones(n,1);

Perm = sparse(id_i,id_j,values,n,n);
same_fact=matbase*Perm;

Fact = cell(1,M);
for i=1:M
    Fact{i} = same_fact;
end
    H=dvp(Fact);
end

