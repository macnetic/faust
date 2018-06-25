%% source http://www.cs.cornell.edu/~bindel/class/cs5220-s10/slides/FFT.pdf

function facts = TransStockhamFact(dim)
%  TransStockham Factorization of the Fourier Matrix.
%
%  facts = TransStockhamFact(dim) factorize the Fourier matrix of dimension dim
%
%                    facts is a cell-array of sparse matrix where the product represent the Fourier matrix
%                    
%
%                    
%                    warning : dim must be a power of 2
%                               

%
%% source http://www.cs.cornell.edu/~bindel/class/cs5220-s10/slides/FFT.pdf
% 
%
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot,  Luc Le Magoarou, Remi Gribonval
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
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
%	approximations of matrices and applications", Journal of Selected
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
nb_fact = 2*log2(dim);

facts=cell(nb_fact,1);

for q=1:nb_fact/2
    L=2^q;
    r=dim/L;
    L_star = L/2;
    r_star = 2 * r;
    wL=exp(-2*i*pi/L);
    index=1:L_star;
    OmegaDiag=wL.^(index-1);
    Omega=sparse(index,index,OmegaDiag,L/2,L/2);
    BL=[speye(L_star) Omega;speye(L_star) -Omega];
    Aq = kron(speye(r),BL);
    
    Id_r_star = speye(r_star);
    Pi_r_star=Id_r_star(:,[(1:2:r_star),(2:2:r_star)]);
    Gq = kron(Pi_r_star,speye(L_star));
    %%nb_fact+2-2*i
    %%nb_fact+2-2*i-1
    facts{nb_fact+2-2*q}=Gq;
    facts{nb_fact+2-2*q-1}=Aq;
end


end

