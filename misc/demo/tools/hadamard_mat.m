%% Description speed_up_hadamard
%
%  This demo makes some time comparison between (Hadamard matrix)-vector multiplication and
%  (Hadamard factorisation i.e a FAµST)-vector multiplication for different dimension
%  of the Hadamard matrix.
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
%   Leman Adrien   : adrien.leman@inria.fr
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
%	approximations of matrices and applications", Journal of Selected
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%

nb_mult=500;
Ms=6:11;
ns=2.^Ms;
nb_dim=length(Ms);
threshold=10^(-10);
dense_times=zeros(nb_mult,nb_dim);
faust_times=zeros(nb_mult,nb_dim);

h = waitbar(0,'speed up hadamard : Generation of the data ...');
Hadamard_matrices=cell(1,nb_dim);
Hadamard_facts=cell(1,nb_dim);
for k=1:nb_dim
    waitbar(k/nb_dim);
    M=Ms(k);
    n=ns(k);
    % generation of the hadamard factorisation
    [H,facts] = hadamard_mat(M);
    Hadamard_matrices{k}=H;
    Hadamard_facts{k}=facts;
end
close(h);

%     figure,
%         for i=1:M
%             subplot(2,M,i);
%             imagesc(facts{i});
%             axis image
%             subplot(2,M,M+i);
%             imagesc(cum_Hbis{i});
%             axis image
%         end
hadamard_faust=matlab_faust(facts);
hadamard_dense=dvp(facts);

h = waitbar(0,'speed up hadamard : multiplication time comparison ...');
for i=1:nb_mult
    waitbar(i/nb_mult);
    for k=1:nb_dim
        n=ns(k);
        hadamard_dense=full(Hadamard_matrices{k});
        hadamard_faust=matlab_faust(Hadamard_facts{k});
        
        x=rand(n,1);
        ydense=zeros(n,1);
        yfaust=zeros(n,1);
        
        tic
        ydense=hadamard_dense*x;
        t1=toc;
        
        tic
        yfaust=hadamard_faust*x;
        t2=toc;
        
        if(norm(ydense-yfaust)>threshold)
            error('speed_up hadamard : multiplication problem');
        end
        
        
        dense_times(i,k)=t1;
        faust_times(i,k)=t2;
    end
end
close(h);

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
    H=full(dvp(Fact));
end







