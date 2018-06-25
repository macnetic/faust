%% Description gen_artficial_faust
% facts=gen_artficial_faust(Dim,RCG,nb_fact,constraint)
% create a random Faust factor with :
% a given dimension Dim,
% a relative complexity gain RCG, 
% a number of factor nb_fact, 
% a type of sparsity 'sp' random support, 
%                    'sp_row', fixed number of non-zeros per row
%                    'sp_col', fixed number of nonzeros-per-col   
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

function facts=gen_artficial_faust(Dim,RCG,nb_fact,constraint)
%% cas 1
% RCGs = [2 4 6 8 10];
% Dims = [128 256 512];
% nb_facts = [2,4,6];

%% cas 2

% RCGs = [2 4 8 16];
% Dims = [128 256 512 1024 2048];

% nb_facts = [2,4,8,16];

% RCGs = [2 4];
% Dims = [32 64 128];
% nb_facts = [2,4];


% constraint='sp_row';% choix possible 'sp_row' et 'sp_col_'











if ((strcmp(constraint,'sp') + strcmp(constraint,'sp_row') + strcmp(constraint,'sp_col')) == 0)
    error('the constraint parameter must be equal to sp or sp_row or sp_col');
end


dim1=Dim;
dim2=Dim;
densite_per_fact = 1/(nb_fact*RCG);


facts=cell(1,nb_fact);


if (strcmp(constraint,'sp_row'))
    nb_elt_per_row = round(dim2*densite_per_fact);
    nb_elt = nb_elt_per_row * dim1;
    id_i=reshape(repmat((1:dim1),nb_elt_per_row,1),dim1*nb_elt_per_row,1);
    id_j=zeros(nb_elt,1);
    value=zeros(nb_elt,1);
    
end
if (strcmp(constraint,'sp_col'))
    nb_elt_per_col = round(dim1*densite_per_fact);
    nb_elt = nb_elt_per_col * dim2;
    id_j=reshape(repmat((1:dim2),nb_elt_per_col,1),dim2*nb_elt_per_col,1);
    id_i=zeros(nb_elt,1);
    value=zeros(nb_elt,1);
end
for k=1:nb_fact
    if (strcmp(constraint,'sp'))
        facts{k} = sprand(dim1,dim2,densite_per_fact);
    end
    if (strcmp(constraint,'sp_row'))
        value=rand(nb_elt,1);
        for ll=0:dim1-1
            id_j(ll*nb_elt_per_row+1:(ll+1)*nb_elt_per_row)=randperm(dim2,nb_elt_per_row)';
        end
        facts{k}=sparse(id_i,id_j,value,dim1,dim2);
    end
    if (strcmp(constraint,'sp_col'))
        value=rand(nb_elt,1);
        for ll=0:dim2-1
            id_i(ll*nb_elt_per_col+1:(ll+1)*nb_elt_per_col)=randperm(dim1,nb_elt_per_col)';
        end
        facts{k}=sparse(id_i,id_j,value,dim1,dim2);
    end
    
    
    
    
   
    
    
    
    
    
    
end


end
