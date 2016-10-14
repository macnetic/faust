%% Description generate_params
% Generates parameter of hierarchical_fact from a Matrix, the number of
% factor of the factorization and its RCG (Relative Complexity Gain)
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
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
function params = generate_params(dim1,dim2,nb_factor,RCG)

params.nrow=dim1;
params.ncol=dim2;

params.nfacts=nb_factor;



[min_dim,id_min_dim]=min([dim1 dim2]);


if (id_min_dim == 1)
    params.fact_side=1; % left-side is factorized iteratively
    dim1_residu=dim1;
    dim2_residu=min_dim;
    id_factor=nb_factor:-1:1;
    line_residu=1;
    line_fact=2;
else
    params.fact_side=0; % right-side is factorized iteratively
    dim1_residu=min_dim;
    dim2_residu=dim2;
    id_factor=1:nb_factor;
    line_residu=2;
    line_fact=1;
end




averaging_coeff=zeros(1,nb_factor);

dims_factor=zeros(2,nb_factor);% list of the dimension of the factors :
                               % -dim_factor(1,i) is the number of row of
                               % the ith factor
                               % -dim_factor(2,i) is the number of columns
                               % of the ith factor

dims_factor(:,1)=[dim1;min_dim];
for i=2:nb_factor-1
   dims_factor(:,i)=[min_dim;min_dim];  
end
dims_factor(:,nb_factor)=[min_dim;dim2];

nb_coeff_per_fact=prod(dims_factor);
total_nb_coeff=sum(nb_coeff_per_fact);
density_fact=(dim1*dim2)/(total_nb_coeff*RCG);% if Matrix is squared, alpha = 1/(nb_factor*RCG)



% each factor has the same density equal to 1/(RCG*nb_factor)
factor_constraint=cell(1,nb_factor);



% each factor has the same density
for i=id_factor
   % constraint,nb_nnz,nb_row,nb_col of the factor
   dim1_factor=dims_factor(1,i);
   dim2_factor=dims_factor(2,i);
   factor_constraint{i}={'sp',density_fact*nb_coeff_per_fact(i),dim1_factor,dim2_factor};
end



decrease_speed=exp(log(density_fact)/(nb_factor-1));


params.cons=cell(2,nb_factor-1);
params.cons(line_fact,1:nb_factor-1)=factor_constraint(id_factor(1:end-1));
%params.cons(line_fact,nb_factor-1)=factor_constraint(id_factor(end)); %last residuum is a factor




for i=1:nb_factor-1
   params.cons{line_residu,i}={'sp',decrease_speed^i*dim1_residu*dim2_residu,dim1_residu,dim2_residu}; 
end





end

