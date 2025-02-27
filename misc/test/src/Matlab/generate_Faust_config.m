%% Description generate_Faust_config
% This functions generate different list of factor representing a Faust
% and a dense matrix representing this Faust
%
% ouput parameters : -factors is a cell-array representing the factors of the Faust
%                      factors{1} is the leftmost factor, factors{end} is the rightmost factor
%                    - dense is the full-storage format matrix representing the Faust
%                       dense = factors{1} * ... * factors{end} 
%
% input parameters : - list_dim is a tab storing the different size of the factor
%		       list_dim(1) is the number of row of the 1st factor
%		       list_dim(end) is the number of column of the last factor
%                      otherwise list_dim(i) is the number of column ot the (ith -1) factor
%                                and the number of row of the ith factor   
% 	
%                    - type_factor cell-array specifying the type of Faust factors, sparse or dense
%                   for instance 
%		    if type_factor{1}='sparse'
%                   -	the first output factor of the Faust will be Sparse
%		    if type_factor{4}='dense',
%                    - the 4th factor will be a full/dense matrix  
%
%  
%
% Example of use :
%	 [factors,dense]=generate_Faust_config([10,25,54],{'dense','sparse'}) will generate a Faust config
%                 of 2 facts, the Faust will be of size 10x54 with the 1st fact 10x25 dense matrixof size  and the 2nd a 25x54 sparse matrix
% 		
%		at the end, you can build the coressponding Faust :
%		 F=Faust(factors);

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
function [factors, dense]=generate_Faust_config(list_dim,type_factor,density)
int_max = 100; 	
	nb_fact=length(type_factor);
	if (length(list_dim)-1) ~= nb_fact
		error('length(list_dim)-1 must be equal to length(type_factor)');
	end
	if (nb_fact < 1)
		error('the faust must have at least one factor');
	end
	
	factors=cell(1,nb_fact);
	dense=eye(list_dim(1));

	nbr_coeff_faust = 0;
	for i=1:nb_fact
		nbr_coeff_faust = nbr_coeff_faust + prod(list_dim(i:i+1));
	end


	nbr_coeff_dense_matrix = list_dim(1)*list_dim(end);

	density_per_fact = density * nbr_coeff_dense_matrix/nbr_coeff_faust;

	

	for i=1:nb_fact
		nb_row = list_dim(i);
		nb_col = list_dim(i+1);
		nnz_fact = round(density_per_fact*nb_row*nb_col);
		if (nnz_fact == 0)
			error('factor with 0 non zeros coefficent ');
		else 
		    if(nnz_fact > nb_row*nb_col)
			error('too much element per factor');
		     end
		end

		type_current_fact =  type_factor{i};
		celltype_current_fact = strsplit(type_current_fact,'.');
		
		scalarType = celltype_current_fact{1};
		sparsityType = celltype_current_fact{2};
		
			
		fact_values = double(randi(int_max,[1,nnz_fact]));

		switch scalarType
			case 'real'
								
			case 'complex'
				fact_values = fact_values + 1i*randi(int_max,[1,nnz_fact]);
			otherwise
				error('invalid scalar type');
		end

		id_1D = randperm(nb_row*nb_col,nnz_fact);
		fact_id_row = mod(id_1D-1,nb_row)+1;
		fact_id_col = floor((id_1D-1)/nb_row)+1;

		
		current_fact = sparse(fact_id_row,fact_id_col,fact_values,nb_row,nb_col,nnz_fact);

		switch sparsityType
			case 'dense'
				current_fact=full(current_fact);
			case 'sparse'
				current_fact=sparse(current_fact);
		end

		dense = dense * full(current_fact);

		
		factors{i}=current_fact;
			
	end


end
