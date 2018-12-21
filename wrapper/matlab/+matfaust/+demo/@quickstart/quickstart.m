%> @package matfaust.demo @brief The matfaust namespace for the demos partly based on research papers

%% Description quick_start
%
%  This demo shows that a Faust is handled as matlab builtin matrix,
% presenting functions that are overloaded by the Faust class
% (size,mtimes,transpose...)
% and ends with a little time comparison to illustrate
% the speed-up of using a Faust for multiplication.
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
%% ---------------------------------------------------------------------------------------
%% Description factorize_matrix
%
%  This demo presents the method to factorize a given matrix into a Faust
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
%%------------------------------------------------------------------------------
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

% =====================================================================
%> The FAµST quickstart script, a good place to look at for a first tour.
% =====================================================================
classdef quickstart
	methods(Static)
		%==========================================================================
		%> This demo shows that a Faust is handled as matlab builtin matrix, presenting functions that are overloaded by the Faust class (size, mtimes, transpose...) and ends with a little time comparison to illustrate the speed-up of using a Faust for multiplication.
		%==========================================================================
		function quick_start()
			import matfaust.Faust


			% loading a Faust A from saved-one
			A=Faust('faust_quick_start.mat')

			% get the size of the faust
			[dim1,dim2] = size(A)

			%  transpose a faust
			A_trans = A'

			% multiplication by A
			x1 = rand(dim2,1);
			y1 = A*x1;

			% multiplication by A'
			x2 = rand(dim1,5);
			y2 = A'*x2;


			% get the 2-norm (spectral norm) of the faust A
			norm_A = norm(A) % equivalent to norm(A,2);

			% convert Faust to full matrix
			A_full=full(A);

			% get the number of non-zeros coefficient
			nz = nnz_sum(A)

			% READING coefficient
			coeff=A(3,4)
			col_2=A(:,2);
			submatrix_A=A(3:5,2:3)
			submatrix_A=A(end-5:end,4980:end-1)

			% WARNING : WRITING coefficient NOT ALLOWED
			% A(i,j)=3; will throw an error with this message :
			%            'function not implemented for Faust class'





			%% speed-up multiplication
			nb_mult=100;
			time_full=0;
			time_faust=0;

			for i=1:nb_mult
				tic
				y=A_full*x1;
				time_full=time_full+toc;

				tic
				y=A*x1;
				time_faust=time_faust+toc;
			end

			disp('multiplication SPEED-UP using Faust');
			disp(['Faust is ' num2str(time_full/time_faust) ' faster than a full matrix']);
		end

		%============================================================
		%> This demo presents the method to factorize a given matrix into a Faust.
		%============================================================
		function factorize_matrix()
			% number of row of the matrix
			dim1 = 100;
			% number of column of the matrix
			dim2 = 200;
			% matrix to factorize
			A = rand(dim1,dim2);



			% Rational complexity Gain (theoretical speed-up) of the Faust
			rcg_ = 100;
			% number of factor of the Faust
			nb_factor = 2;
			%% generate parameters of the factorization
			params = generate_params(dim1,dim2,nb_factor,rcg_);

			%% factorization (create Faust from matrix A)
			faust_A = faust_decompose(A,params);



			%% speed-up multiplication
			y=zeros(dim2,1);
			nb_mult=100;
			time_full=0;
			time_faust=0;
			x1=rand(dim2,1);

			for i=1:nb_mult
				tic
				y=A*x1;
				t1=toc;

				tic
				y=faust_A*x1;
				t2=toc;


				time_full=time_full+t1;
				time_faust=time_faust+t2;

			end

			disp('multiplication SPEED-UP using Faust');
			disp(['Faust is ' num2str(time_full/time_faust) ' faster than a full matrix']);
		end
		%=======================================================================
		%> This demo presents the method to build a Faust from its factors
		%=======================================================================
		function construct_Faust_from_factors()
			import matfaust.Faust
			% number of rows of the Faust
			dim1=300;
			% number of columns of the Faust
			dim2=100;
			% number of factors of the Faust
			nb_factor=4;
			% density of each factor
			density_factor=0.1;

			%cell-array representing its factors
			factors = cell(1,nb_factor);

			% 1st factor is a rectangular sparse matrix of density equal to density_factor
			factors{1} = sprand(dim1,dim2,density_factor);

			% other factors are square sparse matrices of density equal to density_factor
			for i=2:nb_factor
				factors{i} = sprand(dim2,dim2,density_factor);
			end

			%% construct the Faust
			A_faust=Faust(factors);


			% a multiplicative scalar can be taken into account to construct the Faust
			lambda=2;

			% B_faust=lambda*A_faust
			B_faust=Faust(factors,lambda);
		end
	end
end
