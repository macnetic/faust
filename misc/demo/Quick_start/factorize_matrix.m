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
%%

% number of row of the matrix
dim1 = 100;
% number of column of the matrix
dim2 = 200;
% matrix to factorise
A = rand(dim1,dim2);



% Rational complexity Gain (theoretical speed-up) of the Faust
rcg = 100;
% number of factor of the Faust
nb_factor = 2;
%% generate parameters of the factorisation
params = generate_params(dim1,dim2,nb_factor,rcg);

%% factorisation (create Faust from matrix A)
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
