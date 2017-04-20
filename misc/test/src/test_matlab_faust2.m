%% Description test_matlab_faust2
% This script test the Faust class methods with incorrect input argument
% they must throw an exception
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
nb_fact = 3;
dim1 = 5;
dim2 = 4;
dim3 = 10;
int_max= 100;
threshold = 10^(-5);

disp('****** TEST MATLAB_FAUST ******* ');
disp('CONSTRUCTOR ');
disp('test 1 : invalid factor size '); 
test_pass = 0;
expected_err_message='Faust::Transform<FPP,Cpu> : check_factors_validity : dimensions of the factors mismatch';
factors=cell(1,nb_fact);
factors{1}=double(randi(int_max,dim1,dim2+1));%% incorrect 1st factor size(factors{1},2) ~= size(factors{2},1)
for i=2:nb_fact
    factors{i}=double(randi(int_max,dim2,dim2)); 
end


try
	F = Faust(factors); % this must throw an exception
catch ME
	if strcmp(ME.message,expected_err_message)
		test_pass = 1;
	else
		error([ 'error with a wrong message : ' ME.message ' must be : ' expected_err_message]);  	
	end
end

if(~test_pass)
	error('failure');
end	 

disp('Ok');



disp('test 2 : invalid factor empty'); 
test_pass = 0;
expected_err_message='concatMatGeneric : empty matlab matrix';
factors=cell(1,nb_fact); % each cell is empty, must contained a matrix



try
	F = Faust(factors); % this must throw an exception
catch ME
	if strcmp(ME.message,expected_err_message)
		test_pass = 1;
	else
		error([ 'error with a wrong message : ' ME.message ' must be : ' expected_err_message ]);  	
	end
end

if(~test_pass)
	error('failure');
end	 

disp('Ok');




disp('test 3 : invalid factor type (cell)');
test_pass = 0;
expected_err_message='Unknown matlab type.';
factors=cell(1,1);
factors{1}=cell(1); % contains a character and not a matrix 



try
	F = Faust(factors); % this must throw an exception
catch ME
	if strcmp(ME.message,expected_err_message)
		test_pass = 1;
	else
		error([ 'error with a wrong message : ' ME.message ' must be : ' expected_err_message ]);  	
	end
end

if(~test_pass)
	error('failure');
end	 

disp('Ok');

disp('test 4 : subsasgn not implemented');
test_pass = 0;
expected_err_message='function not implemented for Faust class';
F=Faust({ones(5,4),ones(4,7)});

try
	F(1,2)=3
catch ME
	if strcmp(ME.message,expected_err_message)
		test_pass = 1;
	else
		error([ 'error with a wrong message : ' ME.message ' must be : ' expected_err_message ]);  	
	end
end

if(~test_pass)
	error('failure');
end	 
disp('Ok');




disp('test 5 : ctranspose is not yet implemented for complex scalar Faust');
test_pass = 0;
expected_err_message='ctranspose is not yet implemented for complex scalar Faust';
F=Faust({ones(5,4),ones(4,7)+1i*ones(4,7)});

try
	F_ctrans = F';
catch ME
	if strcmp(ME.message,expected_err_message)
		test_pass = 1;
	else
		error([ 'error with a wrong message : ' ME.message ' must be : ' expected_err_message ]);  	
	end
end

if(~test_pass)
	error('failure');
end	 
disp('Ok');




 






