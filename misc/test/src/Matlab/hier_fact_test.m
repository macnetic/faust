%% Description
%  HIER_FACT_TEST.
%  tests the hierarchical_fact algorithm,
%  checks the factorisation and makes time measurement :
%   
%  INPUT :
% - paramsFile : filename where the matrix to be factorized is stored 
% 
% - paramsfile : filename where the params of hierarchical_fact are stored
%
% - expectedLambda : the expected lambda value of the factorisation
%
% - expectedLambdaPrecision : to succeed |expectedLambda - lambda|< expectedLambdaPrecision must be true
%
% - opt : if opt = 'MEX' => mexHierarchical_fact is called (mexfunction)
%            opt = 'MATLAB' => old_hierarchical_fact is called (pure matlab version) 
%
%  
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
function hier_fact_test(matrixFile,paramsFile,expectedLambda, expectedLambdaPrecision,opt)



%% load the hierarchical_fact configuration
disp(['*** LOADING MATRIX FILE ***']);
disp([matrixFile]);
disp(' ');
disp(' ');
 
load(matrixFile);

disp(['*** LOADING PARAMS FILE ***']);
disp([paramsFile]);
disp(' ');
disp(' ');
 
load(paramsFile);




%% factorization of the matrix 
switch(opt)
	case 'MEX' %% factorisation (mexfile)
		disp('*** MEX FACTORISATION ***');  
		tic
			[lambda,fact]=mexHierarchical_fact(matrix,params);
		t=toc;
	case 'MATLAB' %% factorisation (matlab)
		disp(['*** MATLAB FACTORISATION ***']); 
		tic
			[lambda,fact]=old_hierarchical_fact(matrix,params);
		t=toc;
     
	otherwise
		error('hier_fact_test : invalid opt parameter , must be equal to ''MEX'' or ''MATLAB''');
end


% conversion seconds into [hour,min,seconds])
t=mod(t, [0, 3600, 60]) ./ [3600, 60, 1];
t(1:2)=floor(t(1:2));


disp(['time factorization : ' int2str(t(1)) ' hour ' int2str(t(2)) ' min ' num2str(t(3)) ' seconds']);  



%% check if the result are ok
disp(['lambda value : ' num2str(lambda)]);
if (abs(lambda - expectedLambda) > expectedLambdaPrecision)
    disp(' ');
    
    disp([ 'expected lamba value : ' int2str(expectedLambda) ' in the precision of ' int2str(expectedLambdaPrecision) ]);	
    error('invalid lambda value');
end


F=Faust(fact,lambda);
relative_error = norm(matrix - full(F));

disp(['relative error :  ' num2str(relative_error)]);



%% speed-up test for multiplication with a vector
disp(' ');
disp(' ');
disp('*** product data matrix-vector vs product faust-vector ***');
nbiter = 100;
[nl,nc]=size(matrix);
y_faust=zeros(nl,1);
y_dense=zeros(nl,1);
tps_dense=zeros(nbiter,1);
tps_faust=zeros(nbiter,1);

for i=1:nbiter
    x=rand(nc,1);
   
    tic
        y_faust = F*x;
    tps_faust=toc;
    
    tic
        y_dense = matrix*x;
    tps_dense=toc;    
end


tps_dense=mean(tps_dense);
tps_faust=mean(tps_faust);
speed_up = tps_dense/tps_faust;
disp(['dense time : ' num2str(tps_dense) ' faust time : ' num2str(tps_faust)]);
disp(['speed_up : ' num2str(speed_up)]);
