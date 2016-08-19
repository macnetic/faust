%% Description test_matlab_faust2
% This script test the Faust class methods with incorrect input argument
% they must throw an exception
nb_fact = 3;
dim1 = 5;
dim2 = 4;
dim3 = 10;
int_max= 100;
threshold = 10^(-5);

disp('****** TEST MATLAB_FAUST ******* ');
disp('CONSTRUCTOR ');
disp('test 1 : '); 
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


disp('test 2 : '); 
test_pass = 0;
expected_err_message='addSpmat : empty matlab matrix';
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




disp('test 3 : ');
test_pass = 0;
expected_err_message='getFaustMat :input matrix format must be single or double';
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







