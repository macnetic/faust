%% Description test_faust_transform
% This script tests the malab_faust class, i.e the different method
% overloaded for a faust (constructor, size, mtimes, mtimes_trans ...)

nb_fact = 3;
dim1 = 8000;
dim2 = 200;
dim3 = 10;
int_max= 100;
nb_transposition = 100;
nb_multiplication_vector = 100;

disp('****** TEST MATLAB_FAUST ******* '); 
%% creation du faust
disp(' TEST CONSTRUCTOR : '); 
factors=cell(1,nb_fact);
factors{1}=double(randi(int_max,dim1,dim2));
for i=2:nb_fact
    factors{i}=double(randi(int_max,dim2,dim2)); 
end

F = matlab_faust(factors);

disp('Ok');





%






%%% transpose test

disp('TEST TRANSPOSE : ');
disp('operation F_trans = F'' ');
t_trans=zeros(nb_transposition,1);

for i=1:nb_transposition
	disp([int2str(i) '/' int2str(nb_transposition)]);
	tic		
	F_trans=F';
	t_trans(i)=toc;	
end


disp(['time trans ' num2str(mean(t_trans)) ]); 




disp('Ok');






















%% test faust multiplication with vector
disp('TEST MULTIPLICATION BY A VECTOR : ');
x=zeros(dim2,1);
x(:)=1:dim2;
x_trans=zeros(dim1,1);
x_trans(:)=1:dim1;

F_dense=get_product(F);
y_expected = F_dense*x;
y_expected_trans = F_dense'*x_trans;

t_times=zeros(nb_multiplication_vector,1);
t_trans_times=zeros(nb_multiplication_vector,1);
t_mult_fun=zeros(nb_multiplication_vector,1);
t_trans_mult_fun=zeros(nb_multiplication_vector,1);





for i=1:nb_multiplication_vector
	disp([int2str(i) '/' int2str(nb_multiplication_vector)]);
		
	%% F*x
	tic
	y_star = F*x;
	t_times(i)=toc;
	
	if (y_expected ~= y_star)
    		error(['multiplication faust-vector : invalid result within the precision ' num2str(threshold)]);
	end

	%% F'*x
	tic
	y_star_trans = F'*x_trans;
	t_trans_times(i) = toc;	
	
	if (y_expected_trans~= y_star_trans)
    		error(['multiplication faust-vector with transposition : invalid result within the precision ' num2str(threshold)]);
	end
	
	%% mtimes_trans(F,x_trans,'T');
	tic
	y_mtimes_trans = mtimes_trans(F,x_trans,'T');
	t_mtimes_trans(i) = toc;	
	
	if (y_expected_trans ~= y_mtimes_trans)
	    error(['multiplication faust-vector with transposition : invalid result within the precision '  num2str(threshold)]);
	end

	
	
	%% mtimes_trans(F,x,'N');
	tic
	y_mtimes = mtimes_trans(F,x,'N');
	t_mtimes(i) = toc;	
	if (y_expected ~= y_mtimes)
	    error(['multiplication faust-vector : invalid result within the precision '  num2str(threshold)]);
	end



end
disp(['tps A=A'' ' num2str(mean(t_trans))]);
disp(['tps  A*x  :  ' num2str(mean(t_times))]);
disp(['tps  A''*x  :  ' num2str(mean(t_trans_times))]);
disp(['tps  mtimes_trans(F,x,''N'')  :  ' num2str(mean(t_mtimes))]);
disp(['tps  mtimes_trans(F,x_trans,''T'')  :  ' num2str(mean(t_mtimes_trans))]);



disp('Ok');



