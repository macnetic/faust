%% Description test_faust_transform
% This script tests the malab_faust class from a time computing point of view

nb_fact = 3;
dim1 = 1000;
dim2 = 500;
dim3 = 10;
int_max= 100;



nb_transposition = 100;
nb_multiplication_vector = 100;
nb_access_coeff = 100;
nb_norm = 100;
threshold = 10^(-5);




disp('****** TEST MATLAB_FAUST ******* '); 
%% creation du faust
disp(' TEST CONSTRUCTOR : '); 
factors=cell(1,nb_fact);
factors{1}=double(randi(int_max,dim1,dim2));
for i=2:nb_fact
    factors{i}=double(randi(int_max,dim2,dim2)); 
end

F = Faust(factors);

disp('Ok');











%%% transpose test

disp('TEST TRANSPOSE : ');
disp('operation F_trans = F'' ');
t_trans=zeros(nb_transposition,1);

for i=1:nb_transposition
	%disp([int2str(i) '/' int2str(nb_transposition)]);
	tic		
	F_trans=F';
	t_trans(i)=toc;	
end


disp(['time trans  : ' num2str(mean(t_trans)) ]); 




disp('Ok');










%% access to coeff
disp('TEST ACCESS ROW OR COL');
F_dense=full(F);
t_access_row=zeros(nb_transposition,1);
t_access_col=zeros(nb_transposition,1);
t_access_trans_row=zeros(nb_transposition,1);
t_access_trans_col=zeros(nb_transposition,1);
for i=1:nb_access_coeff
	%disp([int2str(i) '/' int2str(nb_access_coeff)]);
	tic
	col=F(:,1);
	t_access_col(i)=toc;

	tic	
	row=F(1,:);
	t_access_row(i)=toc;

	tic
	col_trans=F_trans(:,1);
	t_access_trans_col(i)=toc;

	tic
	row_trans=F_trans(1,:);
	t_access_trans_row(i)=toc;

	

	if (col ~=F_dense(:,1))	
		error('access to the col');
	end

	if (row ~=F_dense(1,:))
		error('access to the row');
	end

	if (col_trans ~=F_dense(1,:)')	
		error('access to the col');
	end	

	if (row_trans ~=F_dense(:,1)')	
		error('access to the row');
	end
end

disp(['tps F(1,:) : ' num2str(mean(t_access_row))]);
disp(['tps F(:,1) : ' num2str(mean(t_access_col))]);
disp(['tps F_trans(1,:) : ' num2str(mean(t_access_trans_row))]);
disp(['tps F_trans(:,1) : ' num2str(mean(t_access_trans_col))]);












%% test faust multiplication with vector
disp('TEST MULTIPLICATION BY A VECTOR : ');
istransposed=1;
nontransposed=0;
x=zeros(dim2,1);
x(:)=1:dim2;
x_trans=zeros(dim1,1);
x_trans(:)=1:dim1;


y_expected = F_dense*x;
y_expected_trans = F_dense'*x_trans;

t_times=zeros(nb_multiplication_vector,1);
t_trans_times=zeros(nb_multiplication_vector,1);
t_mult_fun=zeros(nb_multiplication_vector,1);
t_trans_mult_fun=zeros(nb_multiplication_vector,1);





for i=1:nb_multiplication_vector
	%disp([int2str(i) '/' int2str(nb_multiplication_vector)]);
		
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
	
	%% mtimes_trans(F,x_trans,istransposed);
	tic
	y_mtimes_trans = mtimes_trans(F,x_trans,istransposed);
	t_mtimes_trans(i) = toc;	
	
	if (y_expected_trans ~= y_mtimes_trans)
	    error(['multiplication faust-vector with transposition : invalid result within the precision '  num2str(threshold)]);
	end

	
	
	%% mtimes_trans(F,x,nontransposed);
	tic
	y_mtimes = mtimes_trans(F,x,nontransposed);
	t_mtimes(i) = toc;	
	if (y_expected ~= y_mtimes)
	    error(['multiplication faust-vector : invalid result within the precision '  num2str(threshold)]);
	end



end

disp(['tps A=A'' : ' num2str(mean(t_trans))]);
disp(['tps  A*x  :  ' num2str(mean(t_times))]);
disp(['tps  A''*x  :  ' num2str(mean(t_trans_times))]);
disp(['tps  mtimes_trans(F,x,nontransposed)  :  ' num2str(mean(t_mtimes))]);
disp(['tps  mtimes_trans(F,x_trans,istransposed)  :  ' num2str(mean(t_mtimes_trans))]);


disp('Ok');


disp('TEST 2-norm : ');

t_dense_norm=zeros(nb_norm,1);
t_faust_norm=zeros(nb_norm,1);

for i=1:nb_norm
	tic
	norm_F=norm(F);
	t_faust_norm(i) = toc;
	
	tic	
	norm_F_dense=norm(F);
	t_dense_norm = toc;

	if (abs(norm_F-norm_F_dense) > threshold)
		error(['norm : invalid result, expects ' num2str(norm_F_dense) ' but get ' num2str(norm_F)]);
	end	
end

disp(['norm dense  :  ' num2str(mean(t_dense_norm))]);
disp(['norme F  :  ' num2str(mean(t_faust_norm))]);






