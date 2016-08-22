%% Description test_matlab_faust
% This script tests the Faust class, i.e the different method
% overloaded for a faust (constructor, size, mtimes, mtimes_trans ...)

nb_fact = 3;
dim1 = 5;
dim2 = 4;
dim3 = 10;
int_max= 100;
threshold = 10^(-5);

disp('****** TEST MATLAB_FAUST ******* '); 
%% creation du faust
disp(' TEST CONSTRUCTOR : ');
disp(['generate a faust of size ' int2str(dim1)  ' by ' int2str(dim2)]); 
factors=cell(1,nb_fact);
factors{1}=double(randi(int_max,dim1,dim2));
for i=2:nb_fact
    factors{i}=double(randi(int_max,dim2,dim2)); 
end

F = Faust(factors);
empty_F = Faust({});
disp('Ok');




%% size test
disp('TEST SIZE : ');

[dim1_test,dim2_test]=size(F);
if ((dim1_test ~= dim1) | (dim2_test ~= dim2))
    error(['size : invalid output dimension']);
end
    
dim1_test1=size(F,1);
dim2_test1=size(F,2);

if ((dim1_test1 ~= dim1) | (dim2_test1 ~= dim2))
    error(['size : invalid output dimension']);
end




[dim1_empty,dim2_empty]=size(empty_F);

if (dim1_empty ~= 0) | (dim2_empty ~= 0)
	dim1_empty
	dim2_empty
	error('size : invalid dimension for an empty Faust');
end


disp('Ok');












disp('TEST GET_NB_FACTOR : ');
%% get_nb_factor test
nb_fact_test=get_nb_factor(F);
if (nb_fact_test ~= nb_fact)
	error('get_nb_factor : invalid number of factor of the faust');
end



disp('Ok');


%% get_fact test
disp('TEST GET_FACT : ');
for i=1:nb_fact
	A=get_fact(F,i);
	if(A~=factors{i})
		error('get_fact : invalid factor');
	end

end
disp('Ok');




%% get_product test
disp('TEST GET_PRODUCT : ');
F_dense= get_product(F);

[dim1_dense,dim2_dense]=size(F_dense);

if((dim1_dense ~= dim1) | (dim2_dense ~= dim2))
    error('get_product : invalid dimension');
end
disp('Ok');



%% get_product test
disp('TEST NNZ : ');
expected_nz = (nb_fact-1)*dim2^2 + dim1*dim2;

nz = nnz(F);

if(expected_nz ~= nz)
    error('nnz : invalid number of nonzeros');
end

new_facts{1}=eye(dim1,dim1);
new_facts{2}=eye(dim1,dim1);
new_facts{1}(1,1)=0;
expected_nz_2  = 2*dim1-1;

F_nnz = Faust(new_facts);

nz_F = nnz(F_nnz);

if(nz_F ~= expected_nz_2)	
    error('nnz : invalid number of nonzeros');
end


expected_nz_empty_F = 0;
nz_empty_F = nnz(empty_F);

if(nz_empty_F ~= expected_nz_empty_F)	
    error('nnz : invalid number of nonzeros (empty Faust)');
end
 

disp('Ok');



disp('TEST DENSITY : ');
expected_density = expected_nz / (dim1*dim2);
dens = density(F);
if(dens ~= expected_density)
    error('density : invalid value');
end

expected_density2 = expected_nz_2 / dim1^2;
dens2 = density(F_nnz);

if(dens2 ~= expected_density2)
    error('density : invalid value');
end


expected_density_empty_F = -1;
density_empty_F = density(empty_F);
if(density_empty_F ~= expected_density_empty_F)
    error('density : invalid value');
end



disp('Ok');


disp('TEST RCG : ');
expected_RCG = 1/expected_density;
rcg_F = RCG(F);
if(rcg_F ~= expected_RCG)	
	rcg_F
	expected_RCG
    error('RCG : invalid value');
end

expected_RCG2 =  1/expected_density2;
rcg2_F = RCG(F_nnz);

if(rcg2_F ~= expected_RCG2)
    rcg2_F
    expected_RCG2	 	
    error('RCG : invalid value');
end

expected_RCG_empty_F = -1;
RCG_empty_F = RCG(empty_F);
if(RCG_empty_F ~= expected_RCG_empty_F)
    error('RCG : invalid value');
end

disp('Ok');





















%%% transpose test
disp('TEST TRANSPOSE : ');
disp('operation F_trans = F'' ');
F_trans=F';
[dim1_trans,dim2_trans]=size(F_trans);
if ((dim1_trans ~= dim2) | (dim2_trans ~= dim1))
    error(['transpose : invalid dimension']);
end

F_dense_trans = get_product(F_trans);
if (F_dense_trans ~= F_dense')
    error(['transpose : invalid transpose matrix']);
end

%% verification de la non modification du faust
[new_dim1,new_dim2]=size(F); 
if ((new_dim1 ~= dim1) | (dim2 ~= new_dim2))
    error(['transpose : modification du faust de depart']);
end

new_F_dense=get_product(F);
if((new_F_dense ~= F_dense))
	error('transpose : modification du faust de depart');
end 


disp('operation F_trans_trans = F_trans'' ');
F_trans_trans=F_trans';
[dim1_trans_trans,dim2_trans_trans]=size(F_trans_trans);

if ((dim1_trans_trans ~= dim1) | (dim2_trans_trans ~= dim2))
    error(['transpose : invalid dimension']);
end

F_dense_trans_trans = get_product(F_trans_trans);
if (F_dense_trans_trans ~= F_dense)  
    error(['transpose : invalid transpose matrix']);
end


%% verification de la non modification du faust
[new_dim1_trans,new_dim2_trans]=size(F_trans); 
if ((new_dim1_trans ~= dim1_trans) | (new_dim2_trans ~= new_dim2_trans))
    error(['transpose : modification du faust de depart']);
end

new_F_dense_trans=get_product(F_trans);
if((new_F_dense_trans ~= F_dense_trans))
	error('transpose : modification du faust de depart');
end 





disp('Ok');

%% load_faust and save_faust test
disp('TEST LOAD AND SAVE : ');
filename = [ '@FAUST_BIN_TEST_OUTPUT_DIR@' filesep 'faust.mat'];
disp(['save faust into the file : ' filename]); 
save(F,filename);
F_loaded = Faust(filename);
[dim1_loaded,dim2_loaded]=size(F_loaded);

if (dim1_loaded ~= dim1) | (dim2_loaded ~= dim2)
	error('load and save faust : invalid dimension');
end

nb_fact_load=get_nb_factor(F);
if (nb_fact_load ~= nb_fact)
	error('load and save faust : invalid number of factor of the loaded faust ');
end

for i=1:nb_fact
	A=get_fact(F_loaded,i);
	if(A~=factors{i})
		error('get_fact : invalid factor');
	end

end



filename_trans = [ '@FAUST_BIN_TEST_OUTPUT_DIR@' filesep 'faust_trans.mat'];
disp(['save transposed faust into the file : ' filename_trans]); 
save(F_trans,filename_trans);

F_trans_loaded = Faust(filename_trans);
[dim1_trans_faust_loaded,dim2_trans_faust_loaded]=size(F_trans_loaded);

if (dim1_trans_faust_loaded ~= dim2) | (dim2_trans_faust_loaded ~= dim1)
	error(['save transposed : invalid dimension to the loaded-saved faust']);
end

F_dense_trans_loaded=get_product(F_trans_loaded);
if (F_dense_trans_loaded ~= F_dense')
	error(['save transposed : invalid faust']);
end





disp('Ok');

%% slicing test
disp('TEST SLICING : ');
for i=1:dim1
	for j=1:dim2
		F_i_j=F(i,j);
		F_trans_j_i=F_trans(j,i);

		if (size(F_i_j) ~= [1 1])
			error('invalid size of F(i,j)');
		end
		if (size(F_trans_j_i) ~= [1 1])
			error('invalid size of F_trans(j,i)');
		end
		if (F_i_j ~= F_dense(i,j))
			error('F(i,j) ~= F_dense(i,j)');
		end
		if (F_trans_j_i ~= F_dense_trans(j,i))
			error('F(j,i) ~= F_dense_trans(j,i)');
		end
	end
end

F_slice_slice=F(:,:);
if (size(F_slice_slice,1) ~= dim1) | (size(F_slice_slice,2) ~= dim2)
	error('invalid dimension');
end
if (F_slice_slice ~= F_dense)
	error('F(:,:) ~= F_dense');
end


F_trans_slice_slice=F_trans(:,:);
if (size(F_trans_slice_slice,1) ~= dim2) | (size(F_trans_slice_slice,2) ~= dim1)
	error('invalid dimension');
end
if (F_trans_slice_slice ~= F_dense')
	error('F_trans(:,:) ~= F_dense''');
end


F_slice_slice_2=F(1:dim1,1:dim2);
if (size(F_slice_slice_2,1) ~= dim1) | (size(F_slice_slice_2,2) ~= dim2)
	error('invalid dimension');
end
if (F_slice_slice_2 ~= F_dense)
	error('F(1:dim1,1:dim2) ~= F_dense');
end

F_inv=F(dim1:-1:1,dim2:-1:1);
if (size(F_inv,1) ~= dim1) | (size(F_inv,2) ~= dim2)
	error('invalid dimension');
end


if (F_inv ~= F_dense(dim1:-1:1,dim2:-1:1))
	error('F(1:dim1,1:dim2) ~= F_dense(dim1:-1:1,dim2:-1:1)');
end 


%% operator end with slicing
F_end = F(1:end,1:end);
if (size(F_end,1) ~= dim1) | (size(F_end,2) ~= dim2)
	error('invalid dimension');
end 

if (F_end ~= F_dense(1:end,1:end))
	error('F(1:end,1:end) ~= F_dense(1:end,1:end)');
end



disp('Ok');































%% test operator=
disp('TEST OPERATOR= : ');

F_eq=F;
F_trans_eq=F_trans;
F_trans_trans_eq=F_trans_trans;

[dim1_eq, dim2_eq]=size(F_eq);
[dim1_trans_eq, dim2_trans_eq]=size(F_trans_eq);
[dim1_trans_trans_eq, dim2_trans_trans_eq]=size(F_trans_trans_eq);
if ((dim1_eq ~= dim1) | (dim2_eq ~= dim2))
    error(['operator = test 1 : invalid dimension']);
end
if ((dim1_trans_eq ~= dim2) | (dim2_trans_eq ~= dim1))
    error(['operator = test 2 : invalid dimension']);
end
if ((dim1_trans_trans_eq ~= dim1) | (dim2_trans_trans_eq ~= dim2))
    error(['operator = test 3 : invalid dimension']);
end
disp('Ok');



%% test faust multiplication with vector
disp('TEST MULTIPLICATION BY A VECTOR : ');
x=zeros(dim2,1);
x(:)=1:dim2;
x_trans=zeros(dim1,1);
x_trans(:)=1:dim1;
istransposed=1;
nontransposed=0;
y_expected = F_dense*x;
y_expected_trans = F_dense'*x_trans;

y_star = F*x;
if (y_expected ~= y_star)
    error(['multiplication faust-vector : invalid result  ' ]);
end



y_star_trans = F'*x_trans;
if (y_expected_trans~= y_star_trans)
    error(['multiplication faust-vector with transposition : invalid result  ' ]);
end


y_mtimes_trans = mtimes_trans(F,x_trans,istransposed);
if (y_expected_trans ~= y_mtimes_trans)
    error(['multiplication faust-vector with transposition : invalid result  '  ]);
end


y_mtimes = mtimes_trans(F,x,nontransposed);
if (y_expected ~= y_mtimes)
    error(['multiplication faust-vector : invalid result  '  ]);
end


y_mtimes_trans_N = mtimes_trans(F_trans,x_trans,nontransposed);
if (y_expected_trans ~= y_mtimes_trans_N)
    error(['multiplication faust-vector with transposition : invalid result  '  ]);
end


y_mtimes_trans_T = mtimes_trans(F_trans,x,istransposed);
if (y_expected ~= y_mtimes_trans_T)
    error(['multiplication faust-vector : invalid result  '  ]);
end


disp('Ok');













%% test multiplication with matrix
disp('TEST MULTIPLICATION BY A MATRIX : ');
X=zeros(dim2,dim3);
X(:)=1:dim2*dim3;
X_trans=zeros(dim1,dim3);
X_trans(:)=1:dim1*dim3;

Y_expected = F_dense*X;
Y_expected_trans = F_dense'*X_trans;



Y_star = F*X;
if (Y_expected ~= Y_star)
    error(['multiplication faust-vector : invalid result  ' ]);
end



Y_star_trans = F'*X_trans;
if (Y_expected_trans~= Y_star_trans)
    error(['multiplication faust-vector with transposition : invalid result  ' ]);
end


Y_mtimes_trans = mtimes_trans(F,X_trans,istransposed);
if (Y_expected_trans ~= Y_mtimes_trans)
    error(['multiplication faust-vector with transposition : invalid result  '  ]);
end


Y_mtimes = mtimes_trans(F,X,nontransposed);
if (Y_expected ~= Y_mtimes)
    error(['multiplication faust-vector : invalid result  '  ]);
end

Y_mtimes_trans_N = mtimes_trans(F_trans,X_trans,nontransposed);
if (Y_expected_trans ~= Y_mtimes_trans_N)
    error(['multiplication faust-vector with transposition : invalid result  '  ]);
end


Y_mtimes_trans_T = mtimes_trans(F_trans,X,istransposed);
if (y_expected ~= y_mtimes_trans_T)
    error(['multiplication faust-vector : invalid result  '  ]);
end

disp('Ok');





%% test 2-norm (spectral norm)
disp('TEST 2-norm : ');
real_norm=norm(F_dense);
norm_faust=norm(F);
norm_faust2=norm(F,2);
norm_faust_trans=norm(F_trans);

if (abs(real_norm - norm_faust)>threshold)
	error(['norm : invalid result, expected ' num2str(real_norm) ' get norm_faust' num2str(norm_faust)]);
end
if (norm_faust ~= norm_faust2)
	error(['norm : norm(F) must be equal to norm(F,2)']);
end

if (norm_faust_trans ~= norm_faust)
	error(['norm : norm(F) must be equal to norm(F_trans)']);
end

disp('Ok');





