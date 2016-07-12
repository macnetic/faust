%% Description test_faust_transform
% This script tests the malab_faust class, i.e the different method
% overloaded for a faust (constructor, size, mtimes, mtimes_trans ...)

nb_fact = 3;
dim1 = 5;
dim2 = 3;
dim3 = 10;
int_max= 100;
threshold = 10^(-15); % equality test

%% creation du faust
factors=cell(1,nb_fact);
factors{1}=double(randi(int_max,dim1,dim2));
for i=2:nb_fact
    factors{i}=double(randi(int_max,dim2,dim2)); 
end

F = matlab_faust(factors);





%% size test
[dim1_test,dim2_test]=size(F);
if ((dim1_test ~= dim1) | (dim2_test ~= dim2))
    error(['size : invalid output dimension']);
end
    
dim1_test1=size(F,1);
dim2_test1=size(F,2);

if ((dim1_test1 ~= dim1) | (dim2_test1 ~= dim2))
    error(['size : invalid output dimension']);
end


%% get_nb_factor test
nb_fact_test=get_nb_factor(F);
if (nb_fact_test ~= nb_fact)
	error('get_nb_factor : invalid number of factor of the faust');
end

%% get_fact test
for i=1:nb_fact
	A=get_fact(F,i);
	if(A~=factors{i})
		error('get_fact : invalid factor');
	end

end


%% load_faust and save_faust test
filename = 'faust.mat';
save_faust(F,filename);
F_loaded = matlab_faust(filename);
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



%% get_product test
F_dense= get_product(F);

[dim1_dense,dim2_dense]=size(F_dense);
if((dim1_dense ~= dim1) | (dim2_dense ~= dim2))
    error('get_product : invalid dimension');
end

%% transpose test
F_trans = F';
[dim1_trans,dim2_trans]=size(F_trans);
F_dense_trans=get_product(F_trans);

if ((dim1_trans ~= dim2) | (dim2_trans ~= dim1))
    error(['transpose : invalid dimension']);
end

F_dense_trans=get_product(F_trans);

if (F_dense' ~= F_dense_trans)
   error(['transpose : invalid value']); 
end






%% test faust multiplication with vector

x=zeros(dim2,1);
x(:)=1:dim2;
x_trans=zeros(dim1,1);
x_trans(:)=1:dim1;

y_expected = F_dense*x;
y_expected_trans = F_dense'*x_trans;

y_star = F*x;
if (y_expected ~= y_star)
    error(['multiplication faust-vector : invalid result within the precision ' num2str(threshold)]);
end



y_star_trans = F'*x_trans;
if (y_expected_trans~= y_star_trans)
    error(['multiplication faust-vector with transposition : invalid result within the precision ' num2str(threshold)]);
end


y_mtimes_trans = mtimes_trans(F,x_trans,'T');
if (y_expected_trans ~= y_mtimes_trans)
    error(['multiplication faust-vector with transposition : invalid result within the precision '  num2str(threshold)]);
end


y_mtimes = mtimes_trans(F,x,'N');
if (y_expected ~= y_mtimes)
    error(['multiplication faust-vector : invalid result within the precision '  num2str(threshold)]);
end



%% test multiplication with matrix
X=zeros(dim2,dim3);
X(:)=1:dim2*dim3;
X_trans=zeros(dim1,dim3);
X_trans(:)=1:dim1*dim3;

Y_expected = F_dense*X;
Y_expected_trans = F_dense'*X_trans;



Y_star = F*X;
if (Y_expected ~= Y_star)
    error(['multiplication faust-vector : invalid result within the precision ' num2str(threshold)]);
end



Y_star_trans = F'*X_trans;
if (Y_expected_trans~= Y_star_trans)
    error(['multiplication faust-vector with transposition : invalid result within the precision ' num2str(threshold)]);
end


Y_mtimes_trans = mtimes_trans(F,X_trans,'T');
if (Y_expected_trans ~= Y_mtimes_trans)
    error(['multiplication faust-vector with transposition : invalid result within the precision '  num2str(threshold)]);
end


Y_mtimes = mtimes_trans(F,X,'N');
if (Y_expected ~= Y_mtimes)
    error(['multiplication faust-vector : invalid result within the precision '  num2str(threshold)]);
end

