%% Description test_matlab_faust
% This functions test different method
% overloaded for a faust (constructor, size, mtimes, mtimes_trans ...)
% with a given Faust
%
% input parameter : dim1 number of row of the Faust F
%		    dim2 number of column of the Faust F
%		    dim3 number of column of the matrix that will be multiplied by the Faust
%                   nb_fact number of factor of the Faust
%
%
%  
%
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


%nb_fact = 3;
%dim1 = 5;
%dim2 = 4;
%dim3 = 10;



 function test_matlab_faust(factors,expected_F_dense,dim3,copyOptimized)
%function test_matlab_faust(dim1,dim2,dim3,nb_fact)
int_max= 100;
threshold = 0.2;


nb_fact=length(factors);
if (nb_fact == 0)
	error('empty faust is not taking into account');
end

dim1=size(factors{1},1);
dim2=size(factors{nb_fact},2);

if (size(expected_F_dense,1) ~= dim1) | (size(expected_F_dense,2) ~= dim2)
	error('the factor dimension mismatch expected_F_dense');
end


isComplex=~isreal(expected_F_dense);

if (isComplex)
	scalarType='complex';
else
	scalarType='real';
end


expected_nz = 0;
for i=1:nb_fact
	expected_nz = expected_nz + nnz(factors{i});
end

expected_density = expected_nz/(dim1*dim2);

disp('****** TEST MATLAB_FAUST ******* '); 
disp([' CONFIG OF THE ' scalarType ' scalar FAUST ']);
disp(['number of row of the Faust : ' int2str(dim1)]); 
disp(['number of column of the Faust : ' int2str(dim2)]);
disp(['number of factor of the Faust : ' int2str(nb_fact)]);
disp(['density of the Faust : ' num2str(expected_density)]);
disp('');
disp(['number of column of the matrix that will be multiplied by the Faust : ' int2str(dim3)]);
disp('');

%% creation du faust
disp(' TEST CONSTRUCTOR : ');
F = Faust(factors,1.0,copyOptimized);
empty_F=Faust({});
disp('Ok');

%% size test
disp('TEST DISP : ');
disp(F);
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







%% full test
disp('TEST FULL : ');
F_dense= full(F);
[dim1_dense,dim2_dense]=size(F_dense);

if((dim1_dense ~= dim1) | (dim2_dense ~= dim2))
    error('full : invalid dimension');
end
if ( ~isequal(expected_F_dense,F_dense) ) 
    expected_F_dense
    F_dense	 
    error(['full : invalid full-storage matrix']);
end
disp('Ok');



%% isreal test
disp('TEST ISREAL : ');
expected_bool_isreal=isreal(expected_F_dense);
bool_isreal=isreal(F);

if (expected_bool_isreal ~= bool_isreal)
	error('isreal : invalid output boolean'); 
end
disp('Ok');


%% full test
disp('TEST NNZ : ');


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
rcg_F = rcg(F);
if(rcg_F ~= expected_RCG)	
	rcg_F
	expected_RCG
    error('RCG : invalid value');
end

expected_RCG2 =  1/expected_density2;
rcg2_F = rcg(F_nnz);

if(rcg2_F ~= expected_RCG2)
    rcg2_F
    expected_RCG2	 	
    error('RCG : invalid value');
end

expected_RCG_empty_F = -1;
RCG_empty_F = rcg(empty_F);
if(RCG_empty_F ~= expected_RCG_empty_F)
    error('RCG : invalid value');
end

disp('Ok');





















%%% transpose test
disp('TEST TRANSPOSE : ');
disp('operation F_trans = F'' ');
F_trans=F.';
[dim1_trans,dim2_trans]=size(F_trans);
if ((dim1_trans ~= dim2) | (dim2_trans ~= dim1))
    error(['transpose : invalid dimension']);
end

F_dense_trans = full(F_trans);

%(F_dense_trans ~= F_dense')
if (~isequal(F_dense_trans,F_dense.'))    
	error(['transpose : invalid transpose matrix']);
end

%% verification de la non modification du faust
[new_dim1,new_dim2]=size(F); 
if ((new_dim1 ~= dim1) | (dim2 ~= new_dim2))
    error(['transpose : modification du faust de depart']);
end

new_F_dense=full(F);
%((new_F_dense ~= F_dense))
if (~isequal(new_F_dense,F_dense))
	error('transpose : modification du faust de depart');
end 


disp('operation F_trans_trans = F_trans'' ');
F_trans_trans=F_trans.';
[dim1_trans_trans,dim2_trans_trans]=size(F_trans_trans);

if ((dim1_trans_trans ~= dim1) | (dim2_trans_trans ~= dim2))
    error(['transpose : invalid dimension']);
end

F_dense_trans_trans = full(F_trans_trans);
%(F_dense_trans_trans ~= F_dense)  
if (~isequal(F_dense_trans_trans,F_dense))
    error(['transpose : invalid transpose matrix']);
end


%% verification de la non modification du faust
[new_dim1_trans,new_dim2_trans]=size(F_trans); 
if ((new_dim1_trans ~= dim1_trans) | (new_dim2_trans ~= new_dim2_trans))
    error(['transpose : modification du faust de depart']);
end

new_F_dense_trans=full(F_trans);
%((new_F_dense_trans ~= F_dense_trans))
if (~isequal(new_F_dense_trans,F_dense_trans))
	error('transpose : modification du faust de depart');
end 





disp('Ok');




%%% ctranspose test
if (isreal(F))
	disp('TEST CTRANSPOSE : ');
	disp('operation F_ctrans = F'' ');
	F_ctrans=F';
	[dim1_ctrans,dim2_ctrans]=size(F_ctrans);
	if ((dim1_ctrans ~= dim2) | (dim2_ctrans ~= dim1))
	    error(['transpose : invalid dimension']);
	end

	F_dense_ctrans = full(F_ctrans);

	%(F_dense_ctrans ~= F_dense')
	if (~isequal(F_dense_ctrans,F_dense'))    
		error(['transpose : invalid transpose matrix']);
	end

	%% verification de la non modification du faust
	[new_dim1,new_dim2]=size(F); 
	if ((new_dim1 ~= dim1) | (dim2 ~= new_dim2))
	    error(['transpose : modification du faust de depart']);
	end

	new_F_dense=full(F);
	%((new_F_dense ~= F_dense))
	if (~isequal(new_F_dense,F_dense))
		error('transpose : modification du faust de depart');
	end 


	disp('operation F_ctrans_ctrans = F_ctrans'' ');
	F_ctrans_ctrans=F_ctrans';
	[dim1_ctrans_ctrans,dim2_ctrans_ctrans]=size(F_ctrans_ctrans);

	if ((dim1_ctrans_ctrans ~= dim1) | (dim2_ctrans_ctrans ~= dim2))
	    error(['transpose : invalid dimension']);
	end

	F_dense_ctrans_ctrans = full(F_ctrans_ctrans);
	%(F_dense_trans_trans ~= F_dense)  
	if (~isequal(F_dense_ctrans_ctrans,F_dense))
	    error(['transpose : invalid transpose matrix']);
	end


	%% verification de la non modification du faust
	[new_dim1_ctrans,new_dim2_ctrans]=size(F_ctrans); 
	if ((new_dim1_ctrans ~= dim1_ctrans) | (new_dim2_ctrans ~= new_dim2_ctrans))
	    error(['transpose : modification du faust de depart']);
	end

	new_F_dense_ctrans=full(F_ctrans);
	%((new_F_dense_ctrans ~= F_dense_ctrans))
	if (~isequal(new_F_dense_ctrans,F_dense_ctrans))
		error('transpose : modification du faust de depart');
	end 





	disp('Ok');
end






























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
disp('TEST MULTIPLICATION BY A REAL VECTOR : ');
x=zeros(dim2,1);
x(:)=1:dim2;
x_trans=zeros(dim1,1);
x_trans(:)=1:dim1;

test_matlab_faust_mult(F,F_dense,x,x_trans);
disp('Ok');


%% test multiplication with matrix
disp('TEST MULTIPLICATION BY A REAL MATRIX : ');
X=zeros(dim2,dim3);
X(:)=1:dim2*dim3;
X_trans=zeros(dim1,dim3);
X_trans(:)=1:dim1*dim3;


test_matlab_faust_mult(F,F_dense,X,X_trans);
disp('Ok');






disp('TEST MULTIPLICATION BY A COMPLEX VECTOR : ');
x_cplx = randi(100,dim2,1) + 1i * randi(100,dim2,1);
x_cplx_trans = randi(100,dim1,1) + 1i * randi(100,dim1,1);

test_matlab_faust_mult(F,F_dense,x_cplx,x_cplx_trans);
disp('Ok');


disp('TEST MULTIPLICATION BY A COMPLEX MATRIX : ');
X_cplx = randi(100,dim2,dim3) + 1i * randi(100,dim2,dim3);
X_cplx_trans = randi(100,dim1,dim3) + 1i * randi(100,dim1,dim3);

test_matlab_faust_mult(F,F_dense,X_cplx,X_cplx_trans);
disp('Ok');








%% test 2-norm (spectral norm)
disp('TEST 2-norm : ');
real_norm=norm(F_dense);
norm_faust=norm(F);
norm_faust2=norm(F,2);
norm_faust_trans=norm(F_trans);

if  ( (abs(real_norm - norm_faust)/abs(real_norm)) > threshold)
	error(['norm : invalid result, expected ' num2str(real_norm) ' get norm_faust' num2str(norm_faust)]);
end
if (norm_faust ~= norm_faust2)
	error(['norm : norm(F) must be equal to norm(F,2)']);
end

if (norm_faust_trans ~= norm_faust)
	error(['norm : norm(F) must be equal to norm(F_trans)']);
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
	if(~isequal(A,factors{i}))
        get_fact(F,i)
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

F_dense_trans_loaded=full(F_trans_loaded);
if (~isequal(F_dense_trans_loaded,F_dense.'))
    disp('F_dense_trans_loaded')
    F_dense_trans_loaded
    disp('F_dense_trans=')
    F_dense.'
    disp(F_dense_trans_loaded==F_dense.')
    error(['save transposed : invalid faust']);
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

%% test conj
disp('TEST CONJ : ');
full(F)
F_conj = conj(F)
F_conj_full= full(F_conj)
expected_F_conj_full = conj(full(F))
[dim1,dim2]=size(expected_F_conj_full);
[dim1_conj,dim2_conj]=size(F_conj_full);

if((dim1_conj ~= dim1) | (dim2_conj ~= dim2))
    dim1
    dim2
    dim1_conj
    dim2_conj
    error('full : invalid dimension');
end
if ( ~isequal(expected_F_conj_full,F_conj_full) )
    expected_F_conj_full
    F_conj_full
    error(['conj test 1 failed.']);
end
% test get_fact on conj
for i=1:nb_fact
	A=get_fact(F_conj,i);
	if(~isequal(A,conj(factors{i})))
        get_fact(F_conj,i)
		error('get_fact : invalid factor');
	end
end
% test conj save
save(conj(F),filename)
saved_conj_F=full(Faust(filename))
if ( ~isequal(saved_conj_F,F_conj_full))
    saved_conj_F
    F_conj_full
    error(['conj test 3 failed.']);
end
disp('Ok');


%% test ctranspose
disp('TEST CTRANSPOSE : ');
full(F)
F_ctranspose = ctranspose(F)
F_ctranspose_full= full(F_ctranspose)
expected_F_ctranspose_full = ctranspose(full(F))
[dim1,dim2]=size(expected_F_ctranspose_full);
[dim1_ctranspose,dim2_ctranspose]=size(F_ctranspose_full);

if((dim1_ctranspose ~= dim1) | (dim2_ctranspose ~= dim2))
    dim1
    dim2
    dim1_ctranspose
    dim2_ctranspose
    error('full : invalid dimension');
end
if ( ~isequal(expected_F_ctranspose_full,F_ctranspose_full) )
    expected_F_ctranspose_full
    F_ctranspose_full
    error(['ctranspose test 1 failed.']);
end
% test get_fact on ctranspose
for i=1:nb_fact
	A=get_fact(F_ctranspose,i)
	if(~isequal(A,ctranspose(factors{nb_fact-i+1})))
        get_fact(F_ctranspose,i)
		error('get_fact : invalid factor');
	end
end
% test ctranspose save
save(ctranspose(F),filename)
saved_ctranspose_F=full(Faust(filename))
if ( ~isequal(saved_ctranspose_F,F_ctranspose_full))
    saved_ctranspose_F
    F_ctranspose_full
    error(['ctranspose test 3 failed.']);
end
disp('Ok');








