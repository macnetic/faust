close all;
cpp_test_path = [getenv('FAUST_ROOT_DIR') '/trunk/devcpp/palm4MSA/test'];
M=1:12;
N_run = 10000;
n = 2.^M;


%densite_or_nnzperroworcol = ones(size(M))*0.5;
densite_or_nnzperroworcol = ones(size(M))*1;

for k = numel(M):-1:1
    if densite_or_nnzperroworcol(k) > n(k)
       M(k)=[];
       n(k)=[];
    end
end

matrix_filename = sprintf('%s/Matrix_test_%d-%d_%d.mat', cpp_test_path,n(1),n(end),N_run);
matCells = demo_compTime_matrix_draw(M, n, densite_or_nnzperroworcol, matrix_filename);
[t_faust_matlab, t_dense_matlab]=demo_compTime(M, N_run, n, matrix_filename, matCells);

setenv('LD_LIBRARY_PATH', [ ...
    '/home/tgautrai/local/OpenBLAS-0.2.14:' ...
    '/home/intel/composer_xe_2013.2.146/mkl/lib/intel64:'...
    '/home/intel/composer_xe_2013.2.146/mkl/../compiler/lib/intel64:'...
    '/usr/local/cuda-6.5/lib64' ...
    ]);
system('cd $FAUST_ROOT_DIR/trunk/devcpp/palm4MSA/test && make comptime0');
[a,b]=system(['cd $FAUST_ROOT_DIR/trunk/devcpp/palm4MSA/test && ./comptime0 ' matrix_filename ' ' num2str(N_run)]);
ind_LF = strfind(b, sprintf('\n'));
if length(ind_LF)~=2, error('incorrect number of character LF'); end
t_dense_cpp=load(b(1:ind_LF(1)-1));
t_faust_cpp=load(b(ind_LF(1)+1:ind_LF(2)-1));

if max( size(t_dense_cpp)    ~= size(t_dense_matlab) | ...
        size(t_dense_matlab) ~= size(t_faust_matlab) | ...
        size(t_faust_matlab) ~= size(t_faust_cpp))
    error('dimensions differentes');
end

figures_compTime(M, n, densite_or_nnzperroworcol, t_dense_matlab, t_faust_matlab, t_dense_cpp, t_faust_cpp);

disp('');