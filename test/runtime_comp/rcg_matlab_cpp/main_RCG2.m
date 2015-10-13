%idem que main_RCG.m sauf qu'au lieu de faire N_run sur les memes matrices,
%on fait N_run sur des matrices differentes
close all;
cpp_test_path = [getenv('FAUST_ROOT_DIR') '/trunk/devcpp/palm4MSA/test'];

% M : nbre de facteurs ; M doit etre soit un scalaire, soit un vecteur LIGNE
M=1:7;
% n : dimension de la matrice ; n doit avoir les memes dimensions que M
n = 2.^M;

N_run = 20;



% on boucle le programme sur le nombre de lignes de densite_or_nnzperroworcol
% et a l iteration i, tous les M(k) facteurs auront une densite (ou un nnz par col/ligne)
% de densite_or_nnzperroworcol(i,k)

% densite_or_nnzperroworcol : densite ou nnz par ligne/colonne
% densite_or_nnzperroworcol doit etre une matrice (ou un vecteur) avec le meme nmbre de
% colonne que M ou un vecteur colonne ou un scalaire
densite_or_nnzperroworcol = [2;4;8];

% le nombre de colonne de densite_or_nnzperroworcol doit etre egal a celui
% de M
if(size(densite_or_nnzperroworcol,2)==1)
    densite_or_nnzperroworcol = repmat(densite_or_nnzperroworcol,size(M));
elseif(size(densite_or_nnzperroworcol,2) ~= size(M,2))
    error('dimension de densite_or_nnzperroworcol incorrectes');
end

% on verifie qu on ne fait pas deux fois la meme simu


% setenv('LD_LIBRARY_PATH', [ ...
%     getenv('OPENBLAS_ROOT_DIR') ':' ...
%     getenv('MKLROOT') '/lib/intel64:'...
%     getenv('MKLROOT') '/../compiler/lib/intel64:'...
%     getenv('CUDADIR') '/lib64' ...
%     ]);
system('cd $FAUST_ROOT_DIR/trunk/devcpp/palm4MSA/test && make comptime0');

for ind_nz=1:size(densite_or_nnzperroworcol,1)
    t1=tic;
    disp(['ind_nz = ' num2str(ind_nz) ' sur ' num2str(size(densite_or_nnzperroworcol,1))])
    for k = numel(M):-1:1
        % si densite_or_nnzperroworcol correspond a nz
        if densite_or_nnzperroworcol(ind_nz,:) > 1
            if densite_or_nnzperroworcol(ind_nz,:) > n(k)
                M(k)=[];
                n(k)=[];
            end
            % si densite_or_nnzperroworcol correspond a densite
        else
            % si la matrice creuse est vide, on ne considere pas ce cas
            if round(densite_or_nnzperroworcol(ind_nz,:)*n(k)) < 2
                M(k)=[];
                n(k)=[];
            end
        end
    end
    
    t_faust_matlab = NaN(N_run, numel(M));
    t_dense_matlab = NaN(N_run, numel(M));
    t_faust_cpp = NaN(N_run, numel(M));
    t_dense_cpp = NaN(N_run, numel(M));
    
    
    if min(densite_or_nnzperroworcol(ind_nz,:) == densite_or_nnzperroworcol(ind_nz,1))
      matrix_filename = sprintf('%s/Matrix_test_%d-%d_%d_%g.mat', cpp_test_path,n(1),n(end),N_run,densite_or_nnzperroworcol(ind_nz,1));
    else
      matrix_filename = sprintf('%s/Matrix_test_%d-%d_%d_densitevariable.mat', cpp_test_path,n(1),n(end),N_run);
    end  
    for r=1:N_run
        matCells = demo_compTime_matrix_draw(M, n, densite_or_nnzperroworcol(ind_nz,:), matrix_filename);
        [t_faust_matlab(r,:), t_dense_matlab(r,:)]=demo_compTime(M, 1, n, matrix_filename, matCells);
        
        [a,b]=system(['cd $FAUST_ROOT_DIR/trunk/devcpp/palm4MSA/test && ./comptime0 ' matrix_filename ' 1']);
        ind_LF = strfind(b, sprintf('\n'));
        if length(ind_LF)~=2, error('incorrect number of character LF'); end
        t_dense_cpp(r,:) = load(b(1:ind_LF(1)-1));
        t_faust_cpp(r,:) = load(b(ind_LF(1)+1:ind_LF(2)-1));
        fprintf('.');
    end
    fprintf('\n');
    
    figures_compTime(M, n, densite_or_nnzperroworcol(ind_nz,:), t_dense_matlab, t_faust_matlab, t_dense_cpp, t_faust_cpp);
    %saveas(gcf,strrep(matrix_filename,'.mat','.fig'));
    %close
    %delete(matrix_filename)
    toc(t1)
    fprintf('\n');
end