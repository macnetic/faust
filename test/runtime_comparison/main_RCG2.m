%idem que main_RCG.m sauf qu'au lieu de faire N_run sur les memes matrices,
%on fait N_run sur des matrices differentes
clear
%close all;

% % M : nbre de facteurs ; M doit etre soit un scalaire, soit un vecteur
% LIGNE dont les element sont de valeur croissante
M=1:10;
%M = [12 2 3 4 5 13];
%M = [2 3 4 6 8 12];
% n : dimension de la matrice ; n doit etre un scalaire ou un vecteur de memes dimensions que M
n = 2.^M;
%n = M./M*64;

% nbre de tests effectues avec des matrices differents
N_run = 10;

%nbre de tests effectues avec les memes matrices/home/tgautrai/faust/trunk/devcpp/palm4MSA/test/runtime_comp/output/main_RCG2/test5/test5.2
N_run_inside = 100;

% si matlab_only=1, le test est fait uniquement avec Matlab
% si matlab_only=0, le test est fait  a la fois avec matlab et en C++
matlab_only = 0;

if ~matlab_only
    tmp_test_path = '@FAUST_PALM4MSATESTRUNCOMP_BIN_DIR@/tmp';
    if(~exist(tmp_test_path,'dir'))
        mkdir(tmp_test_path)
    end
end


% on boucle le programme sur le nombre de lignes de densite_or_nnzperroworcol
% et a l iteration i, tous les M(k) facteurs auront une densite (ou un nnz par col/ligne)
% de densite_or_nnzrowcol(i,k)

% densite_or_nnzrowcol : densite(si <=1) ou nnz(si >1) par ligne/colonne
% densite_or_nnzrowcol doit etre une matrice (ou un vecteur) avec le meme nmbre de
% colonne que M ou un vecteur colonne ou un scalaire
densite_or_nnzrowcol = 2;
%densite_or_nnzrowcol = [2 ;3; 4 ;6 ;8 ;12];

if(length(M)==1 && size(densite_or_nnzrowcol,2)>1), M=repmat(M,1,size(densite_or_nnzrowcol,2)); end
if(length(M)==1 && length(n)>1), M=repmat(M,size(n)); end

if(length(n)==1 && size(densite_or_nnzrowcol,2)>1), n=repmat(n,1,size(densite_or_nnzrowcol,2)); end
if(length(n)==1 && length(M)>1), n=repmat(n,size(M)); end


%%
% le nombre de colonne de densite_or_nnzrowcol doit etre egal a celui
% de M
if(size(densite_or_nnzrowcol,2)==1)
    densite_or_nnzrowcol = repmat(densite_or_nnzrowcol,size(M));
elseif(size(densite_or_nnzrowcol,2) ~= size(M,2))
    error('dimension de densite_or_nnzrowcol incorrecte');
end

% on verifie qu on ne fait pas deux fois la meme simu
for k = size(densite_or_nnzrowcol,1):-1:1
    for l = 1:k-1
        if(min(densite_or_nnzrowcol(l,:)==densite_or_nnzrowcol(k,:)))
            densite_or_nnzrowcol(k,:)=[];
            break;
        end
    end
end



% ATTENTION CONVERSION DE densite_or_nnzrowcol DE MATRICE EN CELLULES
densite_or_nnzrowcol = mat2cell(densite_or_nnzrowcol,ones(1,size(densite_or_nnzrowcol,1)),length(M));

if ~matlab_only
    setenv('LD_LIBRARY_PATH', getenv('OPENBLAS_ROOT_DIR'));
    system('cd @PROJECT_BINARY_DIR@ && make comptime0');
end


for ind_nz=1:length(densite_or_nnzrowcol)
    t1=tic;
    disp(['ind_nz = ' num2str(ind_nz) ' sur ' num2str(length(densite_or_nnzrowcol))])
    for k = numel(M):-1:1
        % si densite_or_nnzrowcol correspond a nz
        if densite_or_nnzrowcol{ind_nz}(k) > 1
            if densite_or_nnzrowcol{ind_nz}(k) > n(k)
                M(k)=[];
                n(k)=[];
                densite_or_nnzrowcol{ind_nz}(k)=[];
            end
            % si densite_or_nnzrowcol correspond a densite
        else
            % si la matrice creuse est vide ou ne contient qu un element
            % par ligne ou par colonne, on ne considere pas ce cas
            if round(densite_or_nnzrowcol{ind_nz}(k)*n(k)) < 2
                M(k)=[];
                n(k)=[];
                densite_or_nnzrowcol{ind_nz}(k)=[];
            end
        end
    end
    
    t_faust_matlab = NaN(N_run, numel(M));
    t_dense_matlab = NaN(N_run, numel(M));
    t_faust_cpp = NaN(N_run, numel(M));
    t_dense_cpp = NaN(N_run, numel(M));
    
    
    if min(densite_or_nnzrowcol{ind_nz} == densite_or_nnzrowcol{ind_nz}(1))
        if matlab_only
            matrix_filename = sprintf('%s/Matrix_test_%d-%d_%d_%g.mat', tmp_test_path, n(1),n(end),N_run,densite_or_nnzrowcol{ind_nz}(1));
        else
            matrix_filename = sprintf('%s/Matrix_test_%d-%d_%d_%g.mat', tmp_test_path,n(1),n(end),N_run,densite_or_nnzrowcol{ind_nz}(1));
        end
    else
        if matlab_only
            matrix_filename = sprintf('%s/Matrix_test_%d-%d_%d_densitevariable.mat', tmp_test_path, n(1),n(end),N_run);
        else
            matrix_filename = sprintf('%s/Matrix_test_%d-%d_%d_densitevariable.mat', tmp_test_path,n(1),n(end),N_run);
        end
    end
    cppOutputFile1 = '';
    cppOutputFile2 = '';
    for r=1:N_run
        [matCells,status] = demo_compTime_matrix_draw(M, n, densite_or_nnzrowcol{ind_nz}, matrix_filename);
        [t_faust_matlab(r,:), t_dense_matlab(r,:)]=demo_compTime(M, N_run_inside, n, matrix_filename, matCells, status);
        if ~matlab_only
            [a,b]=system(['cd $FAUST_ROOT_DIR/trunk/devcpp/build/testing/bin && ./comptime0 ' matrix_filename ' ' num2str(N_run_inside) ' ' sprintf('%d_',status)]);
        end
        delete(matrix_filename);
        if ~matlab_only
            ind_LF = strfind(b, sprintf('\n'));
            if length(ind_LF)~=2, error('incorrect number of character LF'); end
            cppOutputFile1 = b(1:ind_LF(1)-1);
            cppOutputFile2 = b(ind_LF(1)+1:ind_LF(2)-1);
            t_dense_cpp(r,:) = mean(load(cppOutputFile1),1);
            t_faust_cpp(r,:) = mean(load(cppOutputFile2),1);
            fprintf('temps matlab faust / temps cpp faust = %g ',t_faust_matlab(r,:)/t_faust_cpp(r,:));
            fprintf('\n');
        end
        %fprintf('.');
    end
    if ~isempty(cppOutputFile1),delete(cppOutputFile1);end
    if ~isempty(cppOutputFile2),delete(cppOutputFile2);end
    fprintf('\n');
    
    f1=figures_compTime(M, n, densite_or_nnzrowcol{ind_nz}, N_run_inside, t_dense_matlab, t_faust_matlab, t_dense_cpp, t_faust_cpp);
    param.M = M;
    param.n = n;
    param.N_run = N_run;
    param.N_run_inside = N_run_inside;
    param.densite_or_nnzrowcol = densite_or_nnzrowcol{ind_nz};
    setappdata(f1,'param',param);
    [pathstr,name,ext] = fileparts(matrix_filename);
    pathstr = [pathstr filesep '..' filesep 'figures'];
    if ~exist(pathstr,'dir'), mkdir(pathstr);end
    saveas(gcf,[pathstr filesep name '.fig']);
    close(f1);
    toc(t1)
    fprintf('\n');
end







