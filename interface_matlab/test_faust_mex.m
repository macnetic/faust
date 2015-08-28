% Si k>1,            DIMS(k) correspond au nombre de colonnes de facts{k-1}
% Si k<length(DIMS), DIMS(k) correspond au nombre de  lignes  de facts{k}
DIMS    = [5000 2100 2140 2004];

% DENSITY(k) correspond a la densite de facts{k}
DENSITY = [0.01 0.004 0.007];


if(length(DENSITY) ~= length(DIMS)-1)
    error('Lengths of array DIMS or DENSITY don''t match');
end

facts_sparse = cell(1,(length(DENSITY)));
facts = cell(size(facts_sparse));
for k=1:length(facts_sparse)
    facts_sparse{k} = sprand(DIMS(k), DIMS(k+1),DENSITY(k));
    facts{k} = full(facts_sparse{k}); 
end


cd ([getenv('FAUST_ROOT_DIR') '/trunk/Code_Luc']);
set_path;
cd ([getenv('FAUST_ROOT_DIR') '/trunk/devcpp/interface_matlab'] );
objectHandle = mexLoadFaust(facts);

v=rand(size(facts{end},2),1);

tic;
w_faust_mex=faust_mex('multiply', objectHandle, v);
t_faust_mex = toc;
fprintf('temps multiplication avec faust mex : %g s\n',t_faust_mex);toc

for k=1:length(facts),facts_sparse{k}=sparse(facts_sparse{k});end
w_faust_matlab=v;
tic;
for k=length(facts):-1:1, w_faust_matlab=facts_sparse{k}*w_faust_matlab;end
t_faust_matlab = toc;
fprintf('temps multiplication avec faust matlab : %g s\n',t_faust_matlab);toc

M = dvp(facts);
tic
w_dense_matlab=M*v;
t_dense_matlab = toc;
fprintf('temps multiplication avec dense matlab : %g s\n\n',t_dense_matlab);


fprintf('erreur max faust_mex - faust_matlab : %g \n',max(max(abs(w_faust_matlab-w_faust_mex))));
fprintf('erreur max faust_mex - dense_matlab : %g \n\n',max(max(abs(w_dense_matlab-w_faust_mex))));

fprintf('rapport temps faust matlab / temps faust mex : %g \n',t_faust_matlab/t_faust_mex);
fprintf('rapport temps dense matlab / temps faust mex : %g \n',t_dense_matlab/t_faust_mex);