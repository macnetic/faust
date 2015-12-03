%% script pour generer des Faust avec un RCG, une Dimension et un nombre de facteur fixe
% lancer l'executable comptime1 ensuite pour effectuer les tests de performances 
% puis drawComptime1 pour afficher les resultats
function gen_artficial_faust(constraint)
%% cas 1
% RCGs = [2 4 6 8 10];
% Dims = [128 256 512];
% nb_facts = [2,4,6];

%% cas 2

RCGs = [4 8 16 32 64];
Dims = [128 256 512 2048];

nb_facts = [2,4,8,16];

% RCGs = [2 4];
% Dims = [32 64 128];
% nb_facts = [2,4];

nb_test = 50;
% constraint='sp_row';% choix possible 'sp_row' et 'sp_col_'


nRCGS = length(RCGs);
nDims = length(Dims);
nfacts = length(nb_facts);

real_RCGs = zeros(nb_test,nRCGS,nDims,nfacts);
ts_dense = zeros(nb_test,nRCGS,nDims,nfacts);
ts_mex = zeros(nb_test,nRCGS,nDims,nfacts);
ts_random_dense = zeros(nb_test,nRCGS,nDims,nfacts);
opt_calcul = 0;


densite = RCGs.^(-1);


addpath('../build/mex')
addpath('./tools');

if ((strcmp(constraint,'sp') + strcmp(constraint,'sp_row') + strcmp(constraint,'sp_col')) == 0)
    error('the constraint parameter must be equal to sp or sp_row or sp_col');
end


for h=1:nfacts
    nb_fact=nb_facts(h);
    disp(['NB_FACT : ' int2str(nb_fact)]);
    for j=1:nDims
        dim1=Dims(j);
        dim2=dim1;
        disp(['  Dim : ' int2str(dim1)]);
        random_dense_mat = rand(dim1,dim2);
        eval(['DMat_Dim_' int2str(dim1) ' = random_dense_mat;']);
        fprintf('RCG : ');
        for i=1:nRCGS
            RCG=RCGs(i);
            fprintf('%d ',RCG);
            densite_per_fact = densite(i)/nb_fact;
            nnz = densite_per_fact*Dims(j)^2;
            facts=cell(1,nb_fact);
             if (strcmp(constraint,'sp_row'))
                 nb_elt_per_row = round(dim2*densite_per_fact);
                 nb_elt = nb_elt_per_row * dim1;
                 id_i=reshape(repmat((1:dim1),nb_elt_per_row,1),dim1*nb_elt_per_row,1);
                 id_j=zeros(nb_elt,1);
                 value=zeros(nb_elt,1);
                 
             end
             if (strcmp(constraint,'sp_col'))
                 nb_elt_per_col = round(dim1*densite_per_fact);
                 nb_elt = nb_elt_per_col * dim2;
                 id_j=reshape(repmat((1:dim2),nb_elt_per_col,1),dim2*nb_elt_per_col,1);
                 id_i=zeros(nb_elt,1);
                 value=zeros(nb_elt,1);
             end
            for k=1:nb_fact
                if (strcmp(constraint,'sp'))
                    facts{k} = sprand(dim1,dim2,densite_per_fact);
                end
                 if (strcmp(constraint,'sp_row'))
                     value=rand(nb_elt,1);
                     for ll=0:dim1-1
                         id_j(ll*nb_elt_per_row+1:(ll+1)*nb_elt_per_row)=randperm(dim2,nb_elt_per_row)';
                     end
                     facts{k}=sparse(id_i,id_j,value,dim1,dim2);
                 end
                 if (strcmp(constraint,'sp_col'))
                     value=rand(nb_elt,1);
                     for ll=0:dim2-1
                         id_i(ll*nb_elt_per_col+1:(ll+1)*nb_elt_per_col)=randperm(dim1,nb_elt_per_col)';
                     end
                     facts{k}=sparse(id_i,id_j,value,dim1,dim2);
                 end
            end
            
            S_prod = facts{1};
            for k=2:nb_fact
                S_prod = S_prod * facts{k};
            end
            S_prod = full(S_prod);
            
            
            eval(['Faust_nfact_' int2str(nb_fact) '_RCG_' int2str(RCG) '_Dim_' int2str(dim1) ' = facts;']);
            if (opt_calcul ~= 0)
                fc = matlab_faust(facts);
                for k=1:nb_test
                    v=rand(dim2,1);
                    
                    %% MEX
                    tic
                    w_faust_mex = fc*v; % surcharge de l'operateur *  de la class matlab_faust
                    t_mex = toc;
                    
                    
                    %% DENSE
                    w_faust_matlab=v;
                    tic;
                    w_matlab = S_prod * v;
                    t_dense = toc;
                    
                    %% RANDOM
                    w_faust_matlab=v;
                    tic;
                    w_matlab_random = random_dense_mat * v;
                    t_random_dense = toc;
                    
                    
                    ts_mex(k,i,j,h) = t_mex;
                    ts_dense(k,i,j,h) = t_dense;
                    ts_random_dense(k,i,j,h) = t_random_dense;
                    
                    
                end
            end
            
            
            
        end
        fprintf('\n');
        
        
    end
end
if (opt_calcul ~= 0)
    mean_t_random_dense=squeeze(mean(ts_random_dense));
    mean_t_dense=squeeze(mean(ts_dense));
    mean_t_mex=squeeze(mean(ts_mex));
    rapport_dense = mean_t_random_dense ./  mean_t_dense;
    real_RCG = mean_t_dense ./  mean_t_mex;
    disp([' comp_dense vs random_dense ']);
    for k=1:nfacts
        disp( ['nb_fact = ' int2str(nb_facts(k))]);
        for i=1:nRCGS
            disp(num2str(rapport_dense(i,:,k)));
        end
        disp('');
        
    end
    disp('');
    disp('');
    disp('');
    disp(['real RCG ']);
    for i=1:nfacts
        disp( ['nb_fact = ' int2str(nb_facts(k))]);
        disp(num2str([0 Dims]));
        for i=1:nRCGS
            disp(num2str([ RCGs(i) real_RCG(i,:,k)]));
        end
    end
end
disp('sauvegarde des donnees');
save(['data/Faust_example_' constraint '.mat']);                     


end
