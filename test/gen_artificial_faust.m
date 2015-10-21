%% script pour generer des Faust avec un RCG, une Dimension et un nombre de facteur fixe
% lancer l'executable comptime1 ensuite pour effectuer les tests de performances 
% puis drawComptime1 pour afficher les resultats
RCGs = [1 2 3 4 5];
Dims = [128 256 512];
nb_facts = [2,4,6];
nb_test = 50;

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
            for k=1:nb_fact
                facts{k} = sprand(dim1,dim2,densite_per_fact);
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
save('data/Faust_example.mat');

