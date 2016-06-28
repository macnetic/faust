%% Description BSL
%  Brain source localization
%
%  This script performs brain source localization using several gain
%  matrices [2], including FAuSTs, and several solvers. It reproduces the
%  source localization experiment of [1].
%  The results are stored in "./output/results_BSL_user.mat".
%  DURATION: Computations should take around 10 minutes. 
%	
%  The MEG gain matrices used are 
%		- those in "./output/M_user.mat" if available
%		- or the precomputed ones in "./precomputed_results_MEG/M_X.mat"
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
%
%% License:
% Copyright (2016):	Luc Le Magoarou, Remi Gribonval
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
%   Nicolas Bellot : nicolas.bellot@inria.fr
%   Leman Adrien   : adrien.leman@inria.fr
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%
% [2]   A. Gramfort, M. Luessi, E. Larson, D. Engemann, D. Strohmeier, 
%	C. Brodbeck, L. Parkkonen, M. Hamalainen, MNE software for processing
%	MEG and EEG data <http://www.ncbi.nlm.nih.gov/pubmed/24161808>, 
%	NeuroImage, Volume 86, 1 February 2014, Pages 446-460, ISSN 1053-8119, 
%	[DOI] <http://dx.doi.org/10.1016/j.neuroimage.2013.10.027>
%%



%% This is an example of a little experience of source localization in magnetoencephalography (MEG) using faust factorisation for speed-up the calculus,
%% the sparse algorithm used is greed_omp_chol from the toolbox smallbox2.0/Sparsify/GreedLab under GNU GPL 2.0 License
%% to show the different performance figure run figure_localisze_src after this script


runPath=which(mfilename);
pathname = fileparts(runPath);
BSL_data_pathName=strcat(pathname,'/data/');




RCG_approxS_MEG=[6,8,16,25];
nb_approx_MEG = length(RCG_approxS_MEG);
MEG_approxS_norm = cell(1,nb_approx_MEG);


%Loading of the MEG matrix
load([BSL_data_pathName 'X_meg.mat' ]);
load ([BSL_data_pathName 'curv.mat']);
points2 = points;
points = points(points_used_idx,:);
X = X_fixed;
%X = randn(size(X_fixed));
%X_norm = X;
X_norm = X./repmat(sqrt(sum(X.^2,1)),size(X,1),1);

%Loading of the MEG matrix approximations
MEG_faustS=cell(1,nb_approx_MEG);
faustS_mult=cell(1,nb_approx_MEG);
trans_faustS_mult=cell(1,nb_approx_MEG);
matlab_faustS_mult=cell(1,nb_approx_MEG);
matlab_trans_faustS_mult=cell(1,nb_approx_MEG);
for i=1:nb_approx_MEG
    

    
    RCG_approxS_MEG(i);
    load([BSL_data_pathName 'M_' int2str(RCG_approxS_MEG(i))  ]);
    facts{1}=lambda*facts{1};
    X_approx =dvp(facts);
    %X_hat = X' + 0.05*randn(size(X'));
    X_approx = X_approx'; X_approx(:,sum(X_approx.^2,1)==0)=1;
    %X_hat_norm = X_hat;
    X_norm_approx = X_approx./repmat(sqrt(sum(X_approx.^2,1)),size(X_approx,1),1);
    MEG_approxS_norm{i}=X_norm_approx;
    
    %% structure faust
    facts{1}=facts{1}./repmat(sqrt(sum(X_approx.^2,1)),size(X_approx,1),1)';% normalisation of the row dvp(facts)
    sp_facts_trans = make_sparse(facts);
    sp_facts = faust_transpose(sp_facts_trans);
    trans_fc=matlab_faust(facts);
    fc=transpose(trans_fc);
%     f=@(x) x/2;
    matlab_trans_faustS_mult{i}=@(x) f_mult(sp_facts_trans,x);
    matlab_faustS_mult{i}=@(x) f_mult(sp_facts,x);
    trans_faustS_mult{i}=@(x) trans_fc*x;
    faustS_mult{i}=@(x) fc*x;
    MEG_faustS{i}=trans_fc;
%     eval(['function y = fc_mult' int2str(i) '(x) y=fc*x;end']);
%     eval(['foncteur{i}=@fc_mult' int2str(i) ';']);

end



Ntraining = 500; % Number of training vectors
Sparsity = 2; % Number of sources per training vector
dist_paliers = [0.01,0.05,0.08]; dist_paliers = [dist_paliers, 0.5];

resDist = zeros(nb_approx_MEG+1,numel(dist_paliers)-1,Sparsity,Ntraining); % (Matrice,m�thode,dist_sources,src_nb,run);
compute_Times = zeros(nb_approx_MEG+1,numel(dist_paliers)-1,Ntraining);
resDist_matlab = zeros(nb_approx_MEG+1,numel(dist_paliers)-1,Sparsity,Ntraining); % (Matrice,m�thode,dist_sources,src_nb,run);
compute_Times_matlab = zeros(nb_approx_MEG+1,numel(dist_paliers)-1,Ntraining);
for k=1:numel(dist_paliers)-1
    disp(['k=' num2str(k) '/' num2str(numel(dist_paliers)-1)])
    %Parameters settings
    Gamma = zeros(size(X_norm,2),Ntraining);
    for ii=1:Ntraining
        %disp(num2str(ii));
        dist_sources = -1;
        while ~((dist_paliers(k)<dist_sources)&&(dist_sources<dist_paliers(k+1)))
            Gamma(:,ii) = sparse_coeffs(X_norm, 1, Sparsity); %
            idx = find(Gamma(:,ii));
            dist_sources = norm(points(idx(1)) - points(idx(2)));
        end
        % dist_sources = norm(points(idx(1)) - points(idx(2)))
    end
    Data = X_norm*Gamma;
    
    
    
    sol_omp = zeros(size(Gamma));
    sol_omp_hat = zeros(size(Gamma));
    sol_omp_hat2 = zeros(size(Gamma));
    err_omp = zeros(2,Ntraining);
    err_omp_hat = zeros(2,Ntraining);
    err_omp_hat2 = zeros(2,Ntraining);
    diff_omp = zeros(2,Ntraining);
    dist_omp = zeros(Sparsity,Ntraining);
    dist_omp_hat = zeros(Sparsity,Ntraining);
    dist_omp_hat2 = zeros(Sparsity,Ntraining);
    
    
    
    for i=1:Ntraining
        disp(['   i=' num2str(i) '/' num2str(Ntraining)])
        idx = find(Gamma(:,i));
        dist_sources = norm(points(idx(1)) - points(idx(2)));
        
        %l1 solving
        tol = 1e-4;
        lambda = 0.3;
                


        %OMP solving
            tic
            [sol_omp(:,i), err_mse_omp, iter_time_omp]=greed_omp_chol(Data(:,i),X_norm,size(X_norm,2),'stopTol',1*Sparsity);
            t1=toc;
            err_omp(1,i) = norm(X_norm*Gamma(:,i)-X_norm*sol_omp(:,i))/norm(X_norm*Gamma(:,i));
            err_omp(2,i) = isequal(find(Gamma(:,i)),find(sol_omp(:,i)>1e-4));
            idx_omp = find(sol_omp(:,i));
            resDist(1,k,1,i) = min(norm(points(idx(1)) - points(idx_omp(1))),norm(points(idx(1)) - points(idx_omp(2))));
            resDist(1,k,2,i) = min(norm(points(idx(2)) - points(idx_omp(1))),norm(points(idx(2)) - points(idx_omp(2))));
            compute_Times(1,k,i)=t1;
            resDist_matlab(1,k,1,i) = min(norm(points(idx(1)) - points(idx_omp(1))),norm(points(idx(1)) - points(idx_omp(2))));
            resDist_matlab(1,k,2,i) = min(norm(points(idx(2)) - points(idx_omp(1))),norm(points(idx(2)) - points(idx_omp(2))));
            compute_Times_matlab(1,k,i)=t1;
            
       for ll=1:nb_approx_MEG
              X_approx_norm = MEG_approxS_norm{ll};
%             [sol_omp_hat(:,i), err_mse_omp_hat, iter_time_omp_hat]=greed_omp_chol(Data(:,i),X_approx_norm,size(X_approx_norm,2),'stopTol',1*Sparsity);
              tic  
              [sol_omp_hat(:,i), err_mse_omp_hat, iter_time_omp_hat]=greed_omp_chol(Data(:,i),faustS_mult{ll},size(X_approx_norm,2),'stopTol',1*Sparsity,'P_trans',trans_faustS_mult{ll});
              t1=toc;  
              err_omp_hat(1,i) = norm(X_norm*Gamma(:,i)-X_approx_norm*sol_omp_hat(:,i))/norm(X_norm*Gamma(:,i));
            err_omp_hat(2,i) = isequal(find(Gamma(:,i)),find(sol_omp_hat(:,i)>1e-4));
            idx_omp = find(sol_omp_hat(:,i));
            resDist(ll+1,k,1,i) = min(norm(points(idx(1)) - points(idx_omp(1))),norm(points(idx(1)) - points(idx_omp(2))));
            resDist(ll+1,k,2,i) = min(norm(points(idx(2)) - points(idx_omp(1))),norm(points(idx(2)) - points(idx_omp(2))));
            compute_Times(ll+1,k,i)=t1;
               tic  
              [sol_omp_hat(:,i), err_mse_omp_hat, iter_time_omp_hat]=greed_omp_chol(Data(:,i),matlab_faustS_mult{ll},size(X_approx_norm,2),'stopTol',1*Sparsity,'P_trans',matlab_trans_faustS_mult{ll});
              t2=toc;
              err_omp_hat(1,i) = norm(X_norm*Gamma(:,i)-X_approx_norm*sol_omp_hat(:,i))/norm(X_norm*Gamma(:,i));
            err_omp_hat(2,i) = isequal(find(Gamma(:,i)),find(sol_omp_hat(:,i)>1e-4));
            idx_omp = find(sol_omp_hat(:,i));
            resDist_matlab(ll+1,k,1,i) = min(norm(points(idx(1)) - points(idx_omp(1))),norm(points(idx(1)) - points(idx_omp(2))));
            resDist_matlab(ll+1,k,2,i) = min(norm(points(idx(2)) - points(idx_omp(1))),norm(points(idx(2)) - points(idx_omp(2))));
            compute_Times_matlab(ll+1,k,i)=t2;
%             % faust 
%             [sol_omp_hat(:,i), err_mse_omp_hat, iter_time_omp_hat]=greed_omp_chol(Data(:,i),faustS_mult{ll},size(X_approx_norm,2),'stopTol',1*Sparsity,'P_trans',trans_faustS_mult{ll});
        end

        save('results_current','resDist','compute_Times');
    end
end
toc
heure = clock ;

matfile = fullfile(pathname, 'output/results_BSL_user');
save(matfile,'resDist','resDist_matlab','RCG_approxS_MEG','nb_approx_MEG','compute_Times','compute_Times_matlab', 'RCG_approxS_MEG');



