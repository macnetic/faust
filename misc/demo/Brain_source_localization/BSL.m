%% Description BSL
%  Brain source localization
%
%  This script performs brain source localization using several gain
%  matrices [2], including FAuSTs, and OMP solver. It reproduces the
%  source localization experiment of [1].
%  The results are stored in "./output/results_BSL_user.mat".
%  DURATION: Computations should take around 3 minutes.
%
%  The MEG gain matrices used are
%         the precomputed ones in "./data/M_X.mat"  
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot, Adrien Leman, Thomas Gautrais, Luc Le Magoarou, Remi Gribonval
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
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr
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






runPath=which(mfilename);
pathname = fileparts(runPath);

% folder_path where the MEG and Faust approximations are stored
BSL_data_pathName=[pathname filesep 'data'];



% files where a precomputed FauST approximation of MEG matrix are stored
MEG_FAuST_list_filename={'M_6.mat','M_8.mat','M_16.mat','M_25.mat'};

nb_FAuST_MEG = length(MEG_FAuST_list_filename);

% list of the MEG matrices (the original one and her Faust approximations)
MEG_list = cell(1,nb_FAuST_MEG+1);
nb_MEG_matrix = length(MEG_list);

%% Loading the MEG matrix and the points in the brain
load([BSL_data_pathName filesep 'X_meg.mat' ]);
MEG_matrix = normalizeCol(X_fixed);% normalization of the columns of the MEG matrix
MEG_list{1}=MEG_matrix;





%% Loading the precomputed FAuST representing the MEG matrix
% RCG (Relative Complexity Gain) : the theoretical speed-up using Faust
RCG_list=zeros(1,nb_FAuST_MEG);
for i=1:nb_FAuST_MEG
    
    load([BSL_data_pathName filesep  MEG_FAuST_list_filename{i}]);
    facts = normalizeCol(facts,lambda); % normalization of the columns of the FAUST
    MEG_FAuST=Faust(facts); % construct the FAuST from its factorscons15_row
    MEG_list{i+1}=MEG_FAuST; % store the different FAuST approximationscons21_col
    RCG_list(i)=RCG(MEG_FAuST); % compute the RCG of the given FAuST
end


M=size(MEG_matrix,2); % Number of points used in the MEG matrix
Ntraining = 500; % Number of training vectorscons22_row
Sparsity = 2; % Number of sources per training vector
dist_paliers = [0.01,0.05,0.08,0.5];



resDist = zeros(nb_MEG_matrix,numel(dist_paliers)-1,Sparsity,Ntraining);% (Matrice,dist_sources,src_nb,run);
compute_Times = zeros(nb_MEG_matrix,numel(dist_paliers)-1,Ntraining);


h = waitbar(0,['Brain Source Localization : MEG matrix and its faust approximations with omp solver']);
nb_palier=numel(dist_paliers)-1;
for k=1:nb_palier;
    
    %Parameters settings
    % generates the different source position
    Gamma = zeros(M,Ntraining);
    for ii=1:Ntraining
        dist_sources = -1;
        while ~((dist_paliers(k)<dist_sources)&&(dist_sources<dist_paliers(k+1)))
            Gamma(:,ii) = sparse_coeffs(MEG_matrix, 1, Sparsity); %
            idx = find(Gamma(:,ii));
            dist_sources = norm(points(idx(1),:) - points(idx(2),:));
        end
        
    end
    
    % compute the data registered by MEG sensor
    Data = MEG_matrix*Gamma;
    
    for i=1:Ntraining
        waitbar(((k-1)*Ntraining+i)/(nb_palier*Ntraining));
        
        % index of the real source localization
        idx = find(Gamma(:,i));

        
        for j=1:nb_MEG_matrix
            
            MEG_FAuST = MEG_list{j};
            
            % find the active source
            tic
            sol_solver=greed_omp_chol(Data(:,i),MEG_FAuST,M,'stopTol',1*Sparsity,'verbose',false);
            compute_Times(j,k,i)=toc;

            % compute the distance between the estimated source and the real one
            idx_solver = find(sol_solver);
            resDist(j,k,1,i) = min(norm(points(idx(1),:) - points(idx_solver(1),:)),norm(points(idx(1),:) - points(idx_solver(2),:)));
            resDist(j,k,2,i) = min(norm(points(idx(2),:) - points(idx_solver(1),:)),norm(points(idx(2),:) - points(idx_solver(2),:)));
        end
    end
end
close(h);

matfile = fullfile(pathname, 'output/results_BSL_user');
save(matfile,'resDist','Sparsity','RCG_list','compute_Times','Ntraining','nb_MEG_matrix');








