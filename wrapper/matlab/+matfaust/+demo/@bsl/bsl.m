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
%	Hakim Hadj-dji. : hakim.hadj-djilani@inria.fr
%	Nicolas Bellot	: nicolas.bellot@inria.fr
%   Adrien Leman	: adrien.leman@inria.fr
%   Thomas Gautrais	: thomas.gautrais@inria.fr
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr

% ======================================================================
%> @brief Brain Source Localization demo.
% ======================================================================
classdef bsl
	properties(SetAccess = public)
	end
	methods(Static)

		%===============================================================================
		%>  Brain source localization
		%===
		%>  This script performs brain source localization using several gain
		%>  matrices <a href="http://www.ncbi.nlm.nih.gov/pubmed/24161808">[2]</a>, including FAuSTs, and OMP solver. It reproduces the
		%>  source localization experiment of <a href="https://hal.archives-ouvertes.fr/hal-01167948v1">[1]</a>.
		%>
		%>  The results are stored in "./output/results_BSL_user.mat".
		%>  DURATION: Computations should take around 3 minutes.
		%>
		%>  The MEG gain matrices used are
		%>         the precomputed ones in "data/faust_MEG_rcg_X.mat" (in the installation directory of the FAuST toolbox)
		%>
		%> References:
		%>
		%> [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
		%>	approximations of matrices and applications", Journal of Selected
		%>	Topics in Signal Processing, 2016.
		%>	<https://hal.archives-ouvertes.fr/hal-01167948v1>
		%>
		%> [2]   A. Gramfort, M. Luessi, E. Larson, D. Engemann, D. Strohmeier,
		%>	C. Brodbeck, L. Parkkonen, M. Hamalainen, MNE software for processing
		%>	MEG and EEG data <http://www.ncbi.nlm.nih.gov/pubmed/24161808>,
		%>	NeuroImage, Volume 86, 1 February 2014, Pages 446-460, ISSN 1053-8119,
		%>
		%===============================================================================
		function BSL()

			import matfaust.Faust

			runPath=which(mfilename);
			pathname = fileparts(runPath);

			% folder_path where the MEG and Faust approximations are stored
			BSL_data_pathName=[pathname '..' filesep '..' filesep '..' filesep 'data'];
			% useless if the faust toolbox has been properly added in the matlab path


			% files where a precomputed FauST approximation of MEG matrix are stored
			MEG_FAuST_list_filename={'faust_MEG_rcg_6.mat','faust_MEG_rcg_8.mat','faust_MEG_rcg_16.mat','faust_MEG_rcg_25.mat'};

			nb_FAuST_MEG = length(MEG_FAuST_list_filename);

			% list of the MEG matrices (the original one and her Faust approximations)
			MEG_list = cell(1,nb_FAuST_MEG+1);
			nb_MEG_matrix = length(MEG_list);

			%% Loading the MEG matrix and the points in the brain
			%load([BSL_data_pathName filesep 'matrix_MEG.mat' ]);
			load('matrix_MEG.mat')
			MEG_matrix = normalizeCol(matrix');% normalization of the columns of the MEG matrix
			MEG_list{1}=MEG_matrix;





			%% Loading the precomputed FAuST representing the MEG matrix
			% RCG (Relative Complexity Gain) : the theoretical speed-up using Faust
			RCG_list=zeros(1,nb_FAuST_MEG);
			for i=1:nb_FAuST_MEG

				%load([BSL_data_pathName filesep  MEG_FAuST_list_filename{i}]);
				load(MEG_FAuST_list_filename{i});
				facts = normalizeCol(facts,lambda); % normalization of the columns of the FAUST
				MEG_FAuST=Faust(facts); % construct the FAuST from its factorscons15_row
				MEG_list{i+1}=MEG_FAuST; % store the different FAuST approximationscons21_col
				RCG_list(i)=rcg(MEG_FAuST); % compute the RCG of the given FAuST
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

			if(~ exist('output'))
				mkdir 'output'
			end
			matfile = [ 'output' filesep 'results_BSL_user' ];
			save(matfile,'resDist','Sparsity','RCG_list','compute_Times','Ntraining','nb_MEG_matrix');

		end

		%===============================================================================
		%>  Brain source localization figures.
		%===
		%>  This script builds the BSL figure (Fig. 9.) used in <a href="https://hal.archives-ouvertes.fr/hal-01167948v1">[1]</a>.
		%>
		%> References:
		%>
		%> [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
		%>	approximations of matrices and applications", Journal of Selected
		%>	Topics in Signal Processing, 2016.
		%>	<https://hal.archives-ouvertes.fr/hal-01167948v1>
		%>
		%===============================================================================
		function Fig_BSL()

			runPath=which(mfilename);
			pathname = fileparts(runPath);
			matfile = ['output' filesep 'results_BSL_user.mat'];
			if (not(exist(matfile)))
				error('run BSL.m before Fig_BSL.m');
			end

			%% figure configuration
			figure_dir = ['.' filesep 'Figures'];
			if(~ exist(figure_dir))
				mkdir(figure_dir)
			end
			format_fig='-dpng';

			load(matfile);
			nb_approx_MEG=nb_MEG_matrix-1;
			solver_choice='omp';
			Ntest=Ntraining*Sparsity;



			%% convergence analysis
			d1 = cat(4,resDist(:,1,1,:),resDist(:,1,2,:));
			d2 = cat(4,resDist(:,2,1,:),resDist(:,2,2,:));
			d3 = cat(4,resDist(:,3,1,:),resDist(:,3,2,:));
			test2 = 100*[squeeze(d1);zeros(1,Ntest);squeeze(d2);zeros(1,Ntest);squeeze(d3)];


			f=figure('color',[1 1 1]);
			T = bplot(test2','linewidth',1.5);
			legend(T);

			ylabel('Distance between true and estimated sources (cm)')
			box on
			ax = gca;
			x_tick=[1:nb_approx_MEG+1;nb_approx_MEG+3:2*nb_approx_MEG+3;2*nb_approx_MEG+5:3*nb_approx_MEG+5];
			mean_x_tick = round(mean(x_tick,2));
			set(ax,'XTick',reshape(x_tick',1,numel(x_tick)));
			set(ax,'xticklabel', [])
			axi = axis;

			axis([0 3*(nb_approx_MEG+2) -0.3 4.7])
			HorizontalOffset = 10;
			verticalOffset = 0.43;
			verticalOffset2 = 0.7;
			yTicks = get(ax,'ytick');
			xTicks = get(ax, 'xtick');
			minY = min(yTicks);

			text(xTicks(1), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
			text(xTicks(2+nb_approx_MEG), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
			text(xTicks(1+2*(nb_approx_MEG+1)), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
			text(xTicks(round(mean(1:(nb_approx_MEG+1)))), minY - verticalOffset2, '$1<d<5$','HorizontalAlignment','center','interpreter', 'latex');
			text(xTicks(round(mean(nb_approx_MEG+2:2*(nb_approx_MEG+1)))), minY - verticalOffset2, '$5<d<8$','HorizontalAlignment','center','interpreter', 'latex');
			text(xTicks(round(mean(2*(nb_approx_MEG+1)+1:3*(nb_approx_MEG+1)))), minY - verticalOffset2, '$d>8$','HorizontalAlignment','center','interpreter', 'latex');





			for i=1:nb_approx_MEG
				text(xTicks(i+1), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_list(i)) '}$' ],'HorizontalAlignment','center','interpreter', 'latex');
				text(xTicks(i+2+nb_approx_MEG), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_list(i)) ' }$' ],'HorizontalAlignment','center','interpreter', 'latex');
				text(xTicks(i+1+2*(nb_approx_MEG+1)), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_list(i)) ' }$' ],'HorizontalAlignment','center','interpreter', 'latex');

			end

			title(['BSL - convergence (C++ wrapper faust) ' solver_choice ' solver']);
			f.Name =['Brain Source Localization : convergence with ' solver_choice '_solver (C++ wrapper)'];
			figure_name=[figure_dir filesep 'BSL-convergence_Cpp_' solver_choice '_solver'];
			print(figure_name, format_fig);





			%% time comparison

			timeS = cat(3,compute_Times(:,1,:),compute_Times(:,2,:),compute_Times(:,3,:));
			timeS = squeeze(1000*timeS)';



			f=figure('color',[1 1 1]);
			T = bplot(timeS,'linewidth',1.5);
			legend(T);
			ylabel('Computed Time (ms)')
			box on
			ax = gca;
			x_tick=[1:2*nb_approx_MEG+1];
			mean_x_tick = round(mean(x_tick,2));
			set(ax,'XTick',x_tick');
			set(ax,'xticklabel', [])
			axi = axis;


			HorizontalOffset = 10;
			verticalOffset = 0.5;
			verticalOffset2 = 1;
			yTicks = get(ax,'ytick');
			xTicks = get(ax, 'xtick');
			minY = min(yTicks);

			text(xTicks(1), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');


			for i=1:nb_approx_MEG
				text(xTicks(i+1), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_list(i))  '}$' ],'HorizontalAlignment','center','interpreter', 'latex');

			end
			title(['BSL - time comparison (FAUST vs dense matrix) ' solver_choice ' solver']);
			f.Name =['Brain Source Localization : time comparison with ' solver_choice 'solver'];
			figure_name=[figure_dir filesep 'BSL-time_comparison_' solver_choice '_solver'];
			print(figure_name, format_fig);







			%% speed-up
			mean_Times=mean(timeS);
			dense_matrix_time=mean_Times(1);
			real_RCG=dense_matrix_time./mean_Times;

			f=figure;
			ax = gca;

			set(ax,'xticklabel', [])
			set(ax,'Visible','off');
			plot(1:nb_approx_MEG,real_RCG(2:end),'linewidth',1.5);
			hold on
			plot(1:nb_approx_MEG,ones(1,nb_approx_MEG),'linewidth',1.5);
			verticalOffset=1.5;
			minY=min([0.9,real_RCG(2:end)]);
			maxY=max([0.9,real_RCG(2:end)]);
			axis([1,nb_approx_MEG,minY,maxY]);
			set(ax,'XTick',[]);
			for i=1:nb_approx_MEG
				text(i, minY -(maxY-minY)/20, ['$\widehat{\mathbf{M}}_{' int2str(RCG_list(i)) '}$' ],'HorizontalAlignment','center','interpreter', 'latex');

			end
			legend('speed up FAuST','neutral speed up');
			title(['BSL - speed up using FAUST ' solver_choice ' solver']);
			f.Name =['Brain Source Localization : speed-up Faust with ' solver_choice 'solver'];
			figure_name = [figure_dir filesep 'BSL-speed_up_' solver_choice ' solver'];
			print(figure_name, format_fig);



			%% console display
			disp(['**** MEG with OMP solver time comparison ****']);
			disp(['M tps : ' num2str(mean(timeS(:,1))) ]);
			for i=1:nb_approx_MEG
				disp([ 'M_' int2str(RCG_list(i)) ' tps : ' num2str(mean(timeS(:,i+1))) ', speed-up : ' num2str(real_RCG(i+1))]);
			end
		end
	end
end
