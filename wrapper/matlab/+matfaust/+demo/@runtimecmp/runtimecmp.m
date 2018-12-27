
%% License:
% Copyright (2018):	Hakim Hadj-dji., Nicolas Bellot, Adrien Leman, Thomas Gautrais, Luc Le Magoarou, Remi Gribonval
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


% =====================================================================
%> Runtime comparison demo: Faust-vector and dense matrix-vector multiplications, differents RCGs, transpose.
% =====================================================================
classdef runtimecmp
	methods(Static)
		%============================================================================
		%>   Runtime comparison
		%===
		%>  This script performs runtime comparison between faust multiplication
		%>  and dense matrix multiplication for various configuration of faust
		%>  (dimension of the faust, number of factors, Relative Complexity Gain
		%>  (RCG), fix type of sparsity (sp, spcol,splin))
		%>
		%============================================================================
		function runtime_comparison()
			import matfaust.Faust


			params.RCGs = [2 4 8];
			params.Dims = [128 256 512];
			params.nb_facts = [2,4,8];

			params.constraint = 'sp';
			%params.constraint = 'sp_col';
			%params.constraint = 'sp_row';
			params.matrix_or_vector='vector';
			%matrix_or_vector='matrix';
			params.Nb_mult=500;
			[RCGs,Dims,nb_facts,constraint,Nb_mult,matrix_or_vector]=deal(params.RCGs,params.Dims,params.nb_facts,params.constraint,params.Nb_mult,params.matrix_or_vector);




			NDims=length(Dims);
			NRCGs=length(RCGs);
			Nnb_facts=length(nb_facts);








			list_faust=cell(NDims,NRCGs,Nnb_facts);
			list_dense=cell(NDims,1);

			%% loading all the different faust and dense matrix
			for j=1:NDims
				dim=Dims(j);
				A=rand(dim);
				list_dense{j}=A;



				for k=1:NRCGs
					RCG=RCGs(k);

					for l=1:Nnb_facts
						nfact=nb_facts(l);
						fact=gen_artificial_faust(dim,RCG,nfact,constraint);
						faust_transform=Faust(fact);
						list_faust{j,k,l}=faust_transform;
						taillefaust=size(faust_transform);
						if((taillefaust(1) ~= dim) + (taillefaust(2) ~= dim))
							error('invalid faust');
						end


					end
				end
			end



			%% time comparison

			t_faust=zeros(Nb_mult+1,NDims,NRCGs,Nnb_facts,2);
			t_dense=zeros(Nb_mult+1,NDims,2);

			h = waitbar(0,'runtime comparison (faust vs dense matrix) for various configuration ...');
			for i=1:Nb_mult+1
				waitbar(i/(Nb_mult+1));

				for j=1:NDims
					dim=Dims(j);

					if strcmp(matrix_or_vector,'matrix')
						dim2 = dim; % multiplication by a square matrix
					elseif strcmp(matrix_or_vector,'vector')
						dim2 = 1; % multiplication by a column-vector
					else
						error('matrix_or_vector string must be equal to matrix or vector');
					end


					for k=1:NRCGs
						RCG=RCGs(k);

						for l=1:Nnb_facts
							nfact=nb_facts(l);
							fact=gen_artificial_faust(dim,RCG,nfact,constraint);
							faust_transform=Faust(fact);
							x=rand(dim,dim2);
							y=zeros(dim,dim2);
							yfaust=zeros(dim,dim2);
							y_trans=zeros(dim,dim2);
							yfaust_trans=zeros(dim,dim2);



							taillefaust=size(faust_transform);
							if((taillefaust(1) ~= dim)+(taillefaust(2) ~= dim))
								error('invalid faust');
							end
							%% multiplication dense
							if(k==1)&&(l==1)
								A=list_dense{j};
								tic;
								y=A*x;
								tdense=toc;
								t_dense(i,j,1)=tdense;

								tic
								y_trans=A'*x;
								tdense_trans=toc;
								t_dense(i,j,2)=tdense_trans;
							end
							%% multiplication par un faust
							faust_transform=list_faust{j,k,l};
							tic;
							yfaust=faust_transform*x;
							tfaust=toc;
							t_faust(i,j,k,l,1)=tfaust;
							tic;
							yfaust_trans=mtimes_trans(faust_transform,x,1);
							tfaust_trans=toc;
							t_faust(i,j,k,l,2)=tfaust_trans;
						end
					end
				end
			end
			close(h);

			t_faust(1,:,:,:,:)=[];
			t_dense(1,:,:)=[];

			runPath=which(mfilename);
			pathname = 'output'
			if(~ exist(pathname))
				mkdir(pathname)
			end
			matfile = fullfile(pathname, 'runtime_comparison.mat');
			save(matfile,'t_faust','t_dense','params');

		end
		%=========================================================================
		%> Runtime comparison figure
		%===
		%>
		%>  This script displays the result of the runtime comparison between
		%>  faust multiplication and dense matrix multiplication for various
		%>  configuration of faust (dimension of the faust, number of factors,
		%>  Relative Complexity Gain  (RCG), fix type of sparsity (sp, spcol,splin))
		%>
		%=========================================================================
		function Fig_runtime_comparison()
			runPath=which(mfilename);
			pathname = 'output'
			matfile = fullfile(pathname, 'runtime_comparison.mat');
			if (not(exist(matfile)))
				error('run runtime_comparison.m before Fig_runtime_comparison.m');
			end
			load(matfile);

			% figure configuration
			figure_dir = [ '.' filesep 'Figures'];
			if(~ exist(figure_dir))
				mkdir(figure_dir)
			end
			format_fig='-dpng';



			load(matfile);
			[RCGS,DIMS,NB_FACTS,constraint,Nb_mult,constraint,matrix_or_vector]=deal(params.RCGs,params.Dims,params.nb_facts,params.constraint,params.Nb_mult,params.constraint,params.matrix_or_vector);

			nDIMS=length(DIMS);
			nRCGS=length(RCGS);
			nNB_FACTS=length(NB_FACTS);


			%% average the multiplication according to the number of run

			mean_tdense=squeeze(mean(t_dense));
			mean_tfaust=squeeze(mean(t_faust));









			%% plot the time computed in logscale with a fix number of factor in a given figure
			% in each figure we have in x-axis the log2(dimension) of the square matrix
			%                         in y-axis the time
			% all time for faust multiplication with different RCG (density)
			% and the time for dense matrix multiplication are plotted

			legend_curve=cell(1,nRCGS+1);
			thickness_curve = 2;



			%% hold the different figure are in the same box (scale)
			ymin=min([min(mean_tdense(:)),min(mean_tfaust(:))]);
			ymax=max([max(mean_tdense(:)),max(mean_tfaust(:))]);
			fighandle=figure;

			legendy={'Time (A*x)','Time (A''*x)'};
			for h=1:2
				for nfact=1:nNB_FACTS
					subplot(2,nNB_FACTS,(h-1)*nNB_FACTS+nfact);

					for k=1:nRCGS
						semilogy(log2(DIMS),mean_tfaust(:,k,nfact,h),'LineWidth',thickness_curve);
						legend_curve{k}=['Faust RCG ' num2str(RCGS(k))];
						hold on;

					end

					%LA boucle fait remettre les compteurs a zero, que ce soit des
					%numero de courbes ou des couleurs, c'est pourquoi il faut la
					%mettre en premier et
					%forcer la couleur du dense (voire ci dessous).

					semilogy(log2(DIMS),squeeze(mean_tdense(:,h)),'-+','Color',[0 0.8 0.8],'LineWidth',thickness_curve);
					legend_curve{nRCGS+1}=['Dense '];


					% legend the figure,title, font ...
					grid on;
					axis([log2(DIMS(1)) log2(DIMS(end)) ymin ymax]);

					title (['nb factor : ' ,int2str(NB_FACTS(nfact))] );
					if (nfact == 1)
						legend(legend_curve{:},'Location', 'NorthWest');
						xlabel('log2(Dimension)');
						ylabel(legendy{h});
					end

				end

			end

			fighandle.Name =['Faust-' matrix_or_vector ' multiplication (constraint : ' constraint ')'];
			figure_name = [figure_dir filesep 'RuntimeComp-' matrix_or_vector '_multiplication_constraint_' constraint];
			print(figure_name, format_fig);

		end
	end
end
