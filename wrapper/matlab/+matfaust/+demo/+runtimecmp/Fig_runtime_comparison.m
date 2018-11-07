%% Description Fig_runtime_comparison.m
%  Runtime comparison
%
%  This script displays the result of the runtime comparison between
%  faust multiplication and dense matrix multiplication for various
%  configuration of faust (dimension of the faust, number of factors,
%  Relative Complexity Gain  (RCG), fix type of sparsity (sp, spcol,splin))
%
%
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://faust.inria.fr>
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


%% this function displays the resulting C++ computed times between the 2 different
% multiplication (faust vs dense matrix) for different RCG, dimension and
% matrix
%% WARNING this script must be run after @FAUST_SRC_TEST_DIR@/gen_artificial_faust.m
% which store some Faust::Transform example with different dimensions RCG and
% number of factors in
% after run the executable
% @FAUST_BIN_TEST_BIN_DIR@/multiply_compare_time that will do some time
% comparison between dense mat-vector product and faust-vector product
%
% example of inputfile and outputdir
%inputfile='@FAUST_BIN_TEST_OUTPUT_DIR@/multiply_compare_time.mat';
%outputdir='@FAUST_BIN_TEST_FIG_DIR@';
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
