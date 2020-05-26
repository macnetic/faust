%% License:
% Copyright (2019):	Luc Le Magoarou, Remi Gribonval
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
%   Hakim Hadj-dji.: hakim.hadj-djilani@inria.fr
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
%> @package matfaust.demo @brief The matfaust demo namespace.
% ============================================================================
%> Script used to run all demo (brain source localisation, hadamard factorization, …)
%===
%> @fn matfaust.demo.runall
%> The scripts creates corresponding figures of the article <a href="https://hal.archives-ouvertes.fr/hal-01167948v1">[1]</a>.
%>
%> References:
%>
%> [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
%>	approximations of matrices and applications", Journal of Selected
%>	Topics in Signal Processing, 2016.
%>	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%>
%> @Example
%> @code
%> import matfaust.demo *
%> runall()
%> % You'll find the figures in Figures folder located in current directory.
%> % The benchmark files from which figures are based are in the folder output (in the current directory also).
%> @endcode
%>
% ============================================================================
function runall()
	%% Quick start
	disp('*********** Quick Start Demos *************');
	import matfaust.demo.quickstart.*
	quick_start;
	factorize_matrix;
	construct_Faust_from_factors;





	%% brain source localization
	disp('*********** Brain Source Localization *************');
	import matfaust.demo.bsl.*
	BSL;
	Fig_BSL;

	%% Hadamard factorization
	disp('*********** Hadamard Factorization *************');
	import matfaust.demo.hadamard.*
	demo_fact_hadamard;
	speed_up_hadamard;
	norm_hadamard;


	%% Fourier speed-up
	disp('*********** Fourier speed-up *************');
	import matfaust.demo.fft.*
	speed_up_fourier;



	%% Runtime comparison
	disp('*********** Runtime Comparison *************');
	import matfaust.demo.runtimecmp.*
	runtime_comparison;
	Fig_runtime_comparison;
end
