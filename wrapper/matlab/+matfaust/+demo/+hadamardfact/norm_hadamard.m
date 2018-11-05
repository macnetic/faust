%% Description norm_hadamard
%
%  This demo makes some time comparison between the 2-norm of the Hadamard matrix and
%  her Faust representation for different dimension
%  of the Hadamard matrix.
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
%   Nicolas Bellot	: nicolas.bellot@inria.fr
%   Adrien Leman	: adrien.leman@inria.fr
%   Thomas Gautrais : thomas.gautrais@inria.fr
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%
%======================================================================
%> This demo makes some time comparison between the 2-norm of the Hadamard matrix and
%>  her Faust representation for different dimension
%>  of the Hadamard matrix.
%>
%======================================================================
function norm_hadamard()
	import matfaust.Faust
	nb_mult=10;
	Ms=6:11;
	ns=2.^Ms;
	nb_dim=length(Ms);
	threshold=10^(-10);


	h = waitbar(0,'speed up hadamard : Generation of the data ...');
	Hadamard_matrices=cell(1,nb_dim);
	Hadamard_facts=cell(1,nb_dim);





	for k=1:nb_dim
		waitbar(k/nb_dim);
		M=Ms(k);
		n=ns(k);
		% generation of the hadamard factorisation
		[H,facts] = hadamard_mat(M);
		Hadamard_matrices{k}=H;
		Hadamard_facts{k}=facts;
	end
	close(h);




	dense_times=zeros(nb_mult,nb_dim);
	faust_times=zeros(nb_mult,nb_dim);
	norm_dense=zeros(1,nb_dim);
	norm_faust=zeros(1,nb_dim);
	%RCGs=ns./(Ms*2);
	RCGs=zeros(1,nb_dim);

	h = waitbar(0,'2-norm hadamard : multiplication time comparison ...');
	for i=1:nb_mult
		waitbar(i/nb_mult);
		for k=1:nb_dim
			n=ns(k);
			hadamard_dense=Hadamard_matrices{k};
			hadamard_faust=Faust(Hadamard_facts{k});
			RCGs(k)=rcg(hadamard_faust); 


			%% 2-norm of the hadamard matrix
			tic
			norm_dense(k)=norm(hadamard_dense);
			t1=toc;

			tic
			norm_faust(k)=norm(hadamard_faust);
			t2=toc;





			dense_times(i,k)=t1;
			faust_times(i,k)=t2;

		end
	end
	close(h);

	mean_dense_t = mean(dense_times);
	mean_faust_t = mean(faust_times);
	speed_up = mean_dense_t ./ mean_faust_t;
	% expected_norm = ones(1,nb_dim); % version orthonormee
	expected_norm = 2.^(Ms/2); % version non orthonormée 
	err_norm_dense=sqrt((norm_dense - expected_norm).^2);
	err_norm_faust=sqrt((norm_faust - expected_norm).^2);

	disp(['speed-up : ' num2str(speed_up)]);
	disp(['norm dense : ' num2str(norm_dense(1,:))]);
	disp(['norm faust : ' num2str(norm_faust(1,:))]);


	%% Plot the result
	plot_tickness=2.0;
	legend_location = 'NorthWest';
	f=figure;

	% runtime comparison
	subplot(1,3,1);
	semilogy(Ms,mean_faust_t,'linewidth',plot_tickness);
	hold on
	semilogy(Ms,mean_dense_t,'r','linewidth',plot_tickness);
	ymin=min([mean_dense_t,mean_faust_t]);
	ymax=max([mean_dense_t,mean_faust_t]);
	grid on
	axis([Ms(1) Ms(end)  ymin ymax]);
	legend('faust','dense','Location',legend_location);
	ylabel('Computed Time (sec)');
	xlabel('log(dim)');
	title('runtime ');
	set(gca,'XTick',Ms);

	% speed-up
	subplot(1,3,2);
	semilogy(Ms,speed_up,'linewidth',plot_tickness);
	hold on
	semilogy(Ms,ones(1,nb_dim),'k','linewidth',plot_tickness);
	semilogy(Ms,RCGs,'g','linewidth',plot_tickness);

	grid on
	axis([Ms(1) Ms(end)  min([speed_up,1,RCGs]) max([speed_up,1,RCGs])]);
	title('speed-up norm(A)');
	xlabel('log(dim)');
	ylabel('speedup');
	legend('faust','neutral','theoretical','Location',legend_location);
	set(gca,'XTick',Ms);
	%

	subplot(1,3,3);
	id=find(err_norm_faust ~= 0);% semilogy is not compatible with 0
	semilogy(Ms(id),err_norm_faust(id),'r+-','linewidth',plot_tickness);
	hold on
	semilogy(Ms,err_norm_dense,'linewidth',plot_tickness);

	ymin=min([err_norm_dense, err_norm_faust]);
	ymax=max([err_norm_dense, err_norm_faust]);

	grid on
	%axis([Ms(1) Ms(end)  ymin ymax]);
	legend('faust','dense','Location',legend_location);
	ylabel('Computed Time (sec)');
	xlabel('log(dim)');
	title('error ');
	set(gca,'XTick',Ms);









	f.Name =['Hadamard 2-norm'];
	%% save the figure
	runPath=which(mfilename);
	pathname = fileparts(runPath);
	figure_dir = ['.' filesep 'Figures'];
	format_fig='-dpng';
	figure_name=[figure_dir filesep 'Hadamard-norm'];
	print(figure_name, format_fig);







end
