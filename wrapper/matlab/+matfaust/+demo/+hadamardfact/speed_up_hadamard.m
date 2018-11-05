%% Description speed_up_hadamard
%
%  This demo makes some time comparison between (Hadamard matrix)-vector multiplication and
%  (Hadamard factorisation i.e a FAµST)-vector multiplication for different dimension
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
%===============================================================================
%> Description speed_up_hadamard
%>
%>  This demo makes some time comparison between (Hadamard matrix)-vector multiplication and
%>  (Hadamard factorisation i.e a FAµST)-vector multiplication for different dimension
%>  of the Hadamard matrix.
%>
%===============================================================================
function speed_up_hadamard()
	import matfaust.Faust

	nb_mult=500;
	Ms=6:12;
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
	dense_trans_times=zeros(nb_mult,nb_dim);
	faust_trans_times=zeros(nb_mult,nb_dim);
	faust_mtimes_trans_times=zeros(nb_mult,nb_dim);


	h = waitbar(0,'speed up hadamard : multiplication time comparison ...');
	for i=1:nb_mult
		waitbar(i/nb_mult);
		for k=1:nb_dim
			n=ns(k);
			hadamard_dense=Hadamard_matrices{k};
			hadamard_faust=Faust(Hadamard_facts{k});

			x=rand(n,1);
			ydense=zeros(n,1);
			yfaust=zeros(n,1);
			ydense_trans=zeros(n,1);
			yfaust_trans=zeros(n,1);
			yfaust_mtimes_trans=zeros(n,1);

			%% multiplication by the hadamard matrix
			tic
			ydense=hadamard_dense*x;
			t1=toc;

			tic
			yfaust=hadamard_faust*x;
			t2=toc;

			if(norm(ydense-yfaust)>threshold)
				error('speed_up hadamard : multiplication problem');
			end

			%% multiplication by the transpose of the hadamard matrix	
			tic
			ydense_trans=hadamard_dense'*x;
			t3=toc;

			tic
			yfaust_trans=hadamard_faust'*x;
			t4=toc;

			tic
			yfaust_mtimes_trans=mtimes_trans(hadamard_faust,x,1);
			t5=toc;
			tic

			if (norm(ydense_trans - ydense_trans) > threshold)
				error('speed_up hadamard : transpose multiplication problem');
			end

			if (yfaust_trans ~= yfaust_mtimes_trans)
				error('speed_up hadamard : transpose multiplication problem');
			end

			dense_times(i,k)=t1;
			faust_times(i,k)=t2;
			dense_trans_times(i,k)=t3;
			faust_trans_times(i,k)=t4;
			faust_mtimes_trans_times(i,k)=t5;
		end
	end
	close(h);

	mean_dense_t = mean(dense_times);
	mean_faust_t = mean(faust_times);
	speed_up = mean_dense_t ./ mean_faust_t;
	mean_dense_trans_t = mean(dense_trans_times);
	mean_faust_trans_t = mean(faust_trans_times);
	mean_faust_mtimes_trans_t = mean(faust_mtimes_trans_times);
	speed_up_trans = mean_dense_trans_t ./ mean_faust_trans_t;
	speed_up_mtimes_trans = mean_dense_trans_t ./ mean_faust_mtimes_trans_t;

	%% Plot the result
	plot_tickness=2.0;
	legend_location = 'NorthWest';
	f=figure;
	subplot(2,2,1);
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
	title('runtime Hadamard A*x');
	set(gca,'XTick',Ms);


	subplot(2,2,2);
	semilogy(Ms,speed_up,'linewidth',plot_tickness);
	hold on
	semilogy(Ms,ones(1,nb_dim),'k','linewidth',plot_tickness);
	grid on
	axis([Ms(1) Ms(end)  min([speed_up,1]) max([speed_up,1])]);
	title('speed-up Hadamard A*x');
	xlabel('log(dim)');
	ylabel('speedup');
	legend('faust','neutral','Location',legend_location);
	set(gca,'XTick',Ms);


	subplot(2,2,3);
	semilogy(Ms,mean_faust_trans_t,'linewidth',plot_tickness);
	hold on
	semilogy(Ms,mean_dense_trans_t,'r','linewidth',plot_tickness);
	semilogy(Ms,mean_faust_mtimes_trans_t,'g','linewidth',plot_tickness);
	ymin=min([mean_dense_trans_t,mean_faust_trans_t,mean_faust_mtimes_trans_t]);
	ymax=max([mean_dense_trans_t,mean_faust_trans_t,mean_faust_mtimes_trans_t]);
	grid on
	axis([Ms(1) Ms(end)  ymin ymax]);
	legend('faust','dense','mtimes_trans','Location',legend_location);
	ylabel('Computed Time (sec)');
	xlabel('log(dim)');
	title('runtime Hadamard A''*x');
	set(gca,'XTick',Ms);


	subplot(2,2,4);
	semilogy(Ms,speed_up_trans,'linewidth',plot_tickness);
	hold on
	semilogy(Ms,ones(1,nb_dim),'k','linewidth',plot_tickness);
	semilogy(Ms,speed_up_mtimes_trans,'g','linewidth',plot_tickness);
	grid on
	axis([Ms(1) Ms(end)  min([speed_up_trans,1,speed_up_mtimes_trans]) max([speed_up_trans,1,speed_up_mtimes_trans])]);
	title('speed-up Hadamard A''*x');
	xlabel('log(dim)');
	ylabel('speedup');
	legend('faust','neutral','mtimes_trans','Location',legend_location);
	set(gca,'XTick',Ms);






	f.Name =['Hadamard Faust-vector multiplication'];
	%% save the figure
	runPath=which(mfilename);
	pathname = fileparts(runPath);
	figure_dir = [ '.' filesep 'Figures'];
	format_fig='-dpng';
	figure_name=[figure_dir filesep 'Hadamard-speed_up'];
	print(figure_name, format_fig);






end
