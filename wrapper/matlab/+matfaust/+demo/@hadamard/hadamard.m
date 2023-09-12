%% Description demo_fact_hadamard
%
%  This demo hierarchically factorizes the Hadamard dictionary and then
%  plots the results. This essentially reproduces figure 2 from [1].
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2019):	Hakim Hadj-djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais, Luc Le Magoarou, Remi Gribonval
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
%	Hakim Hadj-Dji. : hakim.hadj-djilani@inria.fr
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
%%-------------------------------------------------------------------------
%% Description norm_hadamard
%
%  This demo makes some time comparison between the 2-norm of the Hadamard matrix and her Faust representation for different dimension of the Hadamard matrix.
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
%%-------------------------------------------------------------------------------
%% Description speed_up_hadamard
%
%  This demo makes some time comparison between (Hadamard matrix)-vector multiplication and
%  (Hadamard factorisation i.e a FAuST)-vector multiplication for different dimension
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

% ======================================================================
%> The demo for the hierarchical factorization of Hadamard matrices.
% ======================================================================
classdef hadamard
	methods(Static)
		%======================================================================
		%> This demo hierarchically factorizes the Hadamard dictionary and then plots the results.
		%===
		%>
		%>  This essentially reproduces figure 2 from [1].
		%>
		%>	<a href="https://hal.archives-ouvertes.fr/hal-01167948v1">[1]</a>
		%> 	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
		%>	approximations of matrices and applications", Journal of Selected
		%>	Topics in Signal Processing, 2016.
		%======================================================================
		function demo_fact_hadamard()



			n = 32;
			M = log2(n);

			%% Generating the data
			[matrix,~] = hadamard_mat(M);





			%% Setting of the parameters
			params_hadamard.nrow=n;
			params_hadamard.ncol=n;
			params_hadamard.nfacts = M;
			params_hadamard.cons = cell(2,M-1);

			for j=1:M-1
				params_hadamard.cons{1,j} = {'splincol',2,n,n, true, false}; % bools: normalized, positive
				params_hadamard.cons{2,j} = {'splincol',n/2^j,n,n, true, false};
			end

			params_hadamard.niter1 = 30;
			params_hadamard.niter2 = 30;
			params_hadamard.update_way = 1;
			params_hadamard.verbose = 0;


			hadamard_faust = faust_decompose(matrix,params_hadamard);
			Xhat = full(hadamard_faust);
			relative_error = norm(matrix - Xhat)/norm(matrix);
			fprintf(['\n\n relative error between hadamard matrix and its transform : ' num2str(relative_error) '\n']);

			nb_mult= 500;
			dense_times=zeros(nb_mult,1);
			faust_times=zeros(nb_mult,1);

			for i=1:nb_mult
				x=rand(n,1);
				y_faust=rand(n,1);
				y_X=rand(n,1);

				tic 
				y_X=matrix*x;
				t1=toc;

				tic 
				y_faust=hadamard_faust*x;
				t2=toc;

				dense_times(i)=t1;
				faust_times(i)=t2;


			end


			dense_t=mean(dense_times);
			faust_t=mean(faust_times);
			speed_up=dense_t/faust_t;

			disp(['multiplication by hadamard matrix : ' num2str(dense_t) 'sec']);
			disp(['multiplication by faust : ' num2str(faust_t) 'sec']);
			disp(['multiplication speed-up using faust : ' num2str(speed_up) ]);


			%% Plotting the results

			% figure configuration
			runPath=which(mfilename);
			pathname = fileparts(runPath);
			figure_dir = [ '.' filesep 'Figures'];
			if(~ exist(figure_dir))
				mkdir(figure_dir)
			end
			format_fig='-dpng';


			fighandle=figure;
			hold on;
			subplot(1,params_hadamard.nfacts+1,1);
			imagesc(Xhat); axis square
			set(gca,'xtick',[],'ytick',[])

			%get the factor of the Faust
			facts=cell(1,M);
			for i=1:M
				facts{i}=factors(hadamard_faust,i);
			end

			for kk = 1:params_hadamard.nfacts
				subplot(1,params_hadamard.nfacts+1,1+kk)
				imagesc(facts{kk}); axis square
				set(gca,'xtick',[],'ytick',[])
			end

			fighandle.Name =['Hadamard-factorisation_image'];
			figure_name = [figure_dir filesep 'Hadamard-factorisation_image'];
			print(figure_name, format_fig);


			fighandle=figure;
			hold on;
			subplot(1,params_hadamard.nfacts+1,1);
			imagesc(Xhat); axis square
			set(gca,'xtick',[],'ytick',[])


			for kk = 1:params_hadamard.nfacts
				subplot(1,params_hadamard.nfacts+1,1+kk)
				spy(facts{kk}); axis square
				set(gca,'xtick',[],'ytick',[])
			end
			fighandle.Name =['Hadamard-factorisation_nnz_coeff'];
			figure_name = [figure_dir filesep 'Hadamard-factorisation_nnz_coeff'];
			print(figure_name, format_fig);





		end

		%======================================================================
		%> This demo makes some time comparison between the 2-norm of the Hadamard matrix and its Faust representation for different dimension of the Hadamard matrix.
		%===
		%>
		%======================================================================
		function norm_hadamard()

			import matfaust.Faust
			nb_mult=10;
			Ms=6:11;
			ns=2.^Ms;
			nb_dim=length(Ms);
			threshold=10^(-10);


			if usejava('jvm')
				h = waitbar(0,'speed up hadamard : Generation of the data ...');
			end
			Hadamard_matrices=cell(1,nb_dim);
			Hadamard_facts=cell(1,nb_dim);





			for k=1:nb_dim
				if usejava('jvm')
					waitbar(k/nb_dim);
				end
				M=Ms(k);
				n=ns(k);
				% generation of the hadamard factorisation
				[H,facts] = hadamard_mat(M);
				Hadamard_matrices{k}=H;
				Hadamard_facts{k}=facts;
			end
			if usejava('jvm')
				close(h);
			end




			dense_times=zeros(nb_mult,nb_dim);
			faust_times=zeros(nb_mult,nb_dim);
			norm_dense=zeros(1,nb_dim);
			norm_faust=zeros(1,nb_dim);
			%RCGs=ns./(Ms*2);
			RCGs=zeros(1,nb_dim);

			if usejava('jvm')
				h = waitbar(0,'2-norm hadamard : multiplication time comparison ...');
			end
			for i=1:nb_mult
				if usejava('jvm')
					waitbar(i/nb_mult);
				end
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
			if usejava('jvm')
				close(h);
			end

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
			if(~ exist(figure_dir))
				mkdir(figure_dir)
			end
			format_fig='-dpng';
			figure_name=[figure_dir filesep 'Hadamard-norm'];
			print(figure_name, format_fig);







		end

		%===============================================================================
		%> This demo makes some time comparison between (Hadamard matrix)-vector multiplication and (Hadamard factorisation i.e a FAuST)-vector multiplication for different dimension of the Hadamard matrix.
		%===============================================================================
		function speed_up_hadamard()
			import matfaust.Faust

			nb_mult=500;
			Ms=6:12;
			ns=2.^Ms;
			nb_dim=length(Ms);
			threshold=10^(-10);


			if usejava('jvm')
				h = waitbar(0,'speed up hadamard : Generation of the data ...');
			end
			Hadamard_matrices=cell(1,nb_dim);
			Hadamard_facts=cell(1,nb_dim);





			for k=1:nb_dim
				if usejava('jvm')
					waitbar(k/nb_dim);
				end
				M=Ms(k);
				n=ns(k);
				% generation of the hadamard factorisation
				[H,facts] = hadamard_mat(M);
				Hadamard_matrices{k}=H;
				Hadamard_facts{k}=facts;
			end
			if usejava('jvm')
				close(h);
			end




			dense_times=zeros(nb_mult,nb_dim);
			faust_times=zeros(nb_mult,nb_dim);
			dense_trans_times=zeros(nb_mult,nb_dim);
			faust_trans_times=zeros(nb_mult,nb_dim);
			faust_mtimes_trans_times=zeros(nb_mult,nb_dim);


			if usejava('jvm')
				h = waitbar(0,'speed up hadamard : multiplication time comparison ...');
			end
			for i=1:nb_mult
				if usejava('jvm')
					waitbar(i/nb_mult);
				end
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
			if usejava('jvm')
				close(h);
			end

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
			if(~ exist(figure_dir))
				mkdir(figure_dir)
			end
			format_fig='-dpng';
			figure_name=[figure_dir filesep 'Hadamard-speed_up'];
			print(figure_name, format_fig);






		end
	end

end
