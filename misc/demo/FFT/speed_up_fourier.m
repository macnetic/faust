%% Description speed_up_Fourier
%
%  This demo makes some time comparison between (Fourier matrix)-vector multiplication and
%  (Fourier factorisation i.e a FAµST)-vector multiplication for different dimension
%  of the Fourier matrix.
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot,  Luc Le Magoarou, Remi Gribonval
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
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
%	approximations of matrices and applications", Journal of Selected
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%

nb_mult=100;
Ms=6:12;
ns=2.^Ms;
nb_dim=length(Ms);
threshold=10^(-10);


h = waitbar(0,'speed up Fourier : Generation of the data ...');
Fourier_matrices=cell(1,nb_dim);
Fourier_facts=cell(1,nb_dim);





for k=1:nb_dim
    waitbar(k/nb_dim);
    M=Ms(k);
    n=ns(k);
    % generation of the Fourier factorisation
    [H,facts] = CooleyTukeyFact(n);
    Fourier_matrices{k}=H;
    Fourier_facts{k}=facts;
end
close(h);





nb_type_mult = 3; % dense, faust, fft
mult_times_tab=zeros(nb_mult,nb_dim,nb_type_mult);
dense_trans_times=zeros(nb_mult,nb_dim);
faust_trans_times=zeros(nb_mult,nb_dim);
faust_mtimes_trans_times=zeros(nb_mult,nb_dim);


h = waitbar(0,'speed up Fourier : multiplication time comparison ...');
for i=1:nb_mult
    waitbar(i/nb_mult);
    for k=1:nb_dim
        n=ns(k);
        Fourier_dense=Fourier_matrices{k};
        Fourier_faust=Faust(Fourier_facts{k});
        
        x=rand(n,1);
        ydense=zeros(n,1);
        yfaust=zeros(n,1);
        yfft=zeros(n,1);
        ydense_trans=zeros(n,1);
        yfaust_trans=zeros(n,1);
        yfaust_mtimes_trans=zeros(n,1);
        
        %% multiplication by the Fourier matrix
        tic
        ydense=Fourier_dense*x;
        t1=toc;
        
        tic
        yfaust=Fourier_faust*x;
        t2=toc;
        
        tic
        yfft=fft(x);
        t3=toc;
        
        if(norm(ydense-yfaust)>threshold)
            error('speed_up Fourier : multiplication problem');
        end
        
        if(norm(yfft-yfaust)>threshold)
            error('speed_up Fourier : multiplication problem (fft vs faust)');
        end
        

        
        mult_times_tab(i,k,1)=t1;
        mult_times_tab(i,k,2)=t2;
        mult_times_tab(i,k,3)=t3;

       
    end
end
close(h);

mean_mult_times_tab = squeeze(mean(mult_times_tab));




%% Plot the result
plot_tickness=2.0;
legend_location = 'NorthWest';
f=figure;
listcolor={'ro-','bo-','go-','r+-','b+-','g+-'};
subplot(2,1,1);

for i=1:nb_type_mult
    semilogy(Ms,mean_mult_times_tab(:,i),listcolor{i},'linewidth',plot_tickness);
    
    hold on
end
grid on
axis([Ms(1) Ms(end)  min(mean_mult_times_tab(:)) max(mean_mult_times_tab(:))]);
legend('dense','faust','fft','Location',legend_location);
ylabel('Computed Time (sec)');
title('runtime Fourier A*x');
set(gca,'XTick',Ms);

subplot(2,1,2);
for i=2:nb_type_mult
    semilogy(Ms,mean_mult_times_tab(:,1)./mean_mult_times_tab(:,i),listcolor{i},'linewidth',plot_tickness);
    hold on
end
legend('faust','fft','Location',legend_location);
grid on
title('speed-up Fourier A*x');
xlabel('log(dim)');

f.Name =['Fourier Faust-vector multiplication'];
%% save the figure
runPath=which(mfilename);
pathname = fileparts(runPath);
figure_dir = [pathname filesep '..' filesep 'Figures'];
format_fig='-dpng';
figure_name=[figure_dir filesep 'Fourier-speed_up'];
print(figure_name, format_fig);


