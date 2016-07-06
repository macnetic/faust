%% Description speed_up_hadamard
%
%  This demo makes some time comparison between (Hadamard matrix)-vector multiplication and
%  (Hadamard factorisation i.e a FAµST)-vector multiplication for different dimension
%  of the Hadamard matrix.
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
%%

nb_mult=500;
Ms=6:11;
ns=2.^Ms;
nb_dim=length(Ms);
threshold=10^(-10);
dense_times=zeros(nb_mult,nb_dim);
faust_times=zeros(nb_mult,nb_dim);

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

%     figure,
%         for i=1:M
%             subplot(2,M,i);
%             imagesc(facts{i});
%             axis image
%             subplot(2,M,M+i);
%             imagesc(cum_Hbis{i});
%             axis image
%         end
hadamard_faust=matlab_faust(facts);
hadamard_dense=dvp(facts);

h = waitbar(0,'speed up hadamard : multiplication time comparison ...');
for i=1:nb_mult
    waitbar(i/nb_mult);
    for k=1:nb_dim
        n=ns(k);
        hadamard_dense=Hadamard_matrices{k};
        hadamard_faust=matlab_faust(Hadamard_facts{k});
        
        x=rand(n,1);
        ydense=zeros(n,1);
        yfaust=zeros(n,1);
        
        tic
        ydense=hadamard_dense*x;
        t1=toc;
        
        tic
        yfaust=hadamard_faust*x;
        t2=toc;
        
        if(norm(ydense-yfaust)>threshold)
            error('speed_up hadamard : multiplication problem');
        end
        
        
        dense_times(i,k)=t1;
        faust_times(i,k)=t2;
    end
end
close(h);

mean_dense_t = mean(dense_times);
mean_faust_t = mean(faust_times);
speed_up = mean_dense_t ./ mean_faust_t;

%% Plot the result

f=figure;
subplot(1,2,1);
semilogy(Ms,mean_faust_t,'linewidth',1.5);
hold on
semilogy(Ms,mean_dense_t,'r','linewidth',1.5);
ymin=min([mean_dense_t,mean_faust_t]);
ymax=max([mean_dense_t,mean_faust_t]);
grid on
axis([Ms(1) Ms(end)  ymin ymax]);
legend('faust','dense');
ylabel('Computed Time (sec)');
xlabel('log(dim)');
title('Hadamard-vector multiplication');
set(gca,'XTick',Ms);


subplot(1,2,2);
semilogy(Ms,speed_up,'linewidth',1.5);
hold on
semilogy(Ms,ones(1,nb_dim),'g','linewidth',1.5);
grid on
axis([Ms(1) Ms(end)  min([speed_up,1]) max([speed_up,1])]);
title('Hadamard-vector multiplication');
xlabel('log(dim)');
ylabel('speedup');
legend('speed-up','neutral speed-up');
set(gca,'XTick',Ms);
f.Name =['Hadamard Faust-vector multiplication'];
%% save the figure
runPath=which(mfilename);
pathname = fileparts(runPath);
fig_filename = [pathname filesep 'output' filesep 'speed_up_hadamard.png'];
print(fig_filename,'-dpng','-r300');





