close all
clear
clc
N_run = 100;
M = 1:13;
%M = 1:5;
n = 2.^M;
t_dense = zeros(N_run,numel(M));
t_fact = zeros(N_run,numel(M));
mat_cells = cell(2,numel(M));
density_max = 0.3;

for i=1:numel(M)
    disp(i)
    density = 0;
    while density<=0.5
        facts = cell(1,M(i));
        for j = 1:M(i)
            facts{j} = [];
            facts_i_cell = cell(1,n(i)/2);
            for k=1:n(i)/2
                facts_i_cell{k} = randn(2);
                %facts{j} = blkdiag(facts{j},randn(2));
            end
            facts{j} = blkdiag(facts_i_cell{:});
            facts{j} = facts{j}(randperm(n(i)),:);
            facts{j} = facts{j}(:,randperm(n(i)));% figure; imagesc(facts{j})
            facts{j} = sparse(facts{j});
        end
        mat_cells{1,i} = full(dvp(facts));% figure; imagesc(D)
        mat_cells{2,i} = facts;
        density = nnz(mat_cells{1,i})/n(i)^2;
        disp(['density: ' num2str(density)])
    end
end
%%save(['Matrix_' int2str(n(end)) '.mat'],'mat_cells','-v7.3');
save(['Matrix_' int2str(n(end)) '.mat'],'mat_cells','-v7.3');
%%save('Matrix_32768','mat_cells','-v7.3')%save('Matrices','mat_cells')