%%%% Runtime comparison %%%%
% Script that compares the runtimes of:
%              - Classical dense matrix/vector multiplication
%              - matrix/vector multiplication using a FA�ST
% for different square matrix of power-of-two sizes.

function [t_fact, t_dense]=demo_compTime(M, N_run, n, matrix_filename, mat_cells, status)


t_dense = NaN(N_run,numel(M));
t_fact = NaN(N_run,numel(M));

rand(max(n))*rand(max(n),1);

for run = 1:N_run
    %disp(['run ' num2str(run)])
    %T_start = tic;
    idx_ok = 1:length(status);
    idx_ok(status~=0)=[];
    %for i=numel(M):-1:1
    for i=idx_ok
        D = mat_cells{1,i};
        facts = mat_cells{2,i};
        x = randn(n(i),1);
        
        %%%%%%%%%%%%%% Dense
        tic;
        y1=D*x;
        t_dense(run,i) = toc;
                
        %%%%%%%%%%%%%% Sparse
        tic;
        y2 = facts{M(i)}*x;
        for jj = M(i)-1:-1:1
            y2 = facts{jj} * y2;
        end
        t_fact(run,i) = toc;        
    end
    %toc(T_start);
end
t_fact=mean(t_fact,1);
t_dense=mean(t_dense,1);

%save(time_filename,'t_fact','t_dense')
