function [mat_cells,status]=demo_compTime_matrix_draw(M, n, densite_or_nnzrowcol, matrix_filename)


mat_cells = cell(2,numel(M));
nb_fact = M;
density_min = 0.5;
BOUCLE_MAX = 20;
    
densite = zeros(size(M))-1.0;
nz = zeros(size(M))-1.0;
for k=1:numel(M)
    if densite_or_nnzrowcol(k) <= 1
        % nombre de valeurs non nulles par ligne ou par colonne
        nz(k) = round(densite_or_nnzrowcol(k)*n(k));
        densite(k) = nz(k)/n(k);
        
    else
        % nombre de valeurs non nulles par ligne ou par colonne
        nz(k) = densite_or_nnzrowcol(k);
        densite(k) = nz(k)/n(k);
    end
end
status = -ones(size(M));
for i=1:numel(M)
    density = 0;
    ind_boucle=1;
    while( density<= density_min && ind_boucle<=BOUCLE_MAX)
        facts = cell(1,nb_fact(i));
        for j = 1:nb_fact(i)
            facts{j} = zeros(n(i));
            ind = randperm(n(i), nz(i));
            for l=1:nz(i)
                for k=1:ind(l)
                    facts{j}(ind(l)-k+1,k) = rand();
                end
                for k=1:n(i)-ind(l)
                    facts{j}(n(i)-k+1,ind(l)+k) = rand();
                end
            end
            facts{j} = facts{j}(randperm(n(i)),:);
            facts{j} = facts{j}(:,randperm(n(i)));
            facts{j} = sparse(facts{j});
            %fprintf('-');
        end
        mat_cells{1,i} = full(dvp(facts));% figure; imagesc(D)
        mat_cells{2,i} = facts;
        density = nnz(mat_cells{1,i})/n(i)^2;
        ind_boucle = ind_boucle+1;
    end
    if(ind_boucle>BOUCLE_MAX)
        status(i) = -1;
    else
        status(i) = 0;
    end
    %fprintf('\n');
    %disp(['i=' num2str(i) ' ; M=' num2str(M(i)) ' ; density=' num2str(density)]);
end
save(matrix_filename,'mat_cells','-v7.3');%save('Matrices','mat_cells')




%% version de Luc
% for i=1:numel(M)
%     %disp(i);
%     density = 0;
%     while density<=0.5
%         facts = cell(1,M(i));
%         for j = 1:M(i)
%             facts{j} = [];
%             facts_i_cell = cell(1,n(i)/2);
%             for k=1:n(i)/2
%                 facts_i_cell{k} = randn(2);
%                 %facts{j} = blkdiag(facts{j},randn(2));
%             end
%             facts{j} = blkdiag(facts_i_cell{:});
%             facts{j} = facts{j}(randperm(n(i)),:);
%             facts{j} = facts{j}(:,randperm(n(i)));% figure; imagesc(facts{j})
%             facts{j} = sparse(facts{j});
%         end
%         mat_cells{1,i} = full(dvp(facts));% figure; imagesc(D)
%         mat_cells{2,i} = facts;
%         density = nnz(mat_cells{1,i})/n(i)^2;
%         disp(['M=' num2str(M(i)) ' ; density=' num2str(density)]);
%     end
% end
% save(matrix_filename,'mat_cells','-v7.3');%save('Matrices','mat_cells')