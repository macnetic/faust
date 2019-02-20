% code inspired from FAuST 1.03 /experiment_comparison_PALM_givens.m to extract graph tested matrices
% necessitates the GSP toolbox to generates graphs/Laplacians

DIMtab = [128,256,512,1024]; % Dimension of the graphs (number of nodes)
Nrun=10;

resTab = zeros(4, numel(DIMtab));

graph_types = {'erdos_renyi', 'community', 'sensor', 'path', 'random_ring', 'ring'}

Laps = {}
Us = {}
Ds = {}

for k = 1:numel(DIMtab)
	DIM = DIMtab(k);
	%disp(DIM)
	for n=1:Nrun
		%disp(n)
		for l = 1:6 % l == 3
			disp(['k=' num2str(k) ])%, n=' num2str(n)]) %)', l=' num2str(l)])

			%%%
			%% ==================== Graph generation =========================
			%%%
			if l==1
				%disp(['l = ' num2str(l)])
				G = gsp_erdos_renyi(DIM,0.1);
			elseif l==2
				%disp(['l = ' num2str(l)])
				G = gsp_community(DIM);
			elseif l==3
				%disp(['l = ' num2str(l)])
				G = gsp_sensor(DIM);
			elseif l==4
				%disp(['l = ' num2str(l)])
				G = gsp_path(DIM);
			elseif l==5
				%disp(['l = ' num2str(l)])
				G = gsp_random_ring(DIM);
			elseif l==6
				%disp(['l = ' num2str(l)])
				G = gsp_ring(DIM);
			end

			Lap = full(G.L);
			%save(['Laplacian_' num2str(DIMtab(k)) '_' graph_types{l}], 'Lap')
			Laps = [ Laps Lap ]
			%            G = gsp_compute_fourier_basis(G);
			%%             U = full(G.U);
			%%             D = diag(G.e);
			[U,D]=eig(Lap);
			Us = [ Us U ]
			Ds = [ Ds D ]
		end
	end
end


save('Laps_U_Ds-6_graph_types-dims_128_to_1024.mat', 'Laps', 'Us', 'Ds')
