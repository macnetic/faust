import matfaust.Faust

fpath = mfilename('fullpath');
[path, ~, ~ ] = fileparts(fpath);
fs = filesep;
data_path = [ path fs '..' fs '..' fs '..' fs 'misc' fs 'data' fs 'mat' fs 'Laps_U_Ds-6_graph_types-dims_128_to_1024.mat'];
if(~ exist(data_path))
	error([data_path ' not found. Please run the data generation script misc/test/src/Matlab/get_graph_test_UDLap_matrices.m'])
end

load(data_path)
num_graphs = length(Laps);

hier_errs = {};
givens_errs = {};
par_givens_errs = {};

U_hier_errs = {};
U_givens_errs = {};
U_par_givens_errs = {};

num_graphs = 60; % only first 60 graph Laps which are of 6 types (random sensor, etc.)
for i=1:num_graphs
	disp(['benchmarking graph. Laplacian ' num2str(i) '/' num2str(num_graphs)])
	U = Us{i};
	D = Ds{i};
	Lap = Laps{i};

	dim = size(Lap, 1);
	nfacts = round(log2(dim))-3; % Desired number of factors
	over_sp = 1.5; % Sparsity overhead
	dec_fact = 0.5; % Decrease of the residual sparsity

	params.data = U;
	params.Lap = Lap;
	params.init_D = D;
	params.nfacts = nfacts;
	params.fact_side = 0;


	params.cons = cell(2,nfacts-1);
	for j=1:nfacts-1
		params.cons{1,j} = {'sp',dec_fact^j*dim^2*over_sp,dim,dim};
		params.cons{2,j} = {'sp',2*dim*over_sp,dim,dim};
	end

	% Number of iterations
	params.niter1 = 50;
	params.niter2 = 100;
	params.verbose=0;

	[lambda, facts, errors] = hierarchical_fact(params);
	complexity_global=0;
	for j=1:nfacts
		complexity_global = complexity_global + nnz(facts{j});
	end
	F = Faust(facts);
	RC_PALM = complexity_global/nnz(U);

	[lambda, facts, Dhat, errors] = hierarchical_fact_FFT(params);
	facts{1} = lambda*facts{1};
	Uhat_PALM = Faust(facts);
	hier_err = norm(Uhat_PALM*Dhat*Uhat_PALM' - Lap, 'fro')/norm(Lap, 'fro');
	hier_errs = [ hier_errs hier_err];

	U_hier_errs = [ U_hier_errs norm(Uhat_PALM - U, 'fro')/norm(U, 'fro')];


	J=round(RC_PALM*nnz(U)/4);
	[facts_givens,Dhat,err,L,choices] = diagonalization_givens(Lap,J);
	Uhat_givens = Faust(facts_givens);
	givens_err = norm(Uhat_givens*full(Dhat)*Uhat_givens' - Lap, 'fro')/norm(Lap, 'fro');
	givens_errs = [ givens_errs givens_err];
	U_givens_errs = [ U_givens_errs norm(Uhat_givens - U, 'fro')/norm(U, 'fro')];

	[facts_givens_parall,Dhat,err,L,coord_choices] = diagonalization_givens_parall(Lap,J,dim/2);
	Uhat_givens_par = Faust(facts_givens_parall);
	par_givens_err = norm(Uhat_givens_par*full(Dhat)*Uhat_givens_par' - Lap, 'fro')/norm(Lap, 'fro');
	par_givens_errs = [ par_givens_errs par_givens_err];
	U_par_givens_errs = [ U_par_givens_errs norm(Uhat_givens_par - U, 'fro')/norm(U, 'fro')];

end

save('benchmark_Lap_diag_output.mat', 'hier_errs', 'givens_errs', 'par_givens_errs', 'U_hier_errs', 'U_givens_errs', 'U_par_givens_errs')

hier_errs
givens_errs
par_givens_errs

U_hier_errs
U_givens_errs
U_par_givens_errs

