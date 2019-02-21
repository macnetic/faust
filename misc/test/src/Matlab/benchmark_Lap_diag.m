% this script uses factorization algorithms from FAµST 1.03 even if it uses matfaust version >= 2 to create Faust
% the goal is to benchmark original algorithm implementations to compare their results to recent FAµST'
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

hier_palm_times = {};
hier_fgft_times = {};
givens_fgft_times = {};
par_givens_fgft_times = {};

num_graphs = 45; % only first 60 graph Laps which are of 6 types (random sensor, etc.)
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

	tic
	[lambda, facts, errors] = hierarchical_fact(params);
	F = Faust(facts, lambda);
	hier_palm_times = [ hier_palm_times toc ];
	complexity_global = nnz_sum(F);
	RC_PALM = complexity_global/nnz(U);

	tic
	[lambda, facts, Dhat, errors] = hierarchical_fact_FFT(params);
	facts{1} = lambda*facts{1};
	Uhat_PALM = Faust(facts);
	hier_fgft_times = [ hier_fgft_times toc ];
	hier_err = norm(Uhat_PALM*Dhat*Uhat_PALM' - Lap, 'fro')/norm(Lap, 'fro');
	hier_errs = [ hier_errs hier_err];

	U_hier_errs = [ U_hier_errs norm(Uhat_PALM - U, 'fro')/norm(U, 'fro')];


	J=round(RC_PALM*nnz(U)/4);
	tic
	[facts_givens,Dhat,err,L,choices] = diagonalization_givens(Lap,J);
	givens_fgft_times = [ givens_fgft_times toc];
	Uhat_givens = Faust(facts_givens);
	givens_err = norm(Uhat_givens*full(Dhat)*Uhat_givens' - Lap, 'fro')/norm(Lap, 'fro');
	givens_errs = [ givens_errs givens_err];
	U_givens_errs = [ U_givens_errs norm(Uhat_givens - U, 'fro')/norm(U, 'fro')];

	tic
	[facts_givens_parall,Dhat,err,L,coord_choices] = diagonalization_givens_parall(Lap,J,dim/2);
	par_givens_fgft_times = [ par_givens_fgft_times toc ];
	Uhat_givens_par = Faust(facts_givens_parall);
	par_givens_err = norm(Uhat_givens_par*full(Dhat)*Uhat_givens_par' - Lap, 'fro')/norm(Lap, 'fro');
	par_givens_errs = [ par_givens_errs par_givens_err];
	U_par_givens_errs = [ U_par_givens_errs norm(Uhat_givens_par - U, 'fro')/norm(U, 'fro')];

end

save('benchmark_Lap_diag_output.mat', 'hier_errs', 'givens_errs', 'par_givens_errs', 'U_hier_errs', 'U_givens_errs', 'U_par_givens_errs', 'hier_fgft_times', 'givens_fgft_times', 'par_givens_fgft_times', 'hier_palm_times')

hier_errs
givens_errs
par_givens_errs

U_hier_errs
U_givens_errs
U_par_givens_errs

hier_palm_times
hier_fgft_times
givens_fgft_times
par_givens_fgft_times
