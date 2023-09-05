% assuming we are in faust projet root directory
% but that faust is installed on /opt/local/faust with a sys package (on MacOS or Linux)
%system('git clone https://github.com/MOcov/MOcov.git')
%cd MOcov
%system('make install')
%profile on
%addpath(pwd)
%cd ..
% doctest part
profile on
addpath '/usr/share/octave/packages/doctest-0.7.0/'
setup_FAUST
doctest matfaust.proj.sp
doctest matfaust.proj.splin
doctest matfaust.proj.spcol
doctest matfaust.proj.splincol
doctest matfaust.proj.const
doctest matfaust.proj.supp
doctest matfaust.proj.normlin
doctest matfaust.proj.normcol
doctest matfaust.proj.blockdiag
doctest matfaust.proj.anticirc
doctest matfaust.proj.circ
doctest matfaust.proj.hankel
doctest matfaust.proj.skperm
doctest matfaust.proj.sptriu
doctest matfaust.proj.sptril
doctest matfaust.proj.spsymm
doctest matfaust.proj.toeplitz
doctest matfaust.proj.proj_id
doctest matfaust.poly.poly
doctest matfaust.poly.expm_multiply
doctest matfaust.poly.next
doctest matfaust.poly.basis
doctest matfaust.poly.expm_inv_multiply
doctest matfaust.factparams.ParamsHierarchicalWHTNoResCons
doctest matfaust.factparams.ParamsHierarchicalRectMatNoResCons
doctest matfaust.factparams.ParamsHierarchicalNoResCons
doctest matfaust.factparams.ParamsHierarchicalRectMat
doctest matfaust.fact.butterfly
doctest matfaust.fact.check_fact_mat
doctest matfaust.fact.eigtj
doctest matfaust.fact.fact
doctest matfaust.fact.fgft_palm
doctest matfaust.fact.hierarchical_constends
% doctest matfaust.fact.hierarchical % too heavy
doctest matfaust.fact.hierarchical_mhtp
doctest matfaust.fact.palm4msa_constends
doctest matfaust.fact.palm4msa
doctest matfaust.fact.palm4msa_mhtp
doctest matfaust.fact.pinvtj
doctest matfaust.fact.svdtj
doctest matfaust.tools.omp
doctest matfaust.demo.runall
% unit test part
addpath([ pwd 'misc/test/src/Matlab/'])
cd misc/test/src/Matlab;
FaustTest '../../../../misc/test/src/'
FaustFactoryTest '../../../../misc/test/src/'
if any(strcmp(split(getenv('CI_RUNNER_TAGS'), ':'), 'cuda'))
    matfaust.rand(10, 10, 'dev', 'gpu')
end
% compute coverage
%mocov('-cover','/opt/local/faust/matlab/+matfaust',...
cd '../../../..'
mocov('-cover','/opt/local/faust/matlab/+matfaust',...
	'-profile_info',...
	'-cover_html_dir','coverage_html');

profile off
