%======================================================================
%> @brief This function is an helper to create a FAµST BSR matrix (kind of matrix that is not natively supported by Matlab).
%>
%> In fact it doesn't directly create a BSR matrix but rather an encoding to create the matrix later through the matfaust.Faust constructor.
%>
%> @param M The number of rows of the BSR matrix.
%> @param N The number of columns of the BSR matrix.
%> @param bnnz The number of nonzero blocks of the BSR matrix.
%> @param bdata The horizontal concatenation of the matrix nonzero blocks taken in rowmajor order. Note that all the nonzero blocks have the same size.
%> @param bcolinds The block column indices corresponding (in the same order) to the blocks in bdata. A block column index is not a column index, for example if your matrix contains N columns and each nonzero block contains BN columns, the block column indices must lay in {1, ..., N/BN}.
%> @param brow_count The number of nonzero blocks in each block-row of the matrix (a vector of size M/BM, where brow_count(i) is the number of nonzero blocks in the i-th block-row).
%>
%> @b Examples
%> @code
%> >> bnnz = 3;
%> >> M = 10;
%> >> N = 9;
%> >> nonzero_blocks = [[ 0.6787, 0.7431 0.6555; 0.7577, 0.3922 0.1712] [ 0.7060 0.2769    0.0971; 0.0318    0.0462    0.8235] [0.6948    0.9502    0.4387; 0.3171    0.0344    0.3816]];
%> >> bsr_mat = matfaust.create_bsr(M, N, bnnz, nonzero_blocks, [2, 3, 1], [2, 0, 1, 0, 0]);
%> >> F = matfaust.Faust(bsr_mat)
%>
%> F =
%>
%> Faust size 10x10, density 0.2, nnz_sum 18, 1 factor(s):
%> - FACTOR 0 (double) BSR, size 10x9, density 0.2, nnz 18
%> >> full(F)
%>
%>
%> ans =
%>
%>          0         0         0    0.6787    0.7431    0.6555    0.7060    0.2769    0.0971
%>          0         0         0    0.7577    0.3922    0.1712    0.0318    0.0462    0.8235
%>          0         0         0         0         0         0         0         0         0
%>          0         0         0         0         0         0         0         0         0
%>     0.6948    0.9502    0.4387         0         0         0         0         0         0
%>     0.3171    0.0344    0.3816         0         0         0         0         0         0
%>          0         0         0         0         0         0         0         0         0
%>          0         0         0         0         0         0         0         0         0
%>          0         0         0         0         0         0         0         0         0
%>          0         0         0         0         0         0         0         0         0
%>
%> @endcode
%======================================================================
function bsr_enc = create_bsr(M, N, bnnz, bdata, bcolinds, brow_count)
	[ok, err] = verify_M_N(M, N);
	if(~ ok) error(err); end
	[ok, err] = verify_bnnz(bnnz);
	if(~ ok) error(err); end
	[ok, err] = verify_bdata(bdata, M, N, bnnz);
	if(~ ok) error(err); end
	[ok, err] = verify_bcolinds(bcolinds, bnnz);
	if(~ ok) error(err); end
	[ok, err] = verify_brow_count(brow_count, M, bdata);
	if(~ ok) error(err); end
	bsr_enc = {'bsr', M, N, bnnz, bdata, int32(bcolinds), int32(brow_count)};
	verify_bsr_format(bsr_enc);
end

function [ok, err] = verify_M_N(M, N)
	ok = isnumeric(M) && isnumeric(N) && isreal(M) && isreal(N) && ismatrix(M) && ismatrix(N) && numel(M) == 1 && numel(N) == 1 && floor(M) == M && floor(N) == N;
	err = 'M and N must be integers (defining the respectively the number of rows and number of columns of the BSR matrix).';
end

function [ok, err] = verify_bnnz(bnnz)
	ok = isnumeric(bnnz) && isnumeric(bnnz) && ismatrix(bnnz) && isreal(bnnz) && numel(bnnz) == 1 && floor(bnnz) == bnnz;
	err = 'bnnz must be an integer (the number of nonzero blocks of the BSR matrix)';
end

function [ok, err] = verify_bdata(bdata, M, N, bnnz)
	ok = isnumeric(bdata) && ismatrix(bdata) && mod(M, size(bdata, 1)) == 0 && mod(N, size(bdata, 2)/bnnz) == 0;
	err = 'bdata must be a matrix whose the number of rows (which is the nonzero block number of rows) divides evenly the matrix number of rows. Its number of columns must be equal to the number of nonzero blocks times the nonzero block number of columns (bdata is the horizontal concatenation of the nonzero blocks).';
end

function [ok, err] = verify_bcolinds(bcolinds, bnnz)
	ok = isnumeric(bcolinds) && ismatrix(bcolinds) && isreal(bcolinds) && all(all(floor(bcolinds) == bcolinds)) && size(bcolinds, 1) == 1 && size(bcolinds, 2) == bnnz;
	err = 'bcolinds is the vector of block column indices, each of these indices being an integer between 1 and the number of nonzero blocks per row in the matrix. The size of the vector must be equal to bnnz (one index per nonzero block).';
end

function [ok, err] = verify_brow_count(brow_count, M, bdata)
	BM = size(bdata, 1);
	ok = isnumeric(brow_count) && ismatrix(brow_count) && isreal(brow_count) && all(all(floor(brow_count) == brow_count)) && size(brow_count, 1) == 1 && size(brow_count, 2) == M/BM;
	err = 'brow_count must be a vector containing the number of nonzero blocks on each block-row of the matrix.';
end

function ok = verify_bsr_format(c)
	ok = iscell(c);
	if(~ ok) error('A BSR encoding must be a cell array');end
	ok = length(c) == 7;
	if(~ ok) error('The length of the BSR cell must be 7');end
	ok = strcmp(lower(c{1}), 'bsr');
	if(~ ok) error('The BSR cell''s first entry must be the char array ''bsr''');end
	[ok, err] = verify_M_N(c{2}, c{3});
	if(~ ok) error(err); end
	[ok, err] = verify_bnnz(c{4});
	if(~ ok) error(err);end
	bnnz = c{4};
	[ok, err] = verify_bdata(c{5}, c{2}, c{3}, bnnz);
	if(~ ok) error(err); end
	[ok, err] = verify_bcolinds(c{6}, bnnz);
	if(~ ok) error(err); end
	[ok, err] = verify_brow_count(c{7}, c{2}, c{5});
	if(~ ok) error(err); end
end

