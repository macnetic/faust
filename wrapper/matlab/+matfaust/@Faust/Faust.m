%% class FAUST
% represents a given dense matrix by a product of sparse matrix (i.e Faust)
% in order to speed-up multiplication by this matrix,
% Matlab wrapper class implemented in C++
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://Faust.gforge.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot, Adrien Leman, Thomas Gautrais, Luc Le Magoarou, Remi Gribonval
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
% along with F program.  If not, see <http://www.gnu.org/licenses/>.
%
%% Contacts:
%   Nicolas Bellot	: nicolas.bellot@inria.fr
%   Adrien Leman	: adrien.leman@inria.fr
%   Thomas Gautrais : thomas.gautrais@inria.fr
%   Luc Le Magoarou	: luc.le-magoarou@inria.fr
%   Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
%	approximations of matrices and applications", Journal of Selected
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%

% ======================================================================
%> @brief FAµST Matlab wrapper main class.
%>
%> This class provides a Matlab array-like interface for operations with FAµST data structures, which correspond to matrices that can be written exactly as the product of sparse matrices.
%>
%>A FAµST data structure is designed to allow fast matrix-vector multiplications together with reduced memory storage compared to what would be obtained by manipulating directly the corresponding (dense) Matlab array.
%>
%>A particular example is the matrix associated to the discrete Fourier transform, which can be represented exactly as a FAµST, leading to a fast and compact implementation.
%>
%> Although the sparse matrices are more interesting for optimization it's not forbidden to define a Faust as a product of dense matrices or a mix of dense and sparse matrices.
%>
%> The matrices composing the Faust product, also called the factors, are defined on complex or real fields. Hence a Faust can be a complex Faust or a real Faust.
%>
%> A Faust has ideally to be seen and used as a Matlab native dense matrix, but this matrix exists only virtually and is actually represented by its factors.
%> In order to use a Faust like a Matlab matrix, a certain number of Matlab built-ins available for Matlab matrices are implemented in this class but take note that not all are (e.g. IMAG() and REAL() are not defined for a Faust)
%>
%> It's possible to retrieve the Matlab native dense matrix with the method Faust.full but it will cost the multiplication of the Faust's factors.
%> It's noteworthy that in this documentation the expression 'dense matrix' designates the Matlab native dense matrix corresponding to a Faust, that is the matrix obtained by the multiplication of the previously mentioned Faust's factors.
%>
%> Likewise, other Faust's methods need to calculate the factor product, and they will be indicated with a warning in this documentation. You should avoid to use them if it's not really necessary (for example you might limit their use to test purposes).
%>
%> Other important limitation is that contrary to a Matlab native dense matrix a Faust is immutable. It means that you cannot affect elements of a Faust using the affectation operator `=' like you do with a Matlab matrix (e.g. `M(i,j) = 2').
%> That limitation is the reason why the Matlab built-in `SUBSASGN()' is not implemented in this class.
%>
%> For more information about FAµST take a look at http://faust.inria.fr.
%>
% ======================================================================

classdef Faust
	properties (SetAccess = private, Hidden = true)
		matrix; % Handle to the FaustCore class instance
		isReal;
	end
%	properties (Constant)
%	end
	methods
		%======================================================================
		%> @brief Creates a Faust from a list of factors or alternatively from a file.
		%>
		%> Another easy way to create a Faust is to call the static method FaustFactory.rand().
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; Faust(factors) creates a Faust from a list of factors (1D cell array).<br/><br>
		%> &nbsp;&nbsp;&nbsp; Faust(factors, scale) same as above but multiplying the Faust with a scalar.<br/><br/>
		%> &nbsp;&nbsp;&nbsp; Faust(filepath) creates a Faust from a previously saved Faust file (filepath is a character array).<br/><br/>
		%> &nbsp;&nbsp;&nbsp; Faust(filepath, scale) save as above but multiplying the Faust with a scalar.
		%>
		%>
		%> @param factors (varargin{1}) the 1D cell array of factors to initialize the Faust with.
		%> <br/> The factors must respect the dimensions needed for the product to be defined (for i=1 to size(factors,2), size(factors{i},2) == size(factors{i+1},1)).
		%> <br/> The factors can be sparse or dense matrices.
		%> @param filepath (varargin{1}) the file from which a Faust is created. It must be a character array.<br/>
		%>								The format is Matlab version 5 (.mat extension).<br/>
		%>								The file must have been saved before with Faust.save().
		%> @param scale (optional varargin{2}) a multiplicative scalar (see examples below).
		%>
		%> @b Examples
		%> @code
		%>	import matfaust.FaustFactory
		%>	factors = cell(1,5)
		%>	is_sparse = false
		%>	for i=1:5
		%>		if(is_sparse) % odd index factors are sparse matrices
		%>			factors{i} = sprand(100, 100, 0.1)
		%>		else % even index gives a dense matrix
		%>			factors{i} = rand(100, 100)
		%>		end
		%>		is_sparse = ~ is_sparse
		%>	end
		%>	% define a Faust with those factors
		%>	F = Faust(factors)
		%>
		%>	scale = 2
		%>	G = Faust(factors, scale) % G == scale*F
		%>
		%>	save(F, 'F.mat')
		%>	% define a Faust from file
		%>	H = Faust('F.mat')
		%>	I = Faust('F.mat', scale) % I == scale*H
		%>
		%> @endcode
		%>
		%>
		%> <p>@b See @b also Faust.delete, Faust.save, FaustFactory.rand</p>
		%>
		%======================================================================
		function F = Faust(varargin)
			%% FAUST Constructor - Creates a Faust from various types of input.
			%
			% Examples :
			%
			% F = matfaust.Faust(factors,scale)
			% - factor : 1D cell array of matrix (sparse or
			% dense) representing the factor of the Faust.
			% - scale : (optional) multiplicative scalar.
			%
			% F = matfaust.Faust(filepath, scale)
			% - filepath: the file where a Faust was stored with Faust.save() (in matlab format version 5).
			% - scale: (optional) multiplicative scalar.
			err_msg = 'matfaust.Faust() error: the arguments are not valid.';
			if(nargin < 1 || nargin > 3)
				error([err_msg ' Number of arguments passed is zero or greater than three.'])
			elseif(iscell(varargin{1}))
				% init a Faust from factor list
				% check if the factors are real or complex, one complex factor implies a complex faust
				factors=varargin{1};
				isRealFlag = 1;
				for i=1:length(factors)
					if (~isreal(factors{i}))
						isRealFlag = 0;
						break
					end
				end
				if(nargin == 2 && (~isnumeric(varargin{2}) || ~isscalar(varargin{2})))
					error([err_msg ' The second argument scale must be a scalar.'])
				end
				F.matrix = FaustCore(varargin{:});
				F.isReal = isRealFlag;
			elseif(ischar(varargin{1}))
				% init a Faust from file
				filename=varargin{1};
				load(filename);
				if (~exist('faust_factors','var') )
					error('Faust : invalid file');
				end
				F = matfaust.Faust(faust_factors, varargin{2:end});
			elseif(isa(varargin{1}, 'matfaust.Faust'))
				% create a Faust from another one but not with the same
				% handle to set inside the FaustCore object (matrix)
				oF = varargin{1};
				F.matrix = FaustCore(varargin{2}, oF.isReal);
				F.isReal = oF.isReal;
			elseif(isa(varargin{1}, 'integer') && islogical(varargin{2}))
				% create a Faust directly with the c++ object handler/pointer without any pre-existing Faust
				F.matrix = FaustCore(varargin{1}, varargin{2});
				F.isReal = varargin{2};
			else
				error(err_msg)
			end
		end


		%======================================================================
		%> @brief Deletes the Faust object (destructor).
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b delete(F)
		%>
		%> @param F the Faust to delete.
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand([2, 5], [50, 100], .5)
		%>	delete(F)
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust</p>
		%>
		%======================================================================
		function delete(F)
			%% DELETE Destructor delete the Faust.
			% delete(F)
			%
			% See also Faust
			delete(F.matrix)
		end

		%======================================================================
		%> @brief Multiplies the Faust F by A which is a dense matrix or a Faust object.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b B = mtimes(F, A)<br/>
		%> &nbsp;&nbsp;&nbsp; @b B = F*A
		%>
		%>
		%> @b NOTE: The primary goal of this function is to implement “fast” multiplication by a
		%> MuST, which is the operation performed when A is a standard matrix (dense or
		%> sparse).
		%> In the best cases, F*A is rcg(F) times faster than performing the
		%> equivalent full(F)*A.
		%> <br/>When A is a MuST, F*A is itself a MuST.
		%>
		%> @param F the Faust object.
		%> @param A The dense matrix to multiply or a Faust object.
		%>
		%> @retval B The multiplication result (a dense matrix or a Faust object depending on what A is).
		%>
		%> @b Example
		%> @code
		%>   import matfaust.FaustFactory
		%>   F = FaustFactory.rand([2, 5], [50, 100], .5)
		%>   A = rand(size(F,2), 50)
		%>   B = F*A
		%> % is equivalent to B = mtimes(F, A)
		%>   G = FaustFactory.rand(5,size(F, 2))
		%>   H = F*G
		%> % H is a Faust because F and G are
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust, Faust.rcg.
		%>
		%======================================================================
		function B = mtimes(F,A)
			%% MTIMES * Faust Multiplication (overloaded Matlab built-in function).
			%
			% B=mtimes(F,A) is called for syntax 'B=F*A', when F is a Faust matrix and A a full
			% storage matrix, B is also a full matrix storage.
			%
			% See also mtimes_trans
			B = mtimes_trans(F, A, 0);
		end


		%======================================================================
		%> @brief Multiplies the Faust F by A.' or A which is a dense matrix or a Faust object.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b WARNING: if A is a matrix the function costs get_num_factors(F) matrix multiplications.
		%> In that case its use is discouraged except for test purpose. However if A is a Faust object,
		%> it costs the same that a Faust initialization with a number of factors equal to
		%> F.get_num_factors()+A.get_num_factors() (like you can do directly with Faust.Faust).
		%>
		%> @param F the Faust object.
		%> @param A The dense matrix to multiply or a Faust object.
		%> @param trans equals 1 to calculate C=F'*A
		%> 				or 0 to calculate C=F*A.
		%>
		%> @retval C The multiplication result (a dense matrix or a Faust object depending on what A is).
		%>
		%> <p> @b See @b also mtimes.
		%======================================================================
		function C = mtimes_trans(F,A,trans)
			%% MTIMES_TRANS Multiplication by a Faust or its non-conjugate transposed.
			%
			% C = mtimes_trans(F,A,trans) when F is a Faust,A a full storage
			% matrix and trans a parameter, C a full storage matrix
			% if trans == 0, C=F*A is performed  (multiplication)
			% if trans == 1, C=F'*A is performed (multiplication by transposed)
			%
			% See also mtimes.

			if ~isreal(trans)
				error('invalid argument trans, must be equal to 0 or 1');
			end

			if (trans ~= 1) && (trans ~= 0)
				error('invalid argument trans, must be equal to 0 or 1');
			end
			if(issparse(A))
				error('Faust multiplication to a sparse matrix isn''t supported.')
			elseif(isa(A,'matfaust.Faust'))
				if (F.isReal)
					C = matfaust.Faust(F, mexFaustReal('mul_faust', F.matrix.objectHandle, A.matrix.objectHandle));
				else
					C = matfaust.Faust(F, mexFaustCplx('mul_faust', F.matrix.objectHandle, A.matrix.objectHandle));
				end
			elseif (F.isReal)
				if (isreal(A))
					C = mexFaustReal('multiply', F.matrix.objectHandle, A, trans);
				else
					C_real = mexFaustReal('multiply', F.matrix.objectHandle, real(A), trans);
					C_imag = mexFaustReal('multiply', F.matrix.objectHandle, imag(A), trans);
					C = C_real + 1i * C_imag;
				end
			else
				C = mexFaustCplx('multiply', F.matrix.objectHandle,A, trans);
			end
		end

		%======================================================================
		%> @brief Converts the Faust to a full matrix.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b A = full(F)
		%>
		%> @b Warning: using this function is discouraged except for test purposes,
		%> as it loses the main potential interests of the FAuST structure:
		%> compressed memory storage and faster matrix-vector multiplication compared
		%> to its equivalent dense matrix representation.
		%>
		%>
		%> @retval A the dense matrix resulting from the Faust. A is such that A*x == F*x
		%> for any vector x.
		%>
		%> @b Example
		%> @code
		%>   % in a matlab terminal
		%>   >> import matfaust.FaustFactory
		%>   >> F = FaustFactory.rand(5, 10^6, 10^-5, 'sparse')
		%>   Faust size 1000000x1000000, density 5e-05, nnz_sum 49999995, 5 factor(s):
		%>   - FACTOR 0 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   - FACTOR 1 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   - FACTOR 2 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   - FACTOR 3 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   - FACTOR 4 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   >> % an attempt to convert F to a dense matrix is most likely to raise a memory error
		%>   >> % the sparse format is the only way to handle a so big Faust
		%>   >> full(F)
		%>   Out of Memory
		%> @endcode
		%>
		%======================================================================
		function A = full(F)
			%% FULL  Convert Faust matrix to full matrix (overloaded Matlab
			% built-in function).
			%
			% A=full(F) converts a Faust matrix F to full storage matrix A.
			if (F.isReal)
				A=mexFaustReal('full', F.matrix.objectHandle);
			else
				A=mexFaustCplx('full', F.matrix.objectHandle);
			end

		end

		%======================================================================
		%> @brief Indicates if F is a real Faust or a complex Faust.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b NOTE: if F is a real Faust, all factors are real.
		%> If F is a complex Faust, at least one factor is complex.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b A = isreal(F)
		%>
		%> @retval bool 1 if F is a real Faust, 0 if it's a complex faust.
		%>
		%======================================================================
		function bool = isreal(F)
			%% ISREAL True for real scalar Faust (overloaded Matlab built-in function).
			%
			% isreal(F) returns 1 if Faust F does not have an imaginary part
			% and 0 otherwise.

			bool=F.isReal;

		end

		%======================================================================
		%> @brief Returns the transpose of the Faust F.
		%>
		%> This function overloads a Matlab built-in function/operator.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp;  @b F_trans = @b transpose(F)<br/>
		%> &nbsp;&nbsp;&nbsp;  @b F_trans = @b F.'
		%> @param F the Faust object.
		%>
		%> @retval F_trans F transpose as a Faust object.
		%>
		%>
		%> @b Example
		%> @code
		%>   F_trans = F.'
		%> % is equivalent to
		%>   F_trans = transpose(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.ctranspose
		%======================================================================
		function F_trans=transpose(F)
			%% TRANSPOSE .' Non-conjugate transposed Faust (overloaded Matlab built-in function).
			%
			% F_trans = transpose(F) is called for the syntax F.' when F is Faust.
			%
			% If Faust is a real matrix, the conjugate transposition will be the same as the real one
			%
			% See also ctranspose.
			%F_trans=F; % trans and F point share the same C++ underlying object (objectHandle)
			%F_trans.transpose_flag = xor(1,F.transpose_flag); % inverse the transpose flag
			if (F.isReal)
				F_trans = matfaust.Faust(F, mexFaustReal('transpose', F.matrix.objectHandle));
			else
				F_trans = matfaust.Faust(F, mexFaustCplx('transpose', F.matrix.objectHandle));
			end
		end

		%======================================================================
		%> @brief Returns the conjugate transpose of F.
		%>
		%> This function overloads a Matlab built-in function/operator.
		%>
		%> @param F the Faust object.
		%>
		%> @retval F_ctrans the conjugate transpose of F as a Faust object.
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand(5, [50, 100], .5, 'mixed', false)
		%>	F_ctrans = F'
		%>	F_ctrans2 = ctranspose(F)
		%>	% F_ctrans == F_ctrans2
		%>	F_ctrans2 = transpose(F)
		%>	F_ctrans2 = conj(F_ctrans2)
		%>	% F_ctrans == F_ctrans2
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.transpose, Faust.conj
		%>
		%======================================================================
		function F_ctrans = ctranspose(F)
			%% CTRANSPOSE ' Complex conjugate transposed Faust (overloaded Matlab built-in function).
			%
			% F_trans = ctranspose(F) is called for syntax F' (complex conjugate transpose) when F is a Faust.
			%
			%
			% See also transpose.
			if (isreal(F))
				F_ctrans=transpose(F);
			else
				F_ctrans = matfaust.Faust(F, mexFaustCplx('ctranspose', F.matrix.objectHandle));
			end
		end

		%======================================================================
		%> @brief Returns the conjugate of the Faust.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b F_conj = conj(F)
		%>
		%> @param F the Faust object.
		%>
		%> @retval F_conj a Faust object.
		%>  <br/> If F is a real Faust then F_conj == F.
		%>  <br/> if F is a complex Faust, the Faust object F_conj returned verifies the next assertion for all i=1:get_num_factors(F):
		%> @code
		%> conj(get_factor(F,i)) == get_factor(F_conj,i)
		%> @endcode
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand(5, [50, 100], .5, 'mixed', false)
		%>	F_conj = conj(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.get_factor, Faust.get_num_factors, Faust.ctranspose
		%>
		%======================================================================
		function F_conj = conj(F)
			%% CONJ ' Complex conjugate Faust (overloaded Matlab built-in function).
			%
			%  F_conj = conj(F) For a complex F, F_conj == REAL(F) - i*IMAG(F)
			%
			%
			if (F.isReal)
				F_conj = matfaust.Faust(F, mexFaustReal('conj', F.matrix.objectHandle));
			else
				F_conj = matfaust.Faust(F, mexFaustCplx('conj', F.matrix.objectHandle));
			end
		end


		%======================================================================
		%> @brief Gives the size of the Faust.
		%>
		%> The size is a pair of numbers: the number of rows and the number of columns
		%> of the equivalent dense matrix of F.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp;  [@b NROWS,@b NCOLS] = @b size(F)<br/>
		%> &nbsp;&nbsp;&nbsp; @b N = @b size(F,DIM) with N being the size of the DIM-th dimension of F.<br/>
		%> &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; In other words N= NROWS if DIM == 1, N = NCOLS if DIM == 2.
		%>
		%> @param F the Faust object.
		%> @param varargin can be missing or specifying the index of the dimension to get the size of.
		%>
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand(5, [50, 100], .5, 'mixed', false)
		%>	[nrows, ncols] = size(F)
		%>	nrows = size(F, 1)
		%>	ncols = size(F, 2)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.nnz_sum, Faust.numel
		%======================================================================
		function varargout = size(F,varargin)
			%% SIZE Size of a Faust (overloaded Matlab built-in function).
			%
			% D = size(F), for a Faust F, returns the two-element row vector
			% D = [M,N] containing the number of rows and columns in the Faust.
			%
			% M = size(F,DIM) returns the length of the dimension specified
			% by the scalar DIM.  For example, size(X,1) returns the number
			% of rows and size(F,2) returns the number of columns in the Faust.
			% If DIM > 2, M will be 1.



			nb_input = length(varargin);


			if (nb_input > 1)
				error('Too many input arguments');
			end

			if ((nb_input == 1) && (nargout > 1) | (nargout > 2))
				error('Too many output arguments');
			end
			Size=[-1 -1];
			if (F.isReal)
				Size=mexFaustReal('size',F.matrix.objectHandle);
			else
				Size=mexFaustCplx('size',F.matrix.objectHandle);
			end


			if (nb_input~=0)
				dimension_arg=varargin{1};
				if (floor(dimension_arg) ~= dimension_arg)
					error('Dimension argument must be a positive integer scalar within indexing range');
				end

				if (varargin{1}==1)
					Size=Size(1);
				elseif (varargin{1}==2)
					Size=Size(2);
				else
					Size=1;
				end

			end


			if (nargout < 2)
				varargout{1}=Size;
			else
				varargout{1}=Size(1);
				varargout{2}=Size(2);
			end
		end

		%======================================================================
		%> @brief Gives the number of elements in the Faust (equivalent to prod(size(F)).
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b n = @b numel(F)
		%>
		%> @param F the Faust object.
		%>
		%> @retval n the number of elements of the Faust.
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand(5, [50, 100], .5, 'mixed', false)
		%>	n = numel(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.size
		%>
		%======================================================================
		function n = numel(F)
			n = prod(size(F));
		end

		%======================================================================
		%> @brief Serves as the last index when slicing or indexing a Faust.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b A = full(F)
		%>
	 	%> @b Example
		%> @code
		%>	% in a matlab terminal
		%>	>> import matfaust.FaustFactory
		%>	>>
		%>	>> F = FaustFactory.rand([4, 5], 3, .9);
		%>	>> full(F)
		%>	ans =
		%>
		%>		-0.1006   -0.2041   -0.1878
		%>		-0.1382    0.0400   -0.0954
		%>		-0.1345   -0.1223   -0.1667
		%>
		%>	>> full(F(2:end,1))
		%>	ans =
		%>
		%>		-0.1382
		%>		-0.1345
		%>
		%>	>> full(F(1,1:2:end)
		%>	ans =
		%>
		%>	-0.1006   -0.1878
		%>
		%>	>> full(F(1,1:end-1))
		%>	ans =
		%>
		%>		-0.1006   -0.2041
		%>
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.subsref, Faust.size, Faust.full
		%>
		%======================================================================
		function end_dim = end(F,k,n)
			%% END (useful for slicing) serve as the last index in an indexing expression (overloaded Matlab built-in function).
			%
			% Examples of use for slicing a Faust F are
			% F(3:end,1) : in this case, end=size(F,1)
			%   i.e end equals to the number of row of the Faust F.
			% F(1,1:2:end-1) : in this case, end=size(F,2)
			% end equals to the number of column fo the Faust F.
			%
			% See also subsref, size.

			if (n ~= 2)
				error('invalid slicing : Faust is a 2D array i.e matrix');
			end

			end_dim=size(F,k);


		end




		%=====================================================================
		%> @brief Returns the i-th factor of F.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b factor = get_factor(F, i)
		%>
		%> @param F the Faust object.
		%> @param i the factor index.
		%>
		%> @retval factor the i-th factor as a dense matrix.
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand(5, [50, 100], .5, 'mixed', false)
		%>	f1 = get_factor(F, 1)
		%> @endcode
		%> <p>@b See @b also Faust.get_num_factors
		%=====================================================================
		function factor = get_factor(F, i)
			%% GET_FACT Ith factor of the Faust.
			%
			% A=get_factor(F,i) return the i factor A of the Faust F as a full storage matrix.
			%
			% Example of use :
			% A=get_factor(F,1) returns the 1st factor of the Faust F.
			% A=get_factor(F,4) returns the 4th factor of the Faust F.
			%
			% See also get_num_factors.

			if (~isa(i,'double'))
				error('get_fact second argument (indice) must either be real positive integers or logicals.');
			end

			if (floor(i) ~= i)
				error('get_fact second argument (indice) must either be real positive integers or logicals.');
			end

			if (F.isReal)
				factor = mexFaustReal('get_fact',F.matrix.objectHandle,i);
			else
				factor = mexFaustCplx('get_fact',F.matrix.objectHandle,i);
			end

		end


		%==========================================================================================
		%> @brief Gives the number of factors of F.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b A = get_num_factors(F)
		%>
		%> @param F the Faust object.
		%>
		%> @retval num_factors the number of factors.
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand(5, [50, 100], .5, 'mixed', false)
		%>	nf = get_num_factors(F)
		%> @endcode
		%>
		%> <p>@b See @b also Faust.get_factor.
		%==========================================================================================
		function num_factors = get_num_factors(F)
			%% GET_NB_FACTOR Number of factor of the Faust.
			%
			% num_factors = get_num_factors(F) return the number of factor of the
			% Faust F.
			%
			% See also get_factor.
			if (F.isReal)
				num_factors = mexFaustReal('get_nb_factor', F.matrix.objectHandle);
			else
				num_factors = mexFaustCplx('get_nb_factor', F.matrix.objectHandle);
			end
		end

		%==========================================================================================
		%> @brief Saves the Faust F into file respecting the Matlab format version 5 (.mat file).
		%>
		%> @b NOTE: storing F should typically use rcg(F) times less storage than storing full(F).
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b save(F, filepath)
		%>
		%> @param F the Faust object.
		%> @param filepath the path for saving the Faust. The filename should end with .mat and it must be a character array.
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand(5, [50, 100], .5, 'mixed', false)
		%>	save(F, 'F.mat')
		%>	G = matfaust.Faust('F.mat')
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust, Faust.rcg.
		%==========================================================================================
		function save(F, filepath)
			%% save Saves a Faust into a matfile.
			%
			%  save(F,filepath) saves the Faust F into the .mat file specified by filepath.
			if(F.isReal)
				mexFaustReal('save', F.matrix.objectHandle, filepath)
			else
				mexFaustCplx('save', F.matrix.objectHandle, filepath)
			end
		end

		%==========================================================================================
		%> @brief Returns a Faust representing a submatrix of the dense matrix full(F) or a
		%>         scalar element if that Faust can be reduced to a single element.
		%>
		%>
		%> This function is a Matlab built-in overload.
		%>
		%> @b WARNING:
		%> - This function doesn't handle a slice step different to 1 (e.g. F(i:2:j,:) where slice step is 2.)
		%> - It is not advised to use this function as an element accessor
		%>        (e.g. F(0,0)) because such an use induces to convert the Faust to its
		%>        dense matrix representation and that is a very expensive computation if used
		%>        repetitively.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b F(I,J) extracts the Faust representing the sub-matrix of full(F) defined by the subscript vectors I and J (take a look to examples below).
		%>
		%> @param F the Faust object.
		%> @param S the structure defining the Faust to extract like a submatrix (it's not supposed to be used directly ; usage and examples show how subscript should be used).
		%>
		%> @retval The Faust object requested or just the corresponding scalar if that Faust has
		%>             a size equal to [1,1].
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>
		%>	F = FaustFactory.rand([2, 5], [50, 100], .5)
		%>	i = randi(min(size(F)), 1, 2)
		%>	i1 = i(1);i2 = i(2)
		%>
		%>	F(i1,i2) % is the scalar element located
		%>			 % at row i1, column i2 of the F's dense matrix
		%>
		%>	F(:,i2) % full column i2
		%>
		%>	F(3:4,2:5) % from row 3 to line 4, each row containing only elements from column 2 to 5.
		%>
		%>	F(1:end,5:end-1)  % from row 1 to end row, each one containing only elements from column 5 to column before the last one.
		%> @endcode
		%>
		%> <p>@b See @b also Faust.end.
		%==========================================================================================
		function submatrix = subsref(F,S)
			%% SUBSREF Subscripted reference (overloaded Matlab built-in function).
			%
			%  F(I,J) is an array formed from the elements of the rectangular
			% submatrix of the Faust F specified by the subscript vectors I and J.  The
			% resulting array has LENGTH(I) rows and LENGTH(J) columns.  A colon used
			% as a subscript, as in F(I,:), indicates all columns of those rows
			% indicated by vector I. Similarly, F(:,J) = B means all rows of columns
			%J.
			%
			% Example of use :
			%  A(i,j) A(:,j)  A(3:4,2:5) A(1:end,5:end-1)
			%
			% See also end.

			if (~isfield(S,'type')) | (~isfield(S,'subs'))
				error(' subsref invalid structure S missing field type or subs');
			end

			if (~ischar(S.type)) | (~iscell(S.subs))
				error(' subsref invalid structure S, S.type must be a character array, S.subs must be a cell array');
			end

			if ~strcmp(S.type,'()')
				%submatrix = builtin('subsref', F, S)
				error(' subsref is only overloaded for () operator');
			end

			if (length(S.subs) ~=2)
				invalid(' subsref invalid slicing must have 2 index since F is a 2D-array');
			end



			slicing_row=S.subs{1};
			slicing_col=S.subs{2};

			[dim1 dim2]=size(F);

			if ischar(slicing_row)
				start_row_id = 1;
				end_row_id = dim1;
			else
				start_row_id = slicing_row(1);
				end_row_id = slicing_row(end);
			end

			if ischar(slicing_col)
				start_col_id = 1;
				end_col_id = dim2;
			else
				start_col_id = slicing_col(1);
				end_col_id = slicing_col(end);
			end

			if(start_row_id < 0 || end_row_id < 0 || end_col_id < 0 || start_col_id < 0)
				error(' Subscript indices must be positive integers.')
			end

			if(start_row_id > dim1 || end_row_id > dim1 || end_col_id > dim2 || start_col_id > dim2)
				error(' Index exceeds Faust dimensions.')
			end

			if(F.isReal)
				submatrix = matfaust.Faust(F, mexFaustReal('subsref', F.matrix.objectHandle, start_row_id, end_row_id, start_col_id, end_col_id));
			else
				submatrix = matfaust.Faust(F, mexFaustCplx('subsref', F.matrix.objectHandle, start_row_id, end_row_id, start_col_id, end_col_id));
			end

			if(start_row_id == end_row_id && start_col_id == end_col_id)
				submatrix = full(submatrix);
				submatrix = submatrix(1,1);
			end

		end

		%======================================================================
		%> @brief @b WARNING this function is not implemented because a Faust object is immutable.
		%>
		%> This function overloads a Matlab built-in function.
		%> This function just throws an error.
		%>
		%>
		%> @param F the Faust to display information about.
		%>
		%>
		%======================================================================
		function F = subsasgn(F,S,B)
			%% SUBSASGN (WARNING not implemented) (overloaded Matlab built-in function)
			%
			% This function is not available for Faust class, because a Faust is immutable.
			% This function just throws an error.
			%
			% F(i,j)=1, F(2:5,3:5)=zeros(4,3) will throw
			% a Matlab error with this message :
			% 'function not implemented for Faust class'
			error('Function not implemented for Faust class');
		end

		%======================================================================
		%> @brief Displays information about F.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b disp(F) <br/><br/>
		%> &nbsp;&nbsp;&nbsp; @b F
		%>
		%> @param F the Faust object.
		%>
		%> @b Example
		%> @code
		%>	% in a matlab terminal
		%>	>> import matfaust.FaustFactory
		%>	>>
		%>	>> F = FaustFactory.rand([1, 2], [50, 100], .5)
		%>	>> disp(F)
		%>	Faust size 98x82, density 0.686909, nnz_sum 5520, 2 factor(s):
		%>	- FACTOR 0 (real) SPARSE, size 98x78, density 0.395081, nnz 3020
		%>	- FACTOR 1 (real) SPARSE, size 78x82, density 0.390869, nnz 2500
		%>
		%>	>> F
		%>	Faust size 98x82, density 0.686909, nnz_sum 5520, 2 factor(s):
		%>	- FACTOR 0 (real) SPARSE, size 98x78, density 0.395081, nnz 3020
		%>	- FACTOR 1 (real) SPARSE, size 78x82, density 0.390869, nnz 2500
		%>
		%> @endcode
		%>
		%> <p>@b See @b also Faust.nnz_sum, Faust.density, Faust.size, Faust.get_factor, Faust.get_num_factors
		%>
		%>
		%======================================================================
		function disp(F)
			%% DISP shows the characteristic of the Faust (overloaded Matlab built-in function)
			%
			%
			% This function shows the size of the Faust,
			%  its number of factor, its RCG ...
			%
			if (F.isReal)
				mexFaustReal('disp',F.matrix.objectHandle);
			else
				mexFaustCplx('disp',F.matrix.objectHandle);
			end


		end

		%======================================================================
		%> @brief Computes the norm of F.
		%>
		%> Several types of norm are available: 1-norm, 2-norm and Frobenius norm.
		%>
		%> The norm of F is equal to the norm of its dense matrix full(F).
		%>
		%> @b WARNING: this function costs at least as much as Faust.mtimes.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b n = @b norm(F, 2) returns the 2-norm of F: approximately max(svd(full(F))).<br/><br>
		%> &nbsp;&nbsp;&nbsp; @b n = @b norm(F) is the same as norm(F, 2).<br/><br>
		%> &nbsp;&nbsp;&nbsp; @b n = @b norm(F, 1) returns the 1-norm of F: the maximum absolute column sum of the matrix full(F).<br/><br>
		%> &nbsp;&nbsp;&nbsp; @b n = @b @b norm(F, @b 'fro') returns the Frobenius norm of F.<br/><br>
		%>
		%> @param F the Faust object.
		%> @param p (optional) the norm order or type. Respectively 1 or 2 for the 1-norm and 2-norm or 'fro' for the Frobenius norm (by default the 2-norm is computed).
		%>
		%>
		%> @retval n the norm (real).
		%>
		%>
		%> @b Example
		%> @code
		%> % in a matlab terminal
		%> >> import matfaust.FaustFactory
		%> >> F = FaustFactory.rand([1, 2], [50, 100], .5)
		%> >> norm(F)
		%> ans =
		%> 80.0914
		%> >> norm(F,2)
		%> ans =
		%> 80.0914
		%> >> norm(F, 'fro')
		%> ans =
		%> 133.9623
		%> @endcode
		%>
		%======================================================================
		function n = norm(F,varargin)
			%% NORM Faust norm (overloaded Matlab built-in function).
			%
			% norm(F,1) when F is a Faust returns L1 norm of F (the largest
			% column sum of the absolute values of F).
			% norm(F,2) when F is a Faust returns the L2 norm of F (the largest
			% singular value of A).
			% norm(F,'fro') when F is a Faust returns the Frobenius norm of F.
			% norm(F) is the same as norm(F,2)
			%
			% @b WARNING : norm(F,P) is only supported when P equals 1, 2 or
			% 'fro'.

			nb_input = length(varargin);
			if (nb_input > 1)
				error('Too many input arguments');
			end

			ord = 2;
			if nb_input == 1
				if(varargin{1} == 'fro')
					if (F.isReal)
						n = mexFaustReal('normfro',F.matrix.objectHandle);
					else
						n = mexFaustCplx('normfro',F.matrix.objectHandle);
					end
					return
				end
				if varargin{1} ~= 2 && varargin{1} ~= 1
					error('only L1, L2 or Frobenius norms are supported for the Faust');
				end
				ord = varargin{1};
			end

			if (F.isReal)
				n = mexFaustReal('norm',F.matrix.objectHandle, ord);
			else
				n = mexFaustCplx('norm',F.matrix.objectHandle, ord);
			end

		end

		%==========================================================================================
		%> @brief Gives the total number of non-zero elements in the factors of F.
		%>
		%> The function sums together the number of non-zeros elements of
		%> each factor and returns the result. Note that in fact the sum is
		%>                        computed at Faust creation time and kept in cache.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b nz = nnz_sum(F)
		%> @param F the Faust object.
		%>
		%> @retval nz The number of non-zeros.
		%>
		%> <p>@b See @b also Faust.rcg, Faust.density.
		%==========================================================================================
		function nz = nnz_sum(F)
			%% NNZ Number of nonzero elements in a Faust (overloaded Matlab built-in function).
			%
			% nz = nnz_sum(F) is the number of nonzero elements in the Faust F.
			%
			% See also density, rcg.
			if (F.isReal)
				nz=mexFaustReal('nnz',F.matrix.objectHandle);
			else
				nz=mexFaustCplx('nnz',F.matrix.objectHandle);
			end
		end

		%======================================================================
		%> @brief Calculates the density of F such that nnz_sum(F) == density(F)*numel(F).
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b dens = density(F)
		%>
		%>
		%> @b NOTE: a value of density below one indicates potential memory savings compared
		%> to storing the corresponding dense matrix full(F), as well as potentially
		%> faster matrix-vector multiplication when applying F*x instead of full(F)*x.
		%>
		%> @b NOTE: a density above one is possible but prevents any saving.
		%>
		%>
		%> @param F the Faust object.
		%>
		%> @retval dens density of F.
		%>
		%> @b Example
		%> @code
		%>	import matfaust.FaustFactory
		%>	F = FaustFactory.rand([2, 5], [50, 100], .5)
		%>	dens = density(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.nnz_sum, Faust.rcg, Faust.size, Faust.numel
		%======================================================================
		function dens = density(F)
			%% DENSITY Density of the Faust.
			%
			% dens = density(F) when F is a Faust returns the
			% percentage of nonzero elements of F,
			% dens is a number between 0 and 1.
			% In some degenerated case, dens can be greater than 1.
			% If the Faust is empty, return -1.
			%
			% See also rcg, nnz_sum.

			prod_dim=numel(F);
			if (prod_dim ~= 0)
				dens=nnz_sum(F)/prod_dim;
			else
				dens = -1;
			end
		end

		%==========================================================================================
		%> @brief Computes the Relative Complexity Gain (inverse of Faust.density).
		%>
		%> RCG is the theoretical gain brought by Faust representation relatively to its dense
		%> matrix equivalent. The higher is the RCG, the more computational
		%> savings will be made.
		%> That gain applies both for storage space and computation time.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b gain = rcg(F)
		%>
		%> @param F	the Faust object.
		%>
		%>
		%> @retval gain = the RCG value (real). rcg(F) == 1/density(F).
		%>
		%> <p>@b See @b also Faust.density, Faust.nnz_sum, Faust.size.
		%==========================================================================================
		function gain = rcg(F)
			%% RCG Relative Complexity Gain (inverse of the density)
			%
			% gain =  rcg(F) when F is Faust, returns the
			% inverse of density of the Faust (i.e the theoretical gain
			% both for storage and multiplication computation time between the Faust and its full storage
			% equivalent full(F)).
			%
			% See also density, nnz_sum, size.

			dens=density(F);
			if (dens > 0)
				gain = 1/dens;
			else
				if (dens == 0)
					gain = Inf;
				else
					gain = -1;
				end
			end
		end


	end
	methods(Static)
	end
end

