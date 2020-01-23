%% class FAUST
% represents a given dense matrix by a product of sparse matrix (i.e Faust)
% in order to speed-up multiplication by this matrix,
% Matlab wrapper class implemented in C++
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://Faust.gforge.inria.fr>
%
%% License:
% Copyright (2018):	Hakim Hadj-djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais, Luc Le Magoarou, Remi Gribonval
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
%	Hakim Hadj-dji. : hakim.hadj-djilani@inria.fr
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

%> @package matfaust @brief <b> The FAuST Matlab Wrapper

% ======================================================================
%> @brief <b>FAuST Matlab wrapper main class</b> for using multi-layer sparse transforms.
%>
%> This class provides a Matlab array-like interface for operations with FAuST data structures, which correspond to matrices that can be written exactly as the product of sparse matrices.
%>
%> A FAuST data structure is designed to allow fast matrix-vector multiplications together with reduced memory storage compared to what would be obtained by manipulating directly the corresponding (dense) Matlab array.
%>
%> A particular example is the matrix associated to the discrete Fourier transform, which can be represented exactly as a FAuST, leading to a fast and compact implementation (see matfaust.dft()).
%>
%> Although sparse matrices are more interesting for optimization it's not forbidden to define a Faust as a product of dense matrices or a mix of dense and sparse matrices.
%>
%> The matrices composing the Faust product, also called the factors, are defined on complex or real fields. Hence a Faust can be a complex Faust or a real Faust.
%>
%> Several Matlab builtins have been overloaded to ensure that a Faust is
%> almost handled as a native Matlab matrix.
%>
%> The main exception is that contrary to a Matlab native array a Faust is immutable.
%> It means that you cannot affect elements of a Faust using
%> the affectation operator `=' like you do with a Matlab matrix (e.g. `M(i,j) =
%> 2').
%> That limitation is the reason why the Matlab built-in `SUBSASGN()' is not
%> implemented in this class.
%>
%> Other notable limitations are that one cannot:
%> - compute the real and imaginary parts of a Faust,
%> - perform elementwise operations between two Fausts (e.g. elementwise
%> multiplication), the addition and subtraction are available though,
%> - reshape a Faust.
%>
%> Primarily for convenience and test purposes, a Faust can be converted into
%> the corresponding full matrix using the function Faust.full.
%>
%> \warning using Faust.full is discouraged except for test purposes, as it
%> loses the main potential interests of the FAuST structure: compressed
%> memory storage and faster matrix-vector multiplication compared to its
%> equivalent full matrix representation.
%>
%> In this documentation, the expression 'full matrix' designates the Matlab array
%> Faust.full() obtained by the multiplication of the Faust factors.
%>
%> List of functions that are memory costly: Faust.full(), Faust.pinv(), Faust.mldivide().
%>
%> For more information about FAuST take a look at http://faust.inria.fr.
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
		%> Other easy ways to create a Faust is to call one of the following functions: matfaust.rand(), matfaust.dft() or matfaust.wht().
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; Faust(factors) creates a Faust from a list of factors (1D cell array).<br/><br>
		%> &nbsp;&nbsp;&nbsp; Faust(filepath) creates a Faust from a previously saved Faust file (filepath is a character array).<br/><br/>
		%>
		%>
		%> @param factors (1st arg.) the 1D cell array of factors to initialize the Faust with.
		%> <br/> The factors must respect the dimensions needed for the product to be defined (for i=1 to size(factors,2), size(factors{i},2) == size(factors{i+1},1)).
		%> <br/> The factors can be sparse or dense matrices.
		%> <br/> Passing a single matrix to the constructor instead of
		%> a cell array is equivalent to passing a singleton cell array.
		%> @param filepath (1st arg.) the file from which a Faust is created. It must be a character array.<br/>
		%>								The format is Matlab version 5 (.mat extension).<br/>
		%>								The file must have been saved before with Faust.save().
		%>
		%> @b Examples
		%> @code
		%>	import matfaust.*
		%>	factors = cell(1,5);
		%>	is_sparse = false;
		%>	for i=1:5
		%>		if(is_sparse) % odd index factors are sparse matrices
		%>			factors{i} = sprand(100, 100, 0.1);
		%>		else % even index gives a dense matrix
		%>			factors{i} = rand(100, 100);
		%>		end
		%>		is_sparse = ~ is_sparse;
		%>	end
		%>	% define a Faust with those factors
		%>	F = Faust(factors)
		%>
		%>
		%>	save(F, 'F.mat')
		%>	% define a Faust from file
		%>	H = Faust('F.mat')
		%>
		%>	Faust(rand(10,10)) % creating a Faust with only one factor
		%>
		%> @endcode
		%>
		%>
		%> <p>@b See @b also Faust.delete, Faust.save, matfaust.rand</p>
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
			% scale argument is hidden for user (deprecated) but it's still available
			if(nargin < 1 || nargin > 3)
				error([err_msg ' Number of arguments passed is zero or greater than three.'])
			elseif(iscell(varargin{1}))
				% init a Faust from factor list
				% check if the factors are real or complex, one complex factor implies a complex faust
				factors=varargin{1};
				isRealFlag = 1;
				if(length(factors) == 0)
					error([ err_msg ' Cannot create an empty Faust.'])
				end
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
			elseif(ismatrix(varargin{1}) && nargin == 1)
				c = cell(1, 1);
				c{1} = varargin{1};
				F = matfaust.Faust(c);
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
				% hack to raise an error in case of non-consistent handle and isReal arg.
				try
					n = numfactors(F);
				catch
					error('The Faust handle passed as first argument is not valid or not consistent with the value of isReal (2nd argument).')
				end
			else
				error(err_msg)
			end
		end


		%======================================================================
		%> @brief Deletes the Faust object F (destructor).
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b delete(F)
		%>
		%> @param F the Faust to delete.
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand([2, 5], [50, 100], .5)
		%>	delete(F)
		%>	F = matfaust.rand([2, 5], [50, 100], .5)
		%>	G = matfaust.rand([2, 5], [50, 100], .5)
		%>	clear % equivalent to delete(F);delete(G)
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust, clear (built-in)</p>
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
		%> @brief Plus
		%>
		%===
		%> This function overloads a Matlab built-in function.
		%>
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b plus(F,G) or F+G adds two Faust together, sizes must be compatible.<br/>
		%> &nbsp;&nbsp;&nbsp; @b plus(F,A) or F+A adds a Faust and a matrix A, sizes must be compatible.<br/>
		%> &nbsp;&nbsp;&nbsp; @b plus(F,s) or F+s, with s being a scalar, such that full(F+s) == full(F)+s.<br/>
		%>
		%> @param F (first arg.) The Faust object.
		%> @param G, A, s,… (2nd to n-th args) The variables to add to F; Fausts, matrices or scalars.
		%>
		%> @retval S: the sum as a Faust object.
		%>
		%> <p>@b See @b also Faust.minus
		%======================================================================
		function F = plus(varargin)
			import matfaust.Faust
			F = varargin{1};
			for i=2:nargin
				G = varargin{i};
				if(isa(G,'matfaust.Faust'))
					if(any(size(G) ~= size(F)))
						error('Dimensions must agree.')
					end
					C = [F,G];
					if(~ C.isReal)
						Fid = matfaust.eye(size(C,2)/2, 'complex');
					else
						Fid = matfaust.eye(size(C,2)/2);
					end
					F = C*[Fid;Fid];
				elseif(isscalar(G))
					if(~ isreal(G) && F.isReal)
						F = complex(F);
					elseif(~ isreal(F) && isreal(G))
						G = complex(G);
					end
					if(size(F,1) <= size(F,2))
						F = F+(Faust({speye(size(F,1),size(F,2)), ones(size(F,2), 1)*G, ones(1, size(F,2))}));
					else
						F = F+(Faust({speye(size(F,2),size(F,1)), ones(size(F,1), 1)*G, ones(1, size(F,1))})).';
					end
				elseif(ismatrix(G))
					F = plus(F,Faust(G));
				else
					error('Cannot add a Faust to something that is not a Faust, a matrix/array, or a scalar.')
				end
			end
		end

		%======================================================================
		%> @brief Minus
		%>
		%===
		%> This function overloads a Matlab built-in function.
		%>
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b minus(F,G) or F-G subtracts the Faust G from F, sizes must be compatible.<br/>
		%> &nbsp;&nbsp;&nbsp; @b minus(F,A) or F-A subtracts a matrix A from F, sizes must be compatible.<br/>
		%> &nbsp;&nbsp;&nbsp; @b minus(F,s) or F-s subtracts a scalar s from F, such that full(F-s) == full(F)-s.<br/>
		%>
		%> @param F (first arg.) The Faust object.
		%> @param G, A, s,… (2nd to n-th args) The variables to subtract from F; Fausts, matrices or scalars.
		%>
		%> @retval M: the difference as a Faust object.
		%>
		%> <p>@b See @b also Faust.plus
		%======================================================================
		function M = minus(varargin)
			M = varargin{1};
			for i=2:nargin
				varargin{i} = varargin{i}*-1;
			end
			M  = plus(M, varargin{2:end});
		end

		%======================================================================
		%>  /   Slash or right Faust divide.
		%===
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; G = F/s is the division of the Faust F by the scalar s, such that full(F)/s == full(F/s)<br/>
		%> &nbsp;&nbsp;&nbsp; G = MRDIVIDE(F,s) is called for the syntax 'F / s' when s is a scalar.
		%>
		%> @param F The Faust object to divide.
		%> @param s The scalar to divide F with.
		%>
		%> @retval G: the division result as a Faust object.
		%>
		%> <p>@b See @b also Faust.mtimes
		%======================================================================
		function G = mrdivide(F,s)
			if(~ isscalar(s))
				error('Unsupported operand type(s) for /: a Faust can only be divided by a scalar.')
			end
			G = mtimes(F,1/s);
		end

		%======================================================================
		%> @brief Multiplies the Faust F by A which is a full matrix, a Faust object or a scalar.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b The @b primary @b goal of this function is to implement “fast” multiplication by a
		%> Faust, which is the operation performed when A is a full matrix.<br/>
		%> In the best case, F*A is rcg(F) times faster than performing the
		%> equivalent full(F)*A.<br/>
		%>
		%> @b Other @b use @b cases are available for this function:
		%> - If A is a Faust, no actual multiplication is performed, instead a new Faust is built to implement the multiplication. This Faust verifies that:
		%> @code
		%> full(F*A) == full(F)*full(A)
		%> @endcode
		%> @note you could have an elementwise non-significant absolute difference between the two members (not more than eps(1.0)).
		%>
		%> - If A is a scalar, F*A is also a Faust such that:
		%> @code
		%> factors(F*A,1) == factors(F,1)*A
		%> @endcode
		%>
		%> - Any Faust F (real or complex) can be multiplied by any Faust/array/scalar (real or complex) of appropriate dimensions.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b G = mtimes(F, A)<br/>
		%> &nbsp;&nbsp;&nbsp; @b G = F*A with A and G being full matrices. A could also be a sparse matrix but note that often the Faust-sparse matrix multiplication is slower than performing F*full(A) multiplication. In some cases though, it stays quicker: moreover when the Faust is composed of a small number of factors.<br/>
		%> &nbsp;&nbsp;&nbsp; @b G = F*s with s a scalar and G a Faust.<br/>
		%> &nbsp;&nbsp;&nbsp; @b G = s*F with s a scalar and G a Faust.<br/>
		%> &nbsp;&nbsp;&nbsp; @b G = s*F' with s a scalar multiplying the conjugate-transpose of F.<br/>
		%> &nbsp;&nbsp;&nbsp; @b G = F'*s with s a scalar multiplying the conjugate-transpose of F.<br/>
		%>
		%>
		%>
		%> @param F the Faust object.
		%> @param A either a full matrix, a Faust object or a scalar.
		%>
		%> @retval G the multiplication result:
		%> - When A is a full matrix, G = F*A is also a full matrix.
		%> - When A is a Faust or a scalar, G = F*A is itself a Faust.
		%> - When either F or A is complex, G=F*A is also complex.
		%>
		%>
		%>
		%> @b Example
		%> @code
		%>   F = matfaust.rand([2, 5], [50, 100], .5)
		%>   A = rand(size(F,2), 50)
		%>   G = F*A
		%> % is equivalent to G = mtimes(F, A)
		%>   G = matfaust.rand(5,size(F, 2))
		%>   H = F*G
		%> % H is a Faust because F and G are
		%>   Fs = F*2
		%>   sF = 2*F
		%> % sF == Fs, i.e. the Faust F times 2.
		%> @endcode
		%>
		%> @b Errors
		%>
		%> - F is real but A is a complex scalar.
		%>
		%> @code
		%>>> isreal(F)
		%>
		%>ans =
		%>
		%>  logical
		%>
		%>     1
		%>
		%>>>F*i
		%>Error using mexFaustReal
		%>You cannot multiply a real Faust by a complex scalar (not yet implemented).
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust, Faust.rcg, Faust.ctranspose, Faust.complex
		%>
		%======================================================================
		function G = mtimes(F,A)
			%% MTIMES * Faust Multiplication (overloaded Matlab built-in function).
			%
			% G=mtimes(F,A) is called for syntax 'G=F*A', when F is a Faust matrix and A a full
			% storage matrix, G is also a full matrix storage.
			%
			% See also mtimes_trans
			if(isa(A, 'matfaust.Faust'))
				if(isscalar(F))
					G = mtimes_trans(A, F, 0);
				elseif(ismatrix(F))
					G = mtimes_trans(A', F', 0)';
				end
			else
				G = mtimes_trans(F, A, 0);
			end
		end


		%======================================================================
		%> @brief Multiplies the Faust F by A.' or A which is a dense matrix or a Faust object.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @warning if A is a matrix the function costs numfactors(F) matrix multiplications.
		%> In that case its use is discouraged except for test purpose. However if A is a Faust object,
		%> it costs the same that a Faust initialization with a number of factors equal to
		%> F.numfactors()+A.numfactors() (like you can do directly with Faust.Faust).
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
			% TODO: take trans into account when mul F to a scal or a Faust
			% it's not a serious issue because mtimes_trans() shouln't be called by final user
			% (func's doc is filtered out in doxydoc)
			%if(issparse(A))
			%	error('Faust multiplication to a sparse matrix isn''t supported.')
			%elseif(isa(A,'matfaust.Faust'))
			if(isa(A,'matfaust.Faust'))
				if (F.isReal)
					if(isreal(A))
						C = matfaust.Faust(F, mexFaustReal('mul_faust', F.matrix.objectHandle, A.matrix.objectHandle));
					else
						C = mtimes_trans(complex(F), A, trans);
					end
				else
					if(A.isReal)
						A = complex(A);
					end
					C = matfaust.Faust(F, mexFaustCplx('mul_faust', F.matrix.objectHandle, A.matrix.objectHandle));
				end
			elseif(isscalar(A))
				if (F.isReal)
					if(isreal(A))
						C = matfaust.Faust(F, mexFaustReal('mul_scalar', F.matrix.objectHandle, A));
					else
						C = mtimes_trans(complex(F), A, trans);
					end
				else
					C = matfaust.Faust(F, mexFaustCplx('mul_scalar', F.matrix.objectHandle, A));
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
		%> @brief The full matrix implemented by F.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b A = full(F)
		%>
		%> @warning using this function is discouraged except for test purposes,
		%> as it loses the main potential interests of the FAuST structure:
		%> compressed memory storage and faster matrix-vector multiplication compared
		%> to its equivalent dense matrix representation.
		%>
		%>
		%> @retval A the Matlab native matrix such that A*x == F*x
		%> for any vector x.
		%>
		%> @warning Running the example below is likely to raise a memory error or freeze
		%> your computer for a certain amount of time.
		%>
		%> @b Example
		%> @code
		%>   % in a matlab terminal
		%>   >> F = matfaust.rand(5, 10^6, 10^-5, 'sparse')
		%>   Faust size 1000000x1000000, density 5e-05, nnz_sum 49999995, 5 factor(s):
		%>   - FACTOR 0 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   - FACTOR 1 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   - FACTOR 2 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   - FACTOR 3 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   - FACTOR 4 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999
		%>   >> % an attempt to convert F to a full matrix is most likely to raise a memory error
		%>   >> % the sparse format is the only way to handle such a large Faust
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
		%> @b @note if F is a real Faust, all factors are real.
		%> If F is a complex Faust, at least one factor is complex.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b bool = isreal(F)
		%>
		%> @retval bool 1 if F is a real Faust, 0 if it's a complex Faust.
		%>
		%> <p/>@b See @b also Faust.complex
		%======================================================================
		function bool = isreal(F)
			%% ISREAL True for real scalar Faust (overloaded Matlab built-in function).
			%
			% isreal(F) returns 1 if Faust F does not have an imaginary part
			% and 0 otherwise.

			bool=F.isReal;

		end

		%======================================================================
		%> @brief The transpose of F.
		%>
		%> This function overloads a Matlab built-in function/operator.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp;  @b F_trans = @b transpose(F)<br/>
		%> &nbsp;&nbsp;&nbsp;  @b F_trans = @b F.'
		%>
		%> @param F the Faust object.
		%>
		%> @retval F_trans a Faust object implementing the transpose of full(F) such that:
		%> <code>full(F_trans) == full(F).' == transpose(full(F)).</code>
		%>
		%>
		%> @b Example
		%> @code
		%>   F_trans = F.'
		%> % is equivalent to
		%>   F_trans = transpose(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.conj, Faust.ctranspose
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
		%> @brief The conjugate transpose of F.
		%>
		%> This function overloads a Matlab built-in function/operator.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp;  @b F_ctrans = @b ctranspose(F)<br/>
		%> &nbsp;&nbsp;&nbsp;  @b F_ctrans = @b F'
		%>
		%> @param F the Faust object.
		%>
		%> @retval F_ctrans a Faust object implementing the conjugate transpose of full(F), such that:<br/>
		%> <code>full(F_ctrans) == full(F)' == ctranspose(full(F))</code>
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex')
		%>	F_ctrans = F'
		%>	F_ctrans2 = ctranspose(F)
		%>	% F_ctrans == F_ctrans2
		%>	F_ctrans2 = transpose(F)
		%>	F_ctrans2 = conj(F_ctrans2)
		%>	% F_ctrans == F_ctrans2
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.transpose, Faust.conj, Faust.complex
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
		%> @brief The complex conjugate of F.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b F_conj = conj(F)
		%>
		%> @param F the Faust object.
		%>
		%> @retval F_conj a Faust object implementing the conjugate of full(F), such that:<br/>
		%> <code>full(F_conj) == conj(full(F))</code>
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex')
		%>	F_conj = conj(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.transpose, Faust.ctranspose, Faust.complex
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
		%> @brief Returns a Faust optimized by removing useless zero rows and columns as many times as needed.
		%>
		%> @param F: the Faust to optimize.
		%> @param 'nnz_tres', int: (optional) the treshold of number of nonzeros under what the
		%>            rows/columns are removed.
		%> @param 'only_forward', bool: (optional) True for applying only the forward passes of removal.
		%> @param 'npasses', int: (optional) the number of passes to run, by default it goes until the
		%>            optimal Faust is obtained.
		%>
		%> @retval G The optimized Faust.
		%======================================================================
		function G = pruneout(F, varargin)
			nnz_tres = 0;
			only_forward = false;
			npasses = -1;
			argc = length(varargin);
			if(argc > 0)
				for i=1:argc
					switch(varargin{i})
						case 'only_forward'
							if(argc == i || ~ islogical(varargin{i+1}))
								error('only_forward keyword argument is not followed by a logical')
							else
								only_forward = varargin{i+1};
							end
						case 'nnz_tres'
							if(argc == i || ~ isnumeric(varargin{i+1}) || varargin{i+1}-floor(varargin{i+1}) > 0 || varargin{i+1} <0 )
								error('nnz_tres keyword arg. is not followed by a positive integer')
							end
							nnz_tres = varargin{i+1};
						case 'npasses'
							if(argc == i ||  ~ isnumeric(varargin{i+1}) || varargin{i+1}-floor(varargin{i+1}) > 0)
								error('npasses keyword arg. is not followed by an integer')
							else
								npasses = varargin{i+1};
							end
						otherwise
							if(isstr(varargin{i}))
								error([ varargin{i} ' unrecognized argument'])
							end
					end
				end
			end
			if(F.isReal)
				G = matfaust.Faust(F, mexFaustReal('pruneout', F.matrix.objectHandle, nnz_tres, npasses, only_forward));
			else
				G = matfaust.Faust(F, mexFaustCplx('pruneout', F.matrix.objectHandle, nnz_tres, npasses, only_forward));
			end
		end

		%=====================================================================
		%> @brief Optimizes a Faust by changing the storage format of each factor in
		%> order to optimize the memory size.
		%>
		%=====================================================================
		function OF = optimize_storage(F)
			if(F.isReal)
				OF = matfaust.Faust(F, mexFaustReal('optimize_storage', F.matrix.objectHandle, false));
			else % cplx Faust
				OF = matfaust.Faust(F, mexFaustCplx('optimize_storage', F.matrix.objectHandle, false));
			end
		end

		%=====================================================================
		%> @brief Returns a Faust optimized with pruneout, optimize_storage and configured with the quickest method available to compute a Faust-matrix product.
		%>
		%> NOTE: this function launches a small benchmark on the fly. Basically, the methods
		%> available differ by the order used to compute the matrix chain
		%> multiplication or by the use (or unuse) of threads for the calculation of intermediary
		%> matrix products of the Faust.
		%>
		%=====================================================================
		function OF = optimize(F, varargin)
			len_args = length(varargin);
			transp = false; % default value
			if(len_args > 0)
				if(strcmp(varargin{1}, {'transpose', 'transp'}) &&  len_args > 1 && islogical(varargin{2}))
					transp = varargin{2};
				else
					error('invalid key or value arguments')
				end
			end
			if(F.isReal)
				OF = matfaust.Faust(F, mexFaustReal('optimize', F.matrix.objectHandle, transp));
			else % cplx Faust
				OF = matfaust.Faust(F, mexFaustCplx('optimize', F.matrix.objectHandle, transp));
			end
		end

		function OF = optimize_mul(F, varargin)
			len_args = length(varargin);
			transp = false; % default value
			if(len_args > 0)
				if(strcmp(varargin{1}, {'transpose', 'transp'}) &&  len_args > 1 && islogical(varargin{2}))
					transp = varargin{2};
				else
					error('invalid key or value arguments')
				end
			end
			if(F.isReal)
				mexFaustReal('optimize_mul', F.matrix.objectHandle, transp);
			else % cplx Faust
				mexFaustCplx('optimize_mul', F.matrix.objectHandle, transp);
			end
		end

		%======================================================================
		%> @brief The size of F.
		%>
		%> The size is a pair of numbers: the number of rows and the number of columns
		%> of full(F).
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; [@b nrows, @b ncols] = @b size(F)<br/>
		%> &nbsp;&nbsp;&nbsp; @b n = @b size(F,dim) with n being the size of the dim-th dimension of F.<br/>
		%> &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; In other words n == nrows if dim == 1, n == ncols if dim == 2.
		%>
		%> @param F the Faust object.
		%> @param dim (optional) the index of the dimension to get the size of.
		%>
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex')
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
		%> @brief The number of elements in F.
		%>
		%> It's equivalent to <code>prod(size(F))</code>.
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
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex')
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
		%> @brief The last index when slicing or indexing a Faust.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%>
	 	%> @b Example
		%> @code
		%>	% in a matlab terminal
		%>	>> F = matfaust.rand([4, 5], 3, .9);
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
		%> <p/>@b See @b also Faust.subsref, Faust.size
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
		%> @brief Returns the i-th factor or a range of factors of F.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b factor = factors(F, i) returns the i-th factor of F.<br/>
		%> &nbsp;&nbsp;&nbsp; @b factor = factors(F, i:j) returns a new Faust formed of the F's factors from the i-th to the j-th included.
		%>
		%> @param F the Faust object.
		%> @param varargin the factor indices.
		%>
		%> @retval factors a matrix copy of the i-th factor if i is a single index or a new Faust composed of i-th to the j-th factors of F. The factors copies keep the storage organization of the source matrix (full or sparse).
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex');
		%>	f1 = factors(F, 1);
		%>	G = factors(F, 4:5); % a new Faust composed of the two last factors of F
		%> @endcode
		%> <p>@b See @b also Faust.numfactors
		%=====================================================================
		function factors = factors(F, varargin)
			factors = cell(1, size(varargin{1},2));
			for j=1:length(factors)
				i = varargin{1};
				if(j < length(factors) && i(j+1) - i(j) ~= 1)
					error('Indices must be contiguous.')
				end
				i = i(j);
				if (~isa(i,'double'))
					error('factors second argument (indice) must either be real positive integers or logicals.');
				end

				if (floor(i) ~= i)
					error('factors second argument (indice) must either be real positive integers or logicals.');
				end

				if (F.isReal)
					factors{j} = mexFaustReal('get_fact',F.matrix.objectHandle,i);
				else
					factors{j} = mexFaustCplx('get_fact',F.matrix.objectHandle,i);
				end
			end
			if(length(factors) > 1)
				factors = matfaust.Faust(factors);
			else
				factors = factors{j};
			end
		end
		%================================================================
		%> Returns the left hand side factors of F from index 1 to i included.
		%===
		%>
		%> <p> @b See @b also Faust.factors, Faust.right
		%================================================================
		function lfactors = left(F, i)
			i = check_factor_idx(F, i);
			lfactors = factors(F, 1:i);
		end

		%================================================================
		%> Returns the right hand side factors of F from index i to end.
		%===
		%>
		%> <p> @b See @b also Faust.factors, Faust.left
		%================================================================
		function rfactors = right(F, i)
			i = check_factor_idx(F, i);
			rfactors = factors(F, i:numfactors(F));
		end

		%=====================================================================
		%> @brief DEPRECATED (Use Faust.factors) Returns the i-th factor of F.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b factor = get_factor_nonopt(F, i)
		%>
		%> @param F the Faust object.
		%> @param i the factor index.
		%>
		%> @retval factor the i-th factor as a dense matrix.
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex')
		%>	f1 = get_factor_nonopt(F, 1)
		%> @endcode
		%> <p>@b See @b also Faust.numfactors
		%=====================================================================
		function factor = get_factor_nonopt(F, i)
			%% GET_FACT Ith factor of the Faust.
			%
			% A=get_factor_nonopt(F,i) return the i factor A of the Faust F as a full storage matrix.
			%
			% Example of use :
			% A=get_factor_nonopt(F,1) returns the 1st factor of the Faust F.
			% A=get_factor_nonopt(F,4) returns the 4th factor of the Faust F.
			%
			% See also numfactors.

			if (~isa(i,'double'))
				error('get_fact second argument (indice) must either be real positive integers or logicals.');
			end

			if (floor(i) ~= i)
				error('get_fact second argument (indice) must either be real positive integers or logicals.');
			end

			if (F.isReal)
				factor = mexFaustReal('get_fact_nonopt',F.matrix.objectHandle,i);
			else
				factor = mexFaustCplx('get_fact_nonopt',F.matrix.objectHandle,i);
			end

		end

		%==========================================================================================
		%> @brief The number of factors of F.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b A = numfactors(F)
		%>
		%> @param F the Faust object.
		%>
		%> @retval num_factors the number of factors.
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex')
		%>	nf = numfactors(F)
		%> @endcode
		%>
		%> <p>@b See @b also Faust.factors.
		%==========================================================================================
		function num_factors = numfactors(F)
			%% GET_NB_FACTOR Number of factor of the Faust.
			%
			% num_factors = numfactors(F) return the number of factor of the
			% Faust F.
			%
			% See also factors.
			if (F.isReal)
				num_factors = mexFaustReal('get_nb_factor', F.matrix.objectHandle);
			else
				num_factors = mexFaustCplx('get_nb_factor', F.matrix.objectHandle);
			end
		end

		%==========================================================================================
		%> @brief Saves the Faust F into a file.
		%>
		%> The file is saved in Matlab format version 5 (.mat extension).
		%>
		%> @b @note storing F should typically use rcg(F) times less disk space than storing full(F).
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
		%>	import matfaust.*
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex')
		%>	save(F, 'F.mat')
		%>	G = Faust('F.mat')
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
		%> @brief Subscripted reference of a Faust.
		%>
		%> The function returns a Faust representing a submatrix of full(F) or a
		%> scalar element if that Faust can be reduced to a single element.
		%>
		%>
		%> This function overloads a Matlab built-in.
		%>
		%> @warning
		%> - It is not advised to use this function as an element accessor
		%>        (e.g. F(1,1)) because such a use induces to convert the Faust to its
		%>        dense matrix representation and that is a very expensive computation if used
		%>        repetitively.
		%> - Subindexing a Faust which would create an empty Faust will raise an error.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b G = F(I,J) the Faust representing the sub-matrix of full(F) defined by the subscript vectors I and J (see examples below).
		%>
		%> @param F the Faust object.
		%> @param S the structure defining the Faust to extract like a submatrix (it's not supposed to be used directly ; usage and examples show how subscript should be used).
		%>
		%> @retval The Faust object requested or just the corresponding scalar if that Faust has
		%>             a size equal to [1,1].
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand([2, 5], [50, 100], .5)
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
		%>	F(1:2:end,:) % every row of odd index
		%>	F(end:-2:1,:) %  starts from the last row and goes backward to take one in two rows until the first one (reversing row order of F)
		%>	F([1,18,2],:) % takes in this order the rows 1, 18 and 2.
		%>	F(:,[1,18,2]) % takes in this order the columns 1, 18 and 2
		%>	F([1,18,2], [1,2]) % takes the rows 1, 18 and 2 but keeps only columns 1 and 2 in these rows.
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
				error(' subsref invalid slicing must have 2 index since F is a 2D-array');
			end

			end_ids = zeros(1,2);
			start_ids = zeros(1,2);
			indexing_by_slice = [ true, true ];
			ind_lists = cell(1,2);
			ROW=1;
			COL=2;

			for i=ROW:COL
				ind_list=S.subs{i};
				if ischar(ind_list)
					start_ids(i) = 1;
					end_ids(i) = size(F,i);
					%indexing_by_slice(i) = true
					ind_list = 1:size(F,i); % needed if other dim is not a slice
				else
					if(any(ind_list < 1))
						error(' Subscript indices must be integers >= 1.')
					elseif(any(ind_list > size(F,i)))
						error(' Index exceeds Faust dimensions.')
					elseif(size(ind_list,2) == 0)
						error(' Cannot create empty Faust')
					end
					% check if indices in range are contiguous and not negative step
					sl_sz = size(ind_list,2);
					if(sl_sz >= 2)
						% TODO: couldn't be without a loop by verifiying two shifted views of array ?
						for j=2:sl_sz
							d = ind_list(j)-ind_list(j-1);
							if(abs(d) > 1 || d < 0)
								indexing_by_slice(i) = false;
								break
							end
						end
					end
					if(indexing_by_slice(i))
						start_ids(i) = ind_list(1);
						end_ids(i) = ind_list(end);
					end
				end
				ind_lists{i} = ind_list;
			end

			if(F.isReal)
				if(indexing_by_slice)
					submatrix = matfaust.Faust(F, mexFaustReal('subsref', F.matrix.objectHandle, start_ids(ROW), end_ids(ROW), start_ids(COL), end_ids(COL)));
				else
					% -1 is for converting to 0-base index used in C world
					submatrix =  matfaust.Faust(F, mexFaustReal('subsref_byvec', F.matrix.objectHandle, uint64(ind_lists{1}-1), uint64(ind_lists{2}-1)));
				end
			else
				if(indexing_by_slice)
					submatrix = matfaust.Faust(F, mexFaustCplx('subsref', F.matrix.objectHandle, start_ids(ROW), end_ids(ROW), start_ids(COL), end_ids(COL)));
				else
					submatrix = matfaust.Faust(F, mexFaustCplx('subsref_byvec', F.matrix.objectHandle, uint64(ind_lists{1}-1), uint64(ind_lists{2}-1)));
				end
			end

			if(size(submatrix) == [1,1])
				submatrix = full(submatrix);
				submatrix = submatrix(1,1);
			end
		end

		%======================================================================
		%> @brief @warning this function is not implemented because a Faust object is immutable.
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
		%>	>> F = matfaust.rand([1, 2], [50, 100], .5)
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
		%> <p>@b See @b also Faust.nnz_sum, Faust.density, Faust.size, Faust.factors, Faust.numfactors
		%>
		%>
		%======================================================================
		function disp(F)
			%% DISP shows the characteristic of the Faust (overloaded Matlab built-in function)
			%
			%
			% This function shows the size of the Faust,
			%  its number of factor, its RCG, …
			%
			if (F.isReal)
				mexFaustReal('disp',F.matrix.objectHandle);
			else
				mexFaustCplx('disp',F.matrix.objectHandle);
			end


		end

		%======================================================================
		%> @brief The matrix norm of F.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> Several types of norm are available: 1-norm, 2-norm, inf-norm and Frobenius norm.
		%>
		%> The norm of F is equal to the norm of full(F).
		%>
		%> @warning the computation time can be expected to be of order n*min(size(F)) with n the time for multipliying F by a vector. Nevertheless, the implementation ensures that memory usage remains controlled by avoiding to explicitly compute full(F) (at least for 2-norm).
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b n = @b norm(F, 2) the 2-norm or maximum singular value of F: approximately norm(full(F),2) == max(svd(full(F))).<br/><br/>
		%> &nbsp;&nbsp;&nbsp; @b n = @b norm(F, 2, 'threshold', 0.001, 'max_num_its', 1000).<br/><br/>
		%> &nbsp;&nbsp;&nbsp; @b n = @b norm(F) the same as norm(F, 2).<br/><br>
		%> &nbsp;&nbsp;&nbsp; @b n = @b norm(F, 1) the 1-norm of F: norm(full(F), 1) == max(sum(abs(full(F))))  <br/><br/>
		%> &nbsp;&nbsp;&nbsp; @b n = @b norm(F, inf) the inf-norm of F: norm(full(F), inf) == max(sum(abs(full(F)'))) <br/><br/>
		%> &nbsp;&nbsp;&nbsp; @b n = @b @b norm(F, @b 'fro') the Frobenius norm of F: norm(full(F), 'fro').<br/><br/>
		%>
		%> @param F the Faust object.
		%> @param p (optional) the norm order or type. Respectively 1, 2 or inf for the 1-norm, 2-norm and inf-norm or 'fro' for the Frobenius norm (by default the 2-norm is computed).
		%> @param threshold (optional) power iteration algorithm threshold (default to .001). Used only for norm(2). It's passed in a key-value pair fashion: 'threshold', .001
		%> @param max_num_its (optional) maximum number of iterations for power iteration algorithm. Used only for norm(2). It's passed in a key-value pair fashion: 'max_num_its', 1000.
		%>
		%>
		%> @retval n the norm (real).
		%>
		%>
		%> @b Example
		%> @code
		%> % in a matlab terminal
		%> >> F = matfaust.rand([1, 2], [50, 100], .5)
		%> >> norm(F)
		%> ans =
		%> 7.0151
		%> >> norm(F,2)
		%> ans =
		%> 7.0151
		%> >> norm(F, 'fro')
		%> ans =
		%> 30.1014
		%> >> norm(F, inf)
		%> ans =
		%> 25.9277
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
			% @warning norm(F,P) is only supported when P equals 1, 2 or
			% 'fro'.

			nb_input = length(varargin);
			if (nb_input > 5)
				error('Too many input arguments');
			end

			ord = 2;
			args = {ord};
			ord2_valid_param = false;
			if nb_input >= 1
				if(varargin{1} == 'fro')
					if (F.isReal)
						n = mexFaustReal('normfro',F.matrix.objectHandle);
					else
						n = mexFaustCplx('normfro',F.matrix.objectHandle);
					end
					return
				end
				if varargin{1} ~= 2 && varargin{1} ~= 1 && varargin{1} ~= inf
					error('only 1, 2, inf or Frobenius norms are supported for the Faust');
				end
				ord = varargin{1};
				args = {ord};
				extra_opts = cell(2);
				if(ord == 2)
					keys = {'threshold', 'max_num_its'};
					i = 2;
					while i <= nb_input
						for j=1:length(keys)
							if(strcmp(varargin{i}, keys{j}))
								if(nb_input>i && isnumeric(varargin{i+1}) && isreal(varargin{i+1}))
									eval 'extra_opts{j} = varargin{i+1};';
									ord2_valid_param = true;
									i = i + 1;
								else
									error(['Parameter value (argument index ' num2str(i+2) ') not valid for parameter key ' keys{j}]);
								end
							end
						end
						if(~ ord2_valid_param)
							error(['invalid argument ' num2str(i+1)])
						end
						i = i + 1;
						ord2_valid_param = false;
					end
				end
			end
			if(ord == 2 && nb_input > 2)
				args = [ args, extra_opts{:} ];
			end
			if (F.isReal)
				n = mexFaustReal('norm',F.matrix.objectHandle, args{:});
			else
				n = mexFaustCplx('norm',F.matrix.objectHandle, args{:});
			end

		end

		%===============================================================================
		%> Returns the normalized F.
		%===
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b normalize(F,'norm')  normalizes columns with vector 2-norm.<br/>
		%> &nbsp;&nbsp;&nbsp; @b normalize(F,'norm',p) normalizes columns with vector p-norm (see Faust.norm).<br/>
		%> &nbsp;&nbsp;&nbsp; @b normalize(F,DIM,'norm') normalizes rows (DIM=1) or columns (DIM=2, default) with vector 2-norm.<br/>
		%> &nbsp;&nbsp;&nbsp; @b normalize(F,DIM,'norm',p) normalizes rows (DIM=1) or columns (DIM=2, default) with vector p-norm.<br/>
		%>
		%> <p>@b See @b also Faust.norm
		%===============================================================================
		function NF = normalize(F, varargin)
			nargs = length(varargin);
			dim = 2;
			for i=1:nargs;
				if(strcmp(varargin{i}, 'norm'))
					if(nargs > i)
						ord = varargin{i+1};
					else
						ord = 2;
					%	error('norm argument passed without a norm order in argument just after')
					end
					break
				elseif(isnumeric(varargin{i}) && ~ exist('ord'))
					dim = varargin{i};
					if(dim ~= 1 && dim ~= 2);
						error('dimension is not valid')
					end
				else
					error('invalid argument.')
				end
			end
			if(~ exist('ord'))
				ord = 2;
			end
			switch(ord)
				case {1,2}
					mex_ord = ord;
				case Inf
					mex_ord = -1;
				case 'fro'
					mex_ord = -2;
				otherwise
					error('Invalid type of norm.')
			end
			if(dim == 1)
				F = F.';
				if(mex_ord == -1)
					mex_ord = 1;
				elseif(mex_ord == 1)
					mex_ord = -1;
				end
			end
			if(F.isReal)
				NF = matfaust.Faust(F, mexFaustReal('normalize', F.matrix.objectHandle, mex_ord));
			else
				NF = matfaust.Faust(F, mexFaustCplx('normalize', F.matrix.objectHandle, mex_ord));
			end
			if(dim == 1)
				NF = NF.';
			end
		end


		%==========================================================================================
		%> @brief The total number of non-zero elements in the factors of F.
		%>
		%> The function sums together the number of non-zero elements of
		%> each factor and returns the result. Note that for efficiency the sum is
		%>                        computed at Faust creation time and kept in cache.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b nz = nnz_sum(F)
		%> @param F the Faust object.
		%>
		%> @retval nz the number of non-zeros.
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
		%> @brief The density of F such that nnz_sum(F) == density(F)*numel(F).
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b dens = density(F)
		%>
		%>
		%> @note A value of density below one indicates potential memory savings compared
		%> to storing the corresponding dense matrix full(F), as well as potentially
		%> faster matrix-vector multiplication when applying F*x instead of full(F)*x.
		%>
		%> @note A density above one is possible but prevents any saving.
		%>
		%>
		%> @param F the Faust object.
		%>
		%> @retval dens density of F.
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand([2, 5], [50, 100], .5)
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
		%> @brief The Relative Complexity Gain of F.
		%>
		%>
		%> RCG is the theoretical gain brought by the Faust representation relatively to its dense
		%> matrix equivalent. <br/>The higher is the RCG, the more computational
		%> savings will be made.
		%> This gain applies both for storage space and computation time.
		%>
		%> @b @note rcg(F) == 1/density(F)
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b gain = rcg(F)
		%>
		%> @param F	the Faust object.
		%>
		%>
		%> @retval gain the RCG value (real).
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

		%======================================================================
		%> @brief Concatenates F with n Faust objects or full/sparse matrices.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> The resulting Faust C = cat(DIM,F,G,…) verifies that:
		%> @code
		%> full(C) == full(cat(DIM, full(F), full(G), …))
		%> @endcode
		%> @note you could have an elementwise non-significant absolute difference between the two members (not more than eps(1.0)).
		%>
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b C=CAT(DIM,F,G) concatenates the Faust F and G, which is a Faust or a matrix, along the dimension DIM. The result is the Faust C. <br/>
		%> &nbsp;&nbsp;&nbsp; @b CAT(2,F,G) is the same as [F,G].<br/>
		%> &nbsp;&nbsp;&nbsp; @b CAT(1,F,G) is the same as [F;G].<br/>
		%> &nbsp;&nbsp;&nbsp; @b CAT(@b DIM, @b F, @b G,…) concatenates an arbitrary number of Fausts and matrices along the DIM-th dimension.<br/>
		%>
		%>
		%> @param DIM (1st arg.) the dimension along which to concatenate. DIM==1 means vertical concatenation, DIM==2 means horizontal concatenation.
		%> @param F (2nd arg.) the first Faust object.
		%> @param G,… (3rd to n+3 args) the Faust objects or matrices to concatenate to F. The matrices can be sparse or full.
		%>
		%> @retval C the concatenation result as a new Faust.
		%>
		%>
		%> @b Example
		%> @code
		%>% in a matlab terminal
		%>>> F = matfaust.rand(5,50);
		%>>> G = matfaust.rand(6,50);
		%>>> [F;G] % equiv. to cat(1,F,G)
		%>
		%>ans =
		%>
		%>Faust size 100x50, density 0.5592, nnz_sum 2796, 7 factor(s):
		%>- FACTOR 0 (real) SPARSE, size 100x100, density 0.0481, nnz 481
		%>- FACTOR 1 (real) SPARSE, size 100x100, density 0.0471, nnz 471
		%>- FACTOR 2 (real) SPARSE, size 100x100, density 0.0472, nnz 472
		%>- FACTOR 3 (real) SPARSE, size 100x100, density 0.0508, nnz 508
		%>- FACTOR 4 (real) SPARSE, size 100x100, density 0.0476, nnz 476
		%>- FACTOR 5 (real) SPARSE, size 100x100, density 0.0288, nnz 288
		%>- FACTOR 6 (real) SPARSE, size 100x50, density 0.02, nnz 100
		%>>>[F,G] % equiv. to cat(2,F,G)
		%>
		%>ans =
		%>
		%>Faust size 50x100, density 0.5592, nnz_sum 2796, 7 factor(s):
		%>- FACTOR 0 (real) SPARSE, size 50x100, density 0.02, nnz 100
		%>- FACTOR 1 (real) SPARSE, size 100x100, density 0.0286, nnz 286
		%>- FACTOR 2 (real) SPARSE, size 100x100, density 0.0477, nnz 477
		%>- FACTOR 3 (real) SPARSE, size 100x100, density 0.0476, nnz 476
		%>- FACTOR 4 (real) SPARSE, size 100x100, density 0.0511, nnz 511
		%>- FACTOR 5 (real) SPARSE, size 100x100, density 0.0466, nnz 466
		%>- FACTOR 6 (real) SPARSE, size 100x100, density 0.048, nnz 480
		%>>>[F;rand(100,50)] % vertical concatenation with auto-conversion of the random matrix to a Faust
		%>
		%>ans =
		%>
		%>Faust size 150x50, density 0.865733, nnz_sum 6493, 6 factor(s):
		%>- FACTOR 0 (real) SPARSE, size 150x100, density 0.349667, nnz 5245
		%>- FACTOR 1 (real) SPARSE, size 100x100, density 0.0289, nnz 289
		%>- FACTOR 2 (real) SPARSE, size 100x100, density 0.0285, nnz 285
		%>- FACTOR 3 (real) SPARSE, size 100x100, density 0.0282, nnz 282
		%>- FACTOR 4 (real) SPARSE, size 100x100, density 0.0292, nnz 292
		%>- FACTOR 5 (real) SPARSE, size 100x50, density 0.02, nnz 100
		%>[F;sprand(100,50,.2)] % vertical concatenation with auto-conversion of the random sparse matrix to a Faust
		%>
		%>ans =
		%>
		%>Faust size 150x50, density 0.3204, nnz_sum 2403, 6 factor(s):
		%>- FACTOR 0 (real) SPARSE, size 150x100, density 0.077, nnz 1155
		%>- FACTOR 1 (real) SPARSE, size 100x100, density 0.0289, nnz 289
		%>- FACTOR 2 (real) SPARSE, size 100x100, density 0.0285, nnz 285
		%>- FACTOR 3 (real) SPARSE, size 100x100, density 0.0282, nnz 282
		%>- FACTOR 4 (real) SPARSE, size 100x100, density 0.0292, nnz 292
		%>- FACTOR 5 (real) SPARSE, size 100x50, density 0.02, nnz 100
		%>>>[F;G;F;G] % it's allowed to concatenate an arbitrary number of Fausts
		%>>>[F,G,F,G] % as long as the dimensions agree
		%>>>[F,G;F,G]
		%> @endcode
		%>
		%> @b Errors
		%> - The dimensions of F and G don't agree:
		%>
		%> @code
		%>>> F = matfaust.rand(2,51);
		%>>> G = matfaust.rand(2,25);
		%>>> [F;G]
		%>Error using mexFaustReal
		%>The dimensions of the two Fausts must agree.
		%>
		%> @endcode
		%> - The DIM argument is different from 1 and 2:
		%>
		%> @code
		%>>> F = matfaust.rand(2,4);
		%>>> G = matfaust.rand(2,4);
		%>>> cat(3,F,G)
		%>Error using matfaust.Faust/cat
		%>Wrong first argument: must be an integer between 1 and 2.
		%> @endcode
		%> <p>@b See @b also Faust.vertcat, Faust.horzcat.
		%>
		%======================================================================
		function C = cat(varargin)
			err_1st_arg = 'Wrong first argument: must be an integer between 1 and 2.';
			if(nargin > 0 && isscalar(varargin{1}))
				F = varargin{2}; % we kwnow it's a Faust or we wouldn't be here
				if(varargin{1} == 1)
					mex_func_name = 'vertcat';
				elseif(varargin{1} == 2)
					mex_func_name = 'horzcat';
				else
					error(err_1st_arg)
				end
				C = F;
				for i=3:nargin
					A = varargin{i};
					if(ismatrix(A) && ~ isa(A, 'matfaust.Faust'))
						A = matfaust.Faust({A});
					end
					if(~ isa(A, 'matfaust.Faust'))
						error('Can''t concatenate a Faust to something that is not a Faust or a matrix.')
					end
					if(C.isReal)
						if(~ isreal(A))
							C = complex(C)
							C = matfaust.Faust(C, mexFaustCplx(mex_func_name, C.matrix.objectHandle, A.matrix.objectHandle));
						else
							C = matfaust.Faust(C, mexFaustReal(mex_func_name, C.matrix.objectHandle, A.matrix.objectHandle));
						end
					else
						if(isreal(A))
							A = complex(A)
						end
						C = matfaust.Faust(C, mexFaustCplx(mex_func_name, C.matrix.objectHandle, A.matrix.objectHandle));
					end
				end
			else
				error(err_1st_arg)
			end
		end

		%======================================================================
		%> Horizontal concatenation [F,A, B, … ] where F is a Faust and A, B, … are Faust objects or sparse/full matrix.
		%===
		%> It's equivalent to cat(2, F, A, B, …).
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b C=HORZCAT(F,A) concatenates the Faust F and A horizontally. A is a Faust or a sparse/full matrix. The result is the Faust C. @b HORZCAT(F,A) is the same as [F,G].<br/>
		%> &nbsp;&nbsp;&nbsp; @b C=HORZCAT(@b F,@b A,@b B,…) concatenates horizontally F to A, B, etc. @b HORZCAT(@b F,@b A, @b B,…) is the same as [F,A,B,…].<br/>
%>
		%> <p>@b See @b also Faust.vertcat, Faust.cat.
		%======================================================================
		function HC = horzcat(varargin)
			HC = cat(2, varargin{1}, varargin{2:end});
		end

		%======================================================================
		%> Vertical concatenation [F;A;B;…] where F is a Faust and A, B, … are Faust objects or sparse/full matrices.
		%===
		%>
		%> It's equivalent to cat(1, F, A, B, …).
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b C=VERTCAT(F,A) concatenates the Faust F and A vertically. A is a Faust or a sparse/full matrix. The result is the Faust C. @b VERTCAT(F,A) is the same as [F;G].<br/>
		%> &nbsp;&nbsp;&nbsp; @b C=VERTCAT(@b F, @b A, @b B,…) concatenates vertically F to A, B etc. @b VERTCAT(@b F, @b A, @b B,…) is the same as [F;A;B;…].<br/>
		%>
		%> <p>@b See @b also Faust.horzcat, Faust.cat.
		%======================================================================
		function VC = vertcat(varargin)
			VC = cat(1, varargin{1}, varargin{2:end});
		end

		%======================================================================
		%> Function not implemented in the Faust class.
		%>
		%======================================================================
		function RF = reshape(varargin)
			error('Function not implemented in the Faust class.');
		end


		%=====================================================================
		%> Displays image of F's full matrix and its factors.
		%===
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b imagesc(F)<br/>
		%> &nbsp;&nbsp;&nbsp; @b imagesc(F, 'faust_name') to name the Faust on the plotted figure.
		%>
		%> Data is scaled to use the full colormap.
		%>
		%> @b Example
		%> @code
		%>	F = matfaust.rand(5, [50, 100], .5, 'mixed', 'complex')
		%>	imagesc(F)
		%>	print('test_imsc', '-dpng')
		%>	@endcode
		%>
		%> <p>@b See @b also Faust.disp.
		%======================================================================
		function imagesc(F, varargin)
			numfacts = numfactors(F);
			facs = cell(1,numfacts+1);
			for i=1:numfacts
				facs{i} = factors(F,i);
			end
			facs{i+1} = full(F);
			numfacts = numfacts + 1;
			fig = figure();
			name = 'F';
			if(nargin > 1)
				name = varargin{1};
				if(~ isstr(name))
					error('argument 2 must be a str.')
				end
			end
			facs_per_row = 5;
			set(gca,'XTick',[], 'YTick', []);
			num_rows = numfacts/facs_per_row;
			if(floor(num_rows) < num_rows)
				num_rows = floor(num_rows)+1;
			elseif(num_rows == 0)
				num_rows = numfacts;
			end
			row_widths = zeros(1,num_rows);
			% compute widths of all rows
			for i=1:numfacts
				row_widths(1, floor((i-1)/facs_per_row)+1) = row_widths(1, floor((i-1)/facs_per_row)+1) + size(facs{i}, 2);
			end
			maxw = max(row_widths);
			sumw = maxw;
			sumw = sumw+(facs_per_row)*sumw/100; % extra with for space between subplots
			cumw = (1 - (row_widths(1,floor(1/facs_per_row)+1))/sumw - (facs_per_row-1)*.01)/2;

			figpos = get(fig, 'Position');
			aratio = figpos(3)/figpos(4);
			maxh = 0;
			for i=1:facs_per_row
				fac = facs{i};
				h = size(fac,1)/sumw*aratio;
				if(h > maxh)
					maxh = h;
				end
			end
			b = 1 - maxh - 0.01;
			for i=1:numfacts
				fac = facs{i};
				l = cumw;
				w = size(fac,2)/sumw;
				h = size(fac,1)/sumw*aratio;
				pos = [l b w h];
				subplot('position', pos)
				cumw = cumw + w + .01; % extra percent for space
				imagesc(abs(real(fac)));
				if(i == numfacts)
					text(size(fac,2)/2, size(fac,1)/2, 'full(F)')
				else
					text(size(fac,2)/2, size(fac,1)/2, [int2str(i)])
				end
				set(gca,'XTick',[], 'YTick', []);
				if(mod(i,facs_per_row) == 0 && i < numfacts)
					cumw = (1 - (row_widths(1,floor(i/facs_per_row)+1))/sumw - (facs_per_row-1)*.01)/2;
					maxh = 0;
					for j=i+1:min(facs_per_row+i,numfacts)
						fac = facs{j};
						h = size(fac,1)/sumw*aratio;
						if(h > maxh)
							maxh = h;
						end
					end
					b = b - maxh - 0.01;
				end
			end
			suptitle(['Factor of the Faust ' name ])
		end

		%=====================================================================
		%> \   Backslash or left full(F) divide.
		%===
		%>
		%> @b Usage
		%>
		%>  &nbsp;&nbsp;&nbsp; <b> X = F\ B </b> is the matrix division of full(F) into B, which is roughly the
		%>          same as Faust.pinv(F)*B.
		%>
		%> <p> @b See @b also Faust.pinv, mldivide Matlab built-in.
		%=====================================================================
		function X = mldivide(F,B)
			X = mldivide(full(F),B)
		end

		%=====================================================================
		%> Pseudoinverse matrix of full(F).
		%===
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%>  &nbsp;&nbsp;&nbsp; <b> X = PINV(F) </b> produces the matrix X of the same dimensions as F'
		%>  so that F*(F'*X')' == full(F) (or approximately).
		%>
		%> <p> @b See @b also Faust.mldivide, pinv Matlab built-in.
		%=====================================================================
		function X = pinv(F)
			X = pinv(full(F))
		end

		%================================================================
		%> Converts F to a complex Faust.
		%===
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; <b> cF = COMPLEX(F) </b> for real Faust F returns the complex result cF
		%> with real part F and zero matrix as imaginary part of all cF's factors.
		%>
		%> <p> @b See @b also Faust.isreal, Faust.conj
		%================================================================
		function cF = complex(F)
			if(~ isreal(F))
				cF = F;
			else
				n = numfactors(F);
				facs = cell(1,n);
				for i=1:n
					facs{i} = factors(F,i);
					if(issparse(facs{i}))
						% avoid complex() error: Input A must be numeric and full.
						facs{i} = sparse(complex(full(facs{i})));
					else
						facs{i} = complex(facs{i});
					end
				end
				cF = matfaust.Faust(facs);
			end
		end

	end
	methods(Access = private)
		%================================================================
		%> Returns true if i is a valid factor index for the Faust F.
		%===
		%================================================================
		function  i = check_factor_idx(F,i)
			if(~ isnumeric(i) || ~ isscalar(i))
				error('i must be a scalar.')
			end
			if( i <= 0 || i > numfactors(F))
				error('i is out of range')
			end
			i = floor(i);

		end
	end
	methods(Static)
		%================================================================
		%> Returns true if obj is a Faust object, false otherwise.
		%===
		%>
		%> @b Example
		%> @code
		%> import matfaust.*
		%> Faust.isFaust(1) % returns 0
		%> Faust.isFaust(matfaust.rand(5,10)) % returns 1
		%> @endcode
		%>
		%> <p> @b See @b also Faust.Faust
		%================================================================
		function bool = isFaust(obj)
			bool = isa(obj, 'matfaust.Faust');
		end
	end
end

