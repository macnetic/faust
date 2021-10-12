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


% ======================================================================
%> @brief <b>FAuST Matlab wrapper main class</b> for using multi-layer sparse transforms.
%>
%> This class provides a Matlab array-like interface for operations with FAuST data structures, which correspond to matrices that can be written exactly as the product of sparse matrices.
%>
%> The Faust class is designed to allow fast matrix-vector multiplications together with reduced memory storage compared to what would be obtained by manipulating directly the corresponding (dense) Matlab array.
%>
%> A particular example is the matrix associated to the discrete Fourier transform, which can be represented exactly as a Faust, leading to a fast and compact implementation (see matfaust.dft()).
%>
%> Although sparse matrices are more interesting for optimization it's not forbidden to define a Faust as a product of dense matrices or a mix of dense and sparse matrices.
%>
%> The matrices composing the Faust product, also called the factors, are defined on complex or real fields. Hence a Faust can be a complex Faust or a real Faust.
%>
%> Several Matlab builtins have been overloaded to ensure that a Faust is
%> almost handled as a native Matlab matrix.
%>
%> The main exception is that contrary to a Matlab native array a Faust is immutable.
%> It means that you cannot modify elements of a Faust using
%> the assignment operator `=' like you do with a Matlab matrix (e.g. `M(i,j) =
%> 2').
%> That limitation is the reason why the Matlab built-in `SUBSASGN()' is not
%> implemented in this class.
%> Note however that you can optionally contravene the immuability in certain functions (e.g. with the `inplace` argument of the
%> Faust.optimize_time).
%>
%> Other noticeable limitations are that one cannot:
%> - compute the real and imaginary parts of a Faust,
%> - perform elementwise operations between two Fausts (e.g. elementwise
%> multiplication), the addition and subtraction are available though,
%> - reshape a Faust.
%>
%> Mainly for convenience and test purposes, a Faust can be converted into
%> the corresponding full matrix using the function Faust.full.
%>
%> @warning using Faust.full is discouraged except for test purposes, as it
%> loses the main potential interests of the FAuST structure: compressed
%> memory storage and faster matrix-vector multiplication compared to its
%> equivalent full matrix representation.
%>
%> In this documentation, the expression 'full matrix' designates the Matlab array
%> Faust.full() obtained by the multiplication of the Faust factors.
%>
%> List of functions that are memory costly:
%> - Faust.full(),
%> - Faust.pinv(),
%> - Faust.mldivide().
%> - element indexing (F(i,j) / __getitem__, but note that slicing
%>   is memory efficient through memory views).
%>
%>
%> For more information about FAuST take a look at http://faust.inria.fr.
%>
% ======================================================================

classdef Faust
	properties (SetAccess = protected, Hidden = true)
		matrix; % Handle to the FaustCore class instance
		isReal;
		dev; % cpu or gpu
		dtype; % double (for complex or real), or float
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
			%%
			err_msg = 'matfaust.Faust() error: the arguments are not valid.';
			% scale argument is hidden for user (deprecated) but it's still available
			max_nargin = 5;
			if(nargin < 1 || nargin > max_nargin)
				error([err_msg ' Number of arguments passed is zero or greater than ' str(max_nargin) '.'])
			elseif(iscell(varargin{1}))
				% init a Faust from factor list
				% check if the factors are real or complex, one complex factor implies a complex faust
				factors=varargin{1};
				if(length(factors) == 0)
					error([ err_msg ' Cannot create an empty Faust.'])
				end
				% default argument values
				F.isReal = true;
				scale = 1;
				F.dev = 'cpu';
				F.dtype = 'double';
				optCopy = false;
				for i=1:length(factors)
					if (~isreal(factors{i}))
						F.isReal = false;
						break
					end
				end
				if(strcmp(class(factors{1}), 'single'))
					F.dtype = 'float';
				end
				for i=2:nargin
					switch(varargin{i})
						case 'scale'
							if(nargin < i+1 || ~ isscalar(varargin{i+1}) || ~ isnumeric(varargin{i+1}))
								error([err_msg ' the ''scale'' argument must be followed by a scalar.'])
							end
							scale = varargin{i+1};
							if(~ isreal(scale))
								F.isReal = false;
							end
						case 'dev'
							if(nargin < i+1 || ~ strcmp(varargin{i+1}, 'cpu') && ~ strcmp(varargin{i+1}, 'gpu'))
								error([err_msg ' the ''dev'' argument must be followed by ''cpu'' or ''gpu'''])
							end
							F.dev = varargin{i+1};
						case 'optCopy'
							if(nargin < i+1 || ~ islogical(varargin{i+1}))
								error([err_msg ' the ''optCopy'' argument must be followed by a logical'])
							end
							optCopy = varargin{i+1};
						case 'dtype'
							if(nargin < i+1 || ~ strcmp(varargin{i+1}, 'double') && ~ strcmp(varargin{i+1}, 'float'))
								error([err_msg ' the ''dtype'' argument must be followed by ''float'' or ''double'''])
							end
							F.dtype = varargin{i+1};
						otherwise
							if(isstr(varargin{i}) && ~ strcmp(varargin{i}, 'cpu') && ~ strcmp(varargin{i}, 'gpu') && (~ strcmp(varargin{i}, 'float') && ~ strcmp(varargin{i}, 'double')))
								error([ varargin{i} ' unrecognized argument'])
							end
						end

					end
				F.matrix = FaustCore(factors, scale, optCopy, F.isReal, F.dev, F.dtype);
			elseif(ischar(varargin{1}))
				% init a Faust from file
				filename=varargin{1};
				load(filename);
				if (~exist('faust_factors','var') )
					error('Faust : invalid file');
				end
				F = matfaust.Faust(faust_factors, varargin{2:end});
			elseif(ismatrix(varargin{1}) && nargin == 1)
				% create a Faust from a matrix (single factor)
				c = cell(1, 1);
				c{1} = varargin{1};
				F = matfaust.Faust(c, varargin{2:end});
			elseif(isa(varargin{1}, 'matfaust.Faust'))
				% create a Faust from another one but not with the same
				% handle to set inside the FaustCore object (matrix)
				oF = varargin{1};
				F.matrix = FaustCore(varargin{2}, oF.isReal, oF.dev, oF.dtype);
				F.isReal = oF.isReal;
				F.dev = oF.dev;
			elseif(isa(varargin{1}, 'integer') && islogical(varargin{2}))
				% create a Faust directly with the c++ object handler/pointer without any pre-existing Faust
				F.isReal = varargin{2};
				if(nargin >= 3)
					F.dev = varargin{3};
				else
					F.dev = 'cpu';
				end
				if(nargin >= 4)
					F.dtype = varargin{4};
				else
					F.dtype = 'double';
				end
				F.matrix = FaustCore(varargin{1}, varargin{2}, F.dev, F.dtype);
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
		%>	F = matfaust.rand(5, 10)
		%>	delete(F)
		%>	F = matfaust.rand(5, 10)
		%>	G = matfaust.rand(5, 10)
		%>	clear % equivalent to delete(F);delete(G)
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust, clear (built-in)</p>
		%>
		%======================================================================
		function delete(F)
			%%
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
		%>   F = matfaust.rand(5, 10)
		%>   A = rand(size(F,2), 50)
		%>   G = F*A
		%> % is equivalent to G = mtimes(F, A)
		%>   G = matfaust.rand(size(F, 2), 15)
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
			G = mtimes_trans(F, A, 0);
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
			%%

			if ~isreal(trans)
				error('invalid argument trans, must be equal to 0 or 1');
			end

			if (trans ~= 1) && (trans ~= 0)
				error('invalid argument trans, must be equal to 0 or 1');
			end
			% TODO: take trans into account when mul F to a scal or a Faust
			% it's not a serious issue because mtimes_trans() shouldn't be called by final user
			% (func's doc is filtered out in doxydoc)

			if(isa(A, 'matfaust.Faust') && ~ isa(F, 'matfaust.Faust'))
				if(ismatrix(F))
					C = mtimes_trans(A', F', 0)';
				else
					C = mtimes_trans(A, F, 0); % F is scalar or something else (in which case it'll fail later)
				end
			elseif(isa(A,'matfaust.Faust'))
				if (F.isReal)
					if(isreal(A))
						C = matfaust.Faust(F, call_mex(F, 'mul_faust', A.matrix.objectHandle));
					else
						C = mtimes_trans(complex(F), A, trans);
					end
				else
					if(A.isReal)
						A = complex(A);
					end
					C = matfaust.Faust(F, call_mex(F, 'mul_faust', A.matrix.objectHandle));
				end
			elseif(isscalar(A))
				if (F.isReal)
					if(isreal(A))
						C = matfaust.Faust(F, call_mex(F, 'mul_scalar', A));
					else
						C = mtimes_trans(complex(F), A, trans);
					end
				else
					C = matfaust.Faust(F, call_mex(F, 'mul_scalar', A));
				end
			elseif (F.isReal) % A is not a Faust (should be a matrix)
				if (isreal(A))
					C = call_mex(F, 'multiply', A, trans);
				else
					C_real = call_mex(F, 'multiply', real(A), trans);
					C_imag = call_mex(F, 'multiply', imag(A), trans);
					C = C_real + 1i * C_imag;
				end
			else
				C = call_mex(F, 'multiply', A, trans);
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
		%>   >> F = matfaust.rand(10^5, 10^5, 'num_factors', 2, 'density', 10^-4, 'fac_type', 'sparse')
		%>   Faust size 100000x100000, density 0.00018, nnz_sum 1800000, 2 factor(s):
		%>   - FACTOR 0 (real) SPARSE, size 100000x100000, density 9e-05, nnz 900000
		%>   - FACTOR 1 (real) SPARSE, size 100000x100000, density 9e-05, nnz 900000
		%>   >> % an attempt to convert F to a full matrix is most likely to raise a memory error
		%>   >> % the sparse format is the only way to handle such a large Faust
		%>   >> full(F)
		%>   Out of Memory
		%> @endcode
		%>
		%======================================================================
		function A = full(F)
			A = call_mex(F, 'full');
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
			%% 
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
			F_trans = matfaust.Faust(F, call_mex(F, 'transpose'));
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
		%>	F = matfaust.rand(5, 10)
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
			F_ctrans = matfaust.Faust(F, call_mex(F, 'ctranspose'));
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
		%>	F = matfaust.rand(5, 10)
		%>	F_conj = conj(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.transpose, Faust.ctranspose, Faust.complex
		%>
		%======================================================================
		function F_conj = conj(F)
			F_conj = matfaust.Faust(F, call_mex(F, 'conj'))
		end

		%======================================================================
		%> @brief Returns a Faust optimized by removing useless zero rows and columns as many times as needed.
		%>
		%> @param F: the Faust to optimize.
		%> @param 'thres', int: (optional) the threshold of number of nonzeros under what the
		%>            rows/columns are removed.
		%>
		%> @retval G The optimized Faust.
		%======================================================================
		function G = pruneout(F, varargin)
			thres = 0;
			% hidden parameters: useless except for debug
			% @param 'only_forward', bool: (optional) True for applying only the forward passes of removal.
			% @param 'npasses', int: (optional) the number of passes to run, by default it goes until the
			%            optimal Faust is obtained.
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
						case 'thres'
							if(argc == i || ~ isnumeric(varargin{i+1}) || varargin{i+1}-floor(varargin{i+1}) > 0 || varargin{i+1} <0 )
								error('thres keyword arg. is not followed by a positive integer')
							end
							thres = varargin{i+1};
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
			G = matfaust.Faust(F, call_mex(F, 'pruneout', thres, npasses, only_forward));
		end

		%=====================================================================
		%> @brief Optimizes a Faust by changing the storage format of each factor in
		%> order to optimize the memory size.
		%>
		%>@retval OF The optimized Faust.
		%>
		%> @b See @b also Faust.optimize
		%=====================================================================
		function OF = optimize_memory(F)
			OF = matfaust.Faust(F, call_mex(F, 'optimize_storage', false));
		end

		%=====================================================================
		%> @brief Returns a Faust optimized with Faust.pruneout, Faust.optimize_memory and Faust.optimize_time.
		%>
		%> @b Usage <br/>
		%> &nbsp;&nbsp;&nbsp; @b OF = @b optimize_time(F) returns an optimized Faust object.<br/>
		%> &nbsp;&nbsp;&nbsp; @b OF = @b optimize_time(@b F, @b 'transp', @b true) see Faust.optimize_time
		%>
		%> @param 'transp', bool (optional) default to false. if true optimize the Faust according to its transpose.
		%>
		%> @retval OF The optimized Faust.
		%>
		%> @note this function is still experimental, you might use manually
		%> Faust.optimize_time, Faust.optimize_memory or Faust.pruneout to be
		%> more specific about the optimization to proceed.
		%>
		%>
		%> @b See @b also Faust.optimize_time, Faust.optimize_memory, Faust.pruneout
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
			OF = matfaust.Faust(F, call_mex(F, 'optimize', transp));
		end

		%===============================
		%> @brief Returns a Faust configured with the quickest Faust-matrix multiplication mode (benchmark ran on the fly).
		%> @note this function launches a small benchmark on the fly. Basically, the methods
		%> available differ by the order used to compute the matrix chain
		%> multiplication or by the use (or unuse) of threads for the calculation of intermediary
		%> matrix products of the Faust.
		%>
		%>
		%> @b Usage<br/>
		%> &nbsp;&nbsp;&nbsp; @b OF = @b optimize_time(F) returns a new object with the best product method enabled.<br/>
		%> &nbsp;&nbsp;&nbsp; @b OF = @b optimize_time(@b F, @b 'inplace', @b true) modify directly F instead of creating a new Faust (i.e. OF references the same object as F).<br/>
		%> &nbsp;&nbsp;&nbsp; @b OF = @b optimize_time(@b F, @b 'transp', @b true) try to optimize the Faust transpose instead of the Faust itself. <br/>
		%> &nbsp;&nbsp;&nbsp; @b OF = @b optimize_time(@b F, @b 'nsamples', @b 10) benchmark product methods by computing 10 product instead of one by default.
		%>
		%> @param F the Faust object.
		%> @param 'inplace', bool (optional) default to false. If true the current Faust is modified directly.
		%> @param 'transp', bool (optional) default to false. If true optimize the Faust according to its transpose.
		%> @param 'nsamples', int (optional) default to 1.The number of Faust-Dense matrix products
		%> calculated in order to measure time taken by each method (it could matter
		%> to discriminate methods when the performances are similar). By default,
		%> only one product is computed to evaluate the method.
		%> @param 'mat', matrix (optional) Use this argument to run the benchmark on the Faust multiplication by the matrix mat instead of Faust.full(). Note that mat must be of the same scalar type as F.
		%>
		%> @retval OF The optimized Faust.
		%>
		%> @b See @b also Faust.optimize
		%===============================
		function OF = optimize_time(F, varargin)
			transp = false;
			inplace = false;
			nsamples = 1;
			argc = length(varargin);
			mat = false;
			if(argc > 0)
				i = 1;
				while(i < argc)
					switch(varargin{i})
						case {'transp', 'transpose'}
							if(argc == i || ~ islogical(varargin{i+1}))
								error('transp keyword argument is not followed by a logical')
							else
								transp = varargin{i+1};
							end
						case 'inplace'
							if(argc == i || ~ islogical(varargin{i+1}))
								error('inplace keyword argument is not followed by a logical')
							else
								inplace = varargin{i+1};
							end
						case 'nsamples'
							if(argc == i || ~ isscalar(varargin{i+1}) || floor(varargin{i+1}) < varargin{i+1})
								error('nsamples keyword argument is not followed by an integer')
							else
								nsamples = varargin{i+1}
							end
						case 'mat'
							if(argc == i || ~ isnumeric(varargin{i+1}) || ~ ismatrix(varargin{i+1}))
								error('mat keyword argument is not followed by a matrix.')
							else
								mat = varargin{i+1};
								i = i + 1; % ignore mat from parsing (switch can handle only scalar or char vec
							end
						otherwise
							if(isstr(varargin{i}))
								error([ varargin{i} ' unrecognized argument'])
							end
					end
					i = i + 1;
				end
			end
			args = {transp, inplace, nsamples};
			mex_func = 'optimize_time';
			if(~ islogical(mat) && ismatrix(mat))
				% mat is a matrix on which to run the benchmark
				args = [ args {mat} ];
				mex_func = 'optimize_time_prod';
			end
			if(inplace)
				call_mex(F, mex_func, args{:});
				OF = F
			else
				OF = matfaust.Faust(F, call_mex(F, mex_func, args{:}));
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
		%>	F = matfaust.rand(5, 10)
		%>	[nrows, ncols] = size(F)
		%>	nrows = size(F, 1)
		%>	ncols = size(F, 2)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.nnz_sum, Faust.numel
		%======================================================================
		function varargout = size(F,varargin)
			%%

			nb_input = length(varargin);


			if (nb_input > 1)
				error('Too many input arguments');
			end

			if ((nb_input == 1) && (nargout > 1) | (nargout > 2))
				error('Too many output arguments');
			end
			Size=[-1 -1];
			Size = call_mex(F, 'size');

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
		%>	F = matfaust.rand(5, 10)
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
		%> @brief Gives the memory size of the Faust in bytes.
		%>
		%>
		%> @retval n Faust size in bytes.
		%======================================================================
		function n = nbytes(F)
			n = call_mex(F, 'nbytes');
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
		%>	>> F = matfaust.rand(3, 3);
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
			%%

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
		%>	F = matfaust.rand(5, 10);
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
				factors{j} = call_mex(F, 'factors', i);
			end
			if(length(factors) > 1)
				factors = matfaust.Faust(factors);
			else
				factors = factors{j};
			end
		end

		%================================================================
		%> Returns the left hand side factors of F from index 1 to i included (in 1-base index).
		%===
		%> @Example
		%> @code
		%> %in a matlab terminal
		%> >> F = matfaust.rand(8, 5)
		%>
		%> F =
		%>
		%> Faust size 8x5, density 5.25, nnz_sum 210, 5 factor(s):
		%> - FACTOR 0 (real) SPARSE,  size 8x6, density 0.833333, nnz 40
		%> - FACTOR 1 (real) SPARSE,  size 6x9, density 0.555556, nnz 30
		%> - FACTOR 2 (real) SPARSE,  size 9x10, density 0.5, nnz 45
		%> - FACTOR 3 (real) SPARSE,  size 10x9, density 0.555556, nnz 50
		%> - FACTOR 4 (real) SPARSE,  size 9x5, density 1, nnz 45
		%> >> LF = left(F, 3)
		%>
		%> LF =
		%>
		%> Faust size 8x10, density 1.4375, nnz_sum 115, 3 factor(s):
		%> - FACTOR 0 (real) SPARSE,  size 8x6, density 0.833333, nnz 40
		%> - FACTOR 1 (real) SPARSE,  size 6x9, density 0.555556, nnz 30
		%> - FACTOR 2 (real) SPARSE,  size 9x10, density 0.5, nnz 45
		%>
		%>@endcode
		%> <p> @b See @b also Faust.factors, Faust.right
		%================================================================
		function lfactors = left(F, i)
			i = check_factor_idx(F, i);
			lfactors = factors(F, 1:i);
		end

		%================================================================
		%> Returns the right hand side factors of F from index i to end (in 1-base index).
		%===
		%>
		%> @Example
		%> @code
		%> >> F = matfaust.rand(7,7, 'dim_sizes', [7, 10])
		%>
		%> F =
		%>
		%> Faust size 7x7, density 4.28571, nnz_sum 210, 5 factor(s):
		%> - FACTOR 0 (real) SPARSE, size 7x10, density 0.5, nnz 35
		%> - FACTOR 1 (real) SPARSE,  size 10x9, density 0.555556, nnz 50
		%> - FACTOR 2 (real) SPARSE, size 9x8, density 0.625, nnz 45
		%> - FACTOR 3 (real) SPARSE, size 8x8, density 0.625, nnz 40
		%> - FACTOR 4 (real) SPARSE,  size 8x7, density 0.714286, nnz 40
		%> >> RF = right(F, 2)
		%>
		%> RF =
		%>
		%> Faust size 10x7, density 2.5, nnz_sum 175, 4 factor(s):
		%> - FACTOR 0 (real) SPARSE,  size 10x9, density 0.555556, nnz 50
		%> - FACTOR 1 (real) SPARSE, size 9x8, density 0.625, nnz 45
		%> - FACTOR 2 (real) SPARSE, size 8x8, density 0.625, nnz 40
		%> - FACTOR 3 (real) SPARSE,  size 8x7, density 0.714286, nnz 40
		%>
		%> @endcode
		%>
		%> <p> @b See @b also Faust.factors, Faust.left
		%================================================================
		function rfactors = right(F, i)
			i = check_factor_idx(F, i);
			rfactors = factors(F, i:numfactors(F));
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
		%>	F = matfaust.rand(5, 10)
		%>	nf = numfactors(F)
		%> @endcode
		%>
		%> <p>@b See @b also Faust.factors.
		%==========================================================================================
		function num_factors = numfactors(F)
			num_factors = call_mex(F, 'numfactors');
		end

		%==========================================================================================
		%> @brief Returns true if F factors are all sparse matrices false otherwise.
		%>
		%==========================================================================================
		function is_sparse = issparse(F)
			is_sparse = call_mex(F, 'is_all_sparse');
		end

		%==========================================================================================
		%> @brief Returns true if F factors are all dense matrices/arrays false otherwise.
		%>
		%==========================================================================================
		function is_dense = isdense(F)
			is_dense = call_mex(F, 'is_all_dense');
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
		%>	F = matfaust.rand(5, 10)
		%>	save(F, 'F.mat')
		%>	G = Faust('F.mat')
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust, Faust.rcg.
		%==========================================================================================
		function save(F, filepath)
			call_mex(F, 'save', filepath);
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
		%>	F = matfaust.rand(50, 100)
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
			%%

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

			if(indexing_by_slice)
				submatrix = matfaust.Faust(F, call_mex(F, 'subsref', start_ids(ROW), end_ids(ROW), start_ids(COL), end_ids(COL)));
			else
				% -1 is for converting to 0-base index used in C world
				submatrix =  matfaust.Faust(F, call_mex(F, 'subsref_byvec', uint64(ind_lists{1}-1), uint64(ind_lists{2}-1)));
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
		%>	>> F = matfaust.rand(98, 82)
		%>	>> disp(F)
		%>	Faust size 98x82, density 0.115729, nnz_sum 930, 2 factor(s):
		%>	- FACTOR 0 (real) SPARSE, size 98x88, density 0.0568182, nnz 490
		%>	- FACTOR 1 (real) SPARSE, size 88x82, density 0.0609756, nnz 440
		%>
		%>	>> F
		%>	Faust size 98x82, density 0.115729, nnz_sum 930, 2 factor(s):
		%>	- FACTOR 0 (real) SPARSE, size 98x88, density 0.0568182, nnz 490
		%>	- FACTOR 1 (real) SPARSE, size 88x82, density 0.0609756, nnz 440
		%> @endcode
		%>
		%> <p>@b See @b also Faust.nnz_sum, Faust.density, Faust.size, Faust.factors, Faust.numfactors
		%>
		%>
		%======================================================================
		function disp(F)
			call_mex(F, 'disp');
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
		%> @warning The norm computation time can be expected to be of order
		%> n*min(F.shape) with n the time for multipliying F by a vector.
		%> Nevertheless, the implementation allows that memory usage remains
		%> controlled by avoiding to explicitly compute full(F). Please pay
		%> attention to the full_array (and batch_size) arguments for a better
		%> understanding.

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
		%> @param 'threshold',real (optional) power iteration algorithm threshold (default to .001). Used only for norm(2). It's passed in a key-value pair fashion: 'threshold', .001
		%> @param 'max_num_its',int (optional) maximum number of iterations for power iteration algorithm. Used only for norm(2). It's passed in a key-value pair fashion: 'max_num_its', 1000.
		%> @param 'full_array',bool (optional) this argument applies only for 1-norm,
		%> inf-norm and Frobenius norm. If true the Faust full array
		%> is computed before computing the norm otherwise it is not. By
		%> default it is set to false. Many configurations exist in which
		%> full_array == False can be more efficient but it needs to
		%> finetune the batch_size argument.
		%> @param 'batch_size',int (optional) this argument applies only when the
		%> full_array argument is set to false (for the 1-norm, inf-norm and
		%> Frobenius norm). It determines the number of Faust columns (resp. rows)
		%> that are built in memory in order to compute the Frobenius norm and
		%> the 1-norm (resp. the inf-norm). This parameter is primary in the
		%> efficiency of the computation and memory consumption. By  default,
		%> it is set to 1 (which is certainly not the optimal configuration in
		%> many cases in matter of computation time but always the best in
		%> term of memory cost).
		%>
		%> @retval n the norm (real).
		%>
		%>
		%> @b Example
		%> @code
		%> % in a matlab terminal
		%> >> F = matfaust.rand(50, 100, 'num_factors', 2, 'dim_sizes', [50, 100], 'density', .5)
		%> >> norm(F)
		%> ans =
		%> 436.4094
		%> >> norm(F,2)
		%> ans =
		%> 436.4094
		%> >> norm(F, 'fro')
		%> ans =
		%> 442.3946
		%> >> norm(F, inf)
		%> ans =
		%> 738.0057
		%> @endcode
		%>
		%======================================================================
		function n = norm(F, varargin)
			%%
			nargs = length(varargin);
			% set parameter default values
			ord = 2;
			ord2_valid_param = false;
			threshold = 1e-3;
			max_num_its = 100;
			batch_size = 1; % warning: batch_size and full_array default values are eventually calculated dynamically (see in the next)
			full_array = true;
			full_array_is_default = true;
			batch_size_is_default = true;
			if nargs >= 1
				ord = varargin{1};
				if ~ strcmp(ord, 'fro') && ord ~= 2 && ord ~= 1 && ord ~= inf
					error('only 1, 2, inf or Frobenius norms are supported for a Faust');
				end
				for i=2:nargs
					switch(varargin{i})
						case 'threshold'
							tmparg = varargin{i+1};
							if(nargs == i || ~ isnumeric(tmparg) || ~ isreal(tmparg) || ~ isscalar(tmparg) || tmparg < 0)
								error('threshold keyword arg. is not followed by a positive real number')
							end
							threshold = varargin{i+1};
						case 'max_num_its'
							if(nargs == i ||  ~ isnumeric(varargin{i+1}) || varargin{i+1}-floor(varargin{i+1}) > 0)
								error('max_num_its keyword arg. is not followed by an integer')
							else
								max_num_its = varargin{i+1};
							end
						case 'full_array'
							if(nargs == i || ~ islogical(varargin{i+1}))
								error('full_array keyword arg. is not followed by a logical')
							else
								full_array_is_default = false;
								full_array = varargin{i+1};
							end
						case 'batch_size'
							if(nargs == i ||  ~ isnumeric(varargin{i+1}) || varargin{i+1}-floor(varargin{i+1}) > 0)
								error('batch_size keyword arg. is not followed by an integer')
							else
								batch_size = varargin{i+1};
								batch_size_is_default = false;
							end
						otherwise
							if(isstr(varargin{i}))
								error([ varargin{i} ' unrecognized argument'])
							end
					end
				end
			end
			if(ord == 2)
				extra_opts = {threshold, max_num_its};
			else % fro, 1, inf
				extra_opts = {full_array, batch_size};
			end
			n = call_mex(F, lower([num2str(ord) 'norm']), extra_opts{:});
		end

		%==========================================================
		%> Performs the power iteration algorithm to compute the greatest eigenvalue of the Faust.
		%>
		%> For the algorithm to succeed the Faust should be diagonalizable
		%> (similar to a digonalizable Faust), ideally, a symmetric positive-definite Faust.
		%>
		%> @param 'threshold', real: (optional) the precision required on the eigenvalue. Default value is 1e-3.
		%> @param 'maxiter', integer: (optional) the number of iterations above what the algorithm will stop anyway. Default value is 100.
		%>
		%> @retval lambda the greatest eigenvalue approximate.
		%>
		%> @Example
		%> @code
		%> %in a matlab terminal
		%> >> F = matfaust.rand(8, 5)
		%> >> F = F*F'
		%> >> power_iteration(F)
		%> 1.1795e+04
		%> @endcode
		%==========================================================
		function lambda = power_iteration(F, varargin)
			threshold = 1e-3;
			maxiter = 100;
			argc = length(varargin);
			if(argc > 0)
				for i=1:argc
					switch(varargin{i})
						case 'threshold'
							tmparg = varargin{i+1};
							if(argc == i || ~ isnumeric(tmparg) || ~ isreal(tmparg) || ~ isscalar(tmparg) || tmparg < 0)
								error('threshold keyword arg. is not followed by a positive real number')
							end
							threshold = varargin{i+1};
						case 'maxiter'
							if(argc == i ||  ~ isnumeric(varargin{i+1}) || varargin{i+1}-floor(varargin{i+1}) > 0)
								error('maxiter keyword arg. is not followed by an integer')
							else
								maxiter = varargin{i+1};
							end
						otherwise
							if(isstr(varargin{i}))
								error([ varargin{i} ' unrecognized argument'])
							end
					end
				end
			end
			lambda = call_mex(F, 'power_ite', threshold, maxiter);
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
			NF = matfaust.Faust(F, call_mex(F, 'normalize', mex_ord));
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
			nz = call_mex(F, 'nnz');
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
		%>	F = matfaust.rand(5, 10)
		%>	dens = density(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.nnz_sum, Faust.rcg, Faust.size, Faust.numel
		%======================================================================
		function dens = density(F)
			%%
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
			%%
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
		%>  F = matfaust.rand(50, 50);
		%>  G = matfaust.rand(50, 50);
		%>  [F;G] % equiv. to cat(1,F,G)
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
		%> [F,G] % equiv. to cat(2,F,G)
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
		%> [F;rand(100,50)] % vertical concatenation with auto-conversion of the random matrix to a Faust
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
		%> [F;G;F;G] % it's allowed to concatenate an arbitrary number of Fausts
		%> [F,G,F,G] % as long as the dimensions agree
		%> [F,G;F,G]
		%> @endcode
		%>
		%> @b Errors
		%> - The dimensions of F and G don't agree:
		%>
		%> @code
		%>  F = matfaust.rand(2,51);
		%>  G = matfaust.rand(2,25);
		%>  [F;G]
		%>Error using mexFaustReal
		%>The dimensions of the two Fausts must agree.
		%>
		%> @endcode
		%> - The DIM argument is different from 1 and 2:
		%>
		%> @code
		%>  F = matfaust.rand(2,4);
		%>  G = matfaust.rand(2,4);
		%>  cat(3,F,G)
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
							C = complex(C);
						end
						C = matfaust.Faust(C, call_mex(C, mex_func_name, A.matrix.objectHandle));
					else
						if(isreal(A))
							A = complex(A);
						end
						C = matfaust.Faust(C, call_mex(C, mex_func_name, A.matrix.objectHandle));
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
		%>	F = matfaust.rand(5, 10)
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
			for i=1:min(facs_per_row, numfacts)
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
			suptitle(['Factors of the Faust ' name ])
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
		%> @warning this functions makes a call to Faust.full.
		%>
		%> <p> @b See @b also Faust.pinv, mldivide Matlab built-in.
		%=====================================================================
		function X = mldivide(F,B)
			X = mldivide(full(F),B);
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
		%> @warning this functions makes a call to Faust.full.
		%>
		%> <p> @b See @b also Faust.mldivide, pinv Matlab built-in.
		%=====================================================================
		function X = pinv(F)
			X = pinv(full(F));
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

		%================================================================
		%> Clones the Faust (in a new memory space).
		%===
		%>
		%> @retval Fc: the Faust clone.
		%================================================================
		function Fc = clone(F, varargin)
			%%
			if(nargin > 1)
				dev = varargin{1};
			else
				dev = F.dev;
			end
			if(strcmp(dev, F.dev))
				Fc = matfaust.Faust(F, call_mex(F, 'copy'));
			elseif(strcmp(dev, 'cpu')) % F.device == gpu
				Fc = matfaust.Faust(call_mex(F, 'clone_gpu2cpu'), F.isreal, 'cpu');
			else % dev == gpu
				func_name = 'clone_cpu2gpu';
				if(F.isReal)
					Fc = matfaust.Faust(mexFaustGPUReal(func_name, F.matrix.objectHandle, varargin{:}), F.isreal, 'gpu');
				else
					Fc = matfaust.Faust(mexFaustGPUCplx(func_name, F.matrix.objectHandle, varargin{:}), F.isreal, 'gpu');
				end
			end
		end

		%================================================================
		%> Returns the Faust device ('cpu' or 'gpu').
		%===
		%>
		%> @retval dev: the Faust device.
		%================================================================
		function dev = device(F)
			dev = F.dev
		end

		%================================================================
		%> Returns the Faust class ('single' or 'double').
		%===
		%>
		%> @retval c: the Faust class.
		%================================================================
		function c = class(F)
			if(strcmp(F.dtype, 'float'))
				c = 'single';
			else
				c = 'double';
			end
		end

	end
	methods(Access = public, Hidden = true)
		function set_FM_mul_mode(self, mode)
			set_FM_mul_mode(self.matrix, mode)
		end

		function set_Fv_mul_mode(self, mode)
			set_Fv_mul_mode(self.matrix, mode)
		end

		function H = get_handle(self)
			H = self.matrix.objectHandle;
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

		%================================================================
		%> Helps to call mex functions without the burden of checking the field (real or complex) and the device (CPU or GPU).
		%===
		%================================================================
		function varargout = call_mex(F, func_name, varargin)
			if (strcmp(F.dev, 'cpu'))
				if(F.isReal)
					if(strcmp(F.dtype, 'double'))
						[varargout{1:nargout}] = mexFaustReal(func_name, F.matrix.objectHandle, varargin{:});
					else % float
						[varargout{1:nargout}] = mexFaustRealFloat(func_name, F.matrix.objectHandle, varargin{:});
					end
				else
					[varargout{1:nargout}] = mexFaustCplx(func_name, F.matrix.objectHandle, varargin{:});
				end
			elseif(startsWith(F.dev, 'gpu'))
				if(F.isReal)
					if(strcmp(F.dtype, 'double'))
						[varargout{1:nargout}] = mexFaustGPUReal(func_name, F.matrix.objectHandle, varargin{:});
					else % float
						[varargout{1:nargout}] = mexFaustGPURealFloat(func_name, F.matrix.objectHandle, varargin{:});
					end
				else
					[varargout{1:nargout}] = mexFaustGPUCplx(func_name, F.matrix.objectHandle, varargin{:});
				end
			else
				error('The Faust F has an invalid dev attribute (must be cpu or gpu)')
			end
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

