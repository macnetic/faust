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
%> @brief FAµST class
%>
%> This class represents a given dense matrix by a product of sparse matrices (i.e Faust).
%> The main goal of Faust representation is to speed up operations on that matrix, especially the multiplication. Besides the time optimization, a Faust can reduce the memory space size needed both for storage and loading.
%>
%> Although the sparse matrices are more interesting for optimization it's not forbidden to define a Faust as a product of dense matrices or a mix up of dense and sparse matrices.
%>
%> The matrices composing the Faust product, also called the factors, are defined on complex or real fields. Hence a Faust can be a complex Faust or a real Faust.
%>
% ======================================================================
classdef Faust
	properties (SetAccess = private, Hidden = true)
		matrix; % Handle to the FaustCore class instance
		isReal;
	end
	properties (Constant)
		%> Identifies a complex Faust.
		COMPLEX=3
		%> Identifies a real Faust.
		REAL=4
		% Constants to identify kind of factors to generate through Faust.rand()
		%> Designates a dense factor matrix
		DENSE=0
		%> Designates a dense factor matrix
		SPARSE=1
		%> Means DENSE or SPARSE
		MIXTE=2
	end
	methods
		%======================================================================
		%> @brief Creates a Faust from a list of factors or alternatively from a file.
		%>
		%> Another easy way to create a Faust is to call the static method Faust.rand().
		%>
		%> @param factors (varargin{1}) the 1D cell array of factors to initialize the Faust with.
		%> <br/> The factors must respect the dimensions needed for the product to be defined (for i=1 to size(factors,2), size(factors{i},2) == size(factors{i+1},1)).
		%> <br/> The factors can be sparse or dense matrices.
		%> @param filepath (varargin{1}) the file from which a Faust is created.<br/>
		%>								The format is Matlab version 5 (.mat extension).<br/>
		%>								The file must have been saved before with Faust.save().
		%> @param lambda (optional varargin{2}) multiplicative scalar applied to the factor product before to set the Faust with.
		%>
		%> @b Examples
		%> @code
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
		%>	lambda = 2
		%>	G = Faust(factors, lambda) % G == lambda*F
		%>
		%>	save(F, 'F.mat')
		%>	% define a Faust from file
		%>	H = Faust('F.mat')
		%>	I = Faust('F.mat', lambda) % I == lambda*H
		%>
		%> @endcode
		%>
		%>
		%> <p>@b See @b also Faust.delete, Faust.save, Faust.rand</p>
		%>
		%======================================================================
		function F = Faust(varargin)
			%% FAUST Constructor - Creates a Faust from various types of input.
			%
			% Examples :
			%
			% F = Faust(factors,lambda)
			% - factor : 1D cell array of matrix (sparse or
			% dense) representing the factor of the Faust.
			% - lambda : (optional) multiplicative scalar.
			%
			% F = Faust(filepath, lambda)
			% - filepath: the file where a Faust was stored with Faust.save() (in matlab format version 5).
			% - lambda: (optional) multiplicative scalar.
			err_msg = 'Faust() error: the arguments are not valid.';
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
					error([err_msg ' The second argument lambda must be a scalar.'])
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
				F = Faust(faust_factors, varargin{2:end});
			elseif(isa(varargin{1}, 'Faust'))
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
		%>
		%>
		%> @param F the Faust to delete.
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.REAL, 2, 5, 50, 100, .5)
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
		%> @brief Multiplies the Faust F to the full storage matrix A.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%>
		%> @param F the Faust object.
		%> @param A The matrix to multiply (full storage matrix).
		%>
		%> @retval B The multiplication result (full storage matrix).
		%>
		%> @b Example
		%> @code
		%>   F = Faust.rand(Faust.MIXTE, Faust.REAL, 2, 5, 50, 100, .5)
		%>   A = rand(size(F,2), 50)
		%>   B = F*A
		%> % is equivalent to B = mtimes(F, A)
		%> @endcode
		%>
		%> <p>@b See @b also mtimes_trans.
		%>
		%======================================================================
		function B = mtimes(F,A)
			%% MTIMES * Faust Multiplication (overloaded Matlab built-in function).
			%
			% B=mtimes(F,A) is called for syntax 'B=F*A', when F is a Faust matrix and A a full
			% storage matrix, B is also a full matrix storage.
			%
			% See also mtimes_trans
			B = mtimes_trans(F, A, 0)
		end


		%======================================================================
		%> @brief Multiplies the Faust or its transpose to the A full storage matrix.
		%>
		%>
		%> @param A the matrix to multiply (full storage matrix).
		%> @param trans equals 1 to calculate C=F'*A
		%> 				or 0 to calculate C=F*A.
		%>
		%> @retval C The multiplication result (full storage matrix).
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

			if (F.isReal)
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
		%> @b Warning: this function costs F.get_nb_factor()-1 matrix multiplications.
		%>
		%> @retval A the full storage matrix resulting from the Faust.
		%>
		%======================================================================
		function A = full(F)
			%% FULL  Convert Faust matrix to full matrix (overloaded Matlab
			% built-in function).
			%
			% A=full(F) converts a Faust matrix F to full storage matrix A.
			if (F.isReal)
				A=mexFaustReal('full',F.matrix.objectHandle);
			else
				A=mexFaustCplx('full',F.matrix.objectHandle);
			end

		end

		%======================================================================
		%> @brief Indicates if F is a real Faust or a complex Faust.
		%>
		%> This function overloads a Matlab built-in function.
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
		%> @brief Transposes the Faust F.
		%>
		%> This function overloads a Matlab built-in function/operator.
		%>
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
				F_trans = Faust(F, mexFaustReal('transpose', F.matrix.objectHandle));
			else
				F_trans = Faust(F, mexFaustCplx('transpose', F.matrix.objectHandle));
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
		%>	F = Faust.rand(Faust.MIXTE, Faust.COMPLEX, 2, 5, 50, 100, .5)
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
		function F_ctrans=ctranspose(F)
			%% CTRANSPOSE ' Complex conjugate transposed Faust (overloaded Matlab built-in function).
			%
			% F_trans = ctranspose(F) is called for syntax F' (complex conjugate transpose) when F is a Faust.
			%
			%
			% See also transpose.
			if (isreal(F))
				F_ctrans=transpose(F);
			else
				F_ctrans = Faust(F, mexFaustCplx('ctranspose', F.matrix.objectHandle));
			end
		end

		%======================================================================
		%> @brief Returns the conjugate of the Faust.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @param F the Faust object.
		%>
		%> @retval F_conj = conj(F)
		%> <br/> if F is a complex Faust, the value returned is F_conj == REAL(F) - i*IMAG(F).
		%> <br/> If F is a real Faust then F_conj == F.
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.COMPLEX, 2, 5, 50, 100, .5)
		%>	F_conj = conj(F)
		%> @endcode
		%>
		%======================================================================
		function F_conj=conj(F)
			%% CONJ ' Complex conjugate Faust (overloaded Matlab built-in function).
			%
			%  F_conj = conj(F) For a complex F, F_conj == REAL(F) - i*IMAG(F)
			%
			%
			if (F.isReal)
				F_conj = Faust(F, mexFaustReal('conj', F.matrix.objectHandle));
			else
				F_conj = Faust(F, mexFaustCplx('conj', F.matrix.objectHandle));
			end
		end


		%======================================================================
		%> @brief Gives the size of the Faust.
		%>
		%> @param F the Faust object.
		%> @param varargin can be missing or specifying the index of the dimension to get the size of.
		%>
		%> @retval [NROWS,NCOLS] = size(F)
		%> @retval N = size(F,DIM) with N being the size of DIM-th dimension of F.
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.COMPLEX, 2, 5, 50, 100, .5)
		%>	[nlines, ncols] = size(F)
		%>	nlines = size(F, 1)
		%>	ncols = size(F, 2)
		%> @endcode
		%>
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
		%> @brief Serves as the last index when slicing or indexing a Faust.
		%>
		%> This function overloads a Matlab built-in function.
		%>
		%> @b Example
		%> @code
		%>	% in a matlab terminal
		%>
		%>	>> F = Faust.rand(Faust.MIXTE, Faust.REAL, 4, 5, 3, 3, .9);
		%>	>> full(F)
		%>	ans =
		%>
		%>		-0.1006   -0.2041   -0.1878
		%>		-0.1382    0.0400   -0.0954
		%>		-0.1345   -0.1223   -0.1667
		%>
		%>	>> F(2:end,1)
		%>	ans =
		%>
		%>		-0.1382
		%>		-0.1345
		%>
		%>	>> F(1,1:2:end)
		%>	ans =
		%>
		%>	-0.1006   -0.1878
		%>
		%>	>> F(1,1:end-1)
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
		%> @brief Returns the i-th factor of F.
		%>
		%> @param F the Faust object.
		%> @param i the factor index.
		%>
		%> @retval factor the i-th factor as a full storage matrix.
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.COMPLEX, 2, 5, 50, 100, .5)
		%>	f1 = get_fact(F,1)
		%> @endcode
		%> <p>@b See @b also Faust.get_nb_factor
		%=====================================================================
		function factor = get_fact(F,i)
			%% GET_FACT Ith factor of the Faust.
			%
			% A=get_fact(F,i) return the i factor A of the Faust F as a full storage matrix.
			%
			% Example of use :
			% A=get_fact(F,1) returns the 1st factor of the Faust F.
			% A=get_fact(F,4) returns the 4th factor of the Faust F.
			%
			% See also get_nb_factor.

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


		%===========================================================================================
		%> @brief Gives the number of factors of F.
		%>
		%> @param F the Faust object.
		%>
		%> @retval num_factors the number of factors.
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.COMPLEX, 2, 5, 50, 100, .5)
		%>	nf = get_nb_factor(F)
		%> @endcode
		%>
		%> <p>@b See @b also Faust.get_fact.
		%===========================================================================================
		function num_factors = get_nb_factor(F)
			%% GET_NB_FACTOR Number of factor of the Faust.
			%
			% num_factors = get_num_factors(F) return the number of factor of the
			% Faust F.
			%
			% See also get_fact.
			if (F.isReal)
				num_factors = mexFaustReal('get_nb_factor', F.matrix.objectHandle);
			else
				num_factors = mexFaustCplx('get_nb_factor', F.matrix.objectHandle);
			end
		end

		%===========================================================================================
		%> @brief Saves the Faust F into file respecting the Matlab format version 5 (.mat file).
		%>
		%> @param F the Faust object.
		%> @param filepath the path for saving the Faust (filename should ends with .mat).
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.COMPLEX, 2, 5, 50, 100, .5)
		%>	save(F, 'F.mat')
		%>	G = Faust('F.mat')
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust.
		%===========================================================================================
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

		%===========================================================================================
		%> @brief Gets a submatrix of the full matrix of F.
		%>
		%> This function is a Matlab built-in overload.
		%>
		%> @b WARNING: this function costs as much as Faust.mtimes.
		%>
		%> @param F the Faust object.
		%> @param S the subscript defining the submatrix (see examples below).
		%>
		%> @retval submatrix the full submatrix requested.
		%>
		%> @b Example
        %> @code
        %>		F = Faust.rand(Faust.MIXTE, Faust.REAL, 2, 5, 50, 100, .5)
        %>		i = randi(min(size(F)), 1, 2)
		%>	i1 = i(1);i2 = i(2)
		%>
		%>	F(i1,i2) % element at line i1, column i2
		%>
		%>	F(:,i2) % full column i2
		%>
		%>	F(3:4,2:5) % submatrix from line 3 to line 4, each line containing only elements from column 2 to 5.
		%>
		%>	F(1:end,5:end-1)  % submatrix from line 1 to end line, each line containing only elements from column 5 to column before the last one.
        %> @endcode
		%>
		%> <p>@b See @b also Faust.end.
		%===========================================================================================
		function submatrix=subsref(F,S)
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
				error(' subsref invalid structure S, S.type must be a char, S.subs must be a cell-array');
			end

			if ~strcmp(S.type,'()')
				error(' subsref is only overloaded for () operator');
			end

			if (length(S.subs) ~=2)
				invalid(' subsref invalid slicing must have 2 index since F is a 2D-array');
			end

			slicing_row=S.subs{1};
			slicing_col=S.subs{2};

			[dim1 dim2]=size(F);

			if ischar(slicing_row)
				nb_row_selected = dim1;
			else
				nb_row_selected = length(slicing_row);
			end

			if ischar(slicing_col)
				nb_col_selected = dim2;
			else
				nb_col_selected = length(slicing_col);
			end

			% evaluate the complexity of getting the coeff by doing
			%  A*Identity or A'*Identity
			transpose_evaluation =  (nb_col_selected > nb_row_selected);
			if transpose_evaluation
				identity=eye(dim1);
				transpose_flag=1;

				% switch the 2 different slicing
				tmp=slicing_row;
				slicing_row=slicing_col;
				slicing_col=tmp;

			else
				identity=eye(dim2);
				transpose_flag=0;
			end

			% selects the column of the identity, if slicing_col is a char, all
			% the column are selected
			if ~ischar(slicing_col)
				identity=identity(:,slicing_col);
			end

			% perform A*identity or A'*identity
			submatrix=mtimes_trans(F,identity,transpose_flag);

			% selects the row of the submatrix, if slicing_row is a char, all
			% the row are selected
			if ~ischar(slicing_row)
				submatrix=submatrix(slicing_row,:);
			end

			% transpose if needed
			if transpose_evaluation
				submatrix=submatrix';
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
		%> @param F the Faust object.
		%>
		%> @b WARNING: currently a bug is affecting this function. When the Faust is transposed the dimensions are inverted in the display (like the Faust hasn't been transposed).
		%>
		%> @b Example
		%> @code
		%>	% in a matlab terminal
		%>	>> F = Faust.rand(Faust.MIXTE, Faust.REAL, 1, 2, 50, 100, .5)
		%>	>> disp(F)
		%>	Faust of size : 55x73, nb factor 2, RCG 0.947157,nnz 4239
		%>	- FACTOR 0type : SPARSE size 55x66, density 0.501377, nnz 1820
		%>	- FACTOR 1type : DENSE size 66x73, density 0.502076, nnz 2419
		%>	>> F
		%>	Faust of size : 55x73, nb factor 2, RCG 0.947157,nnz 4239
		%>	- FACTOR 0type : SPARSE size 55x66, density 0.501377, nnz 1820
		%>	- FACTOR 1type : DENSE size 66x73, density 0.502076, nnz 2419
		%>
		%> @endcode
		%>
		%> <p>@b See @b also Faust.nnz, Faust.RCG, Faust.size, Faust.get_fact
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
		%> @b Note: the norm of F is equal to the norm of its dense matrix.
		%>
		%> @b WARNING: this function costs at least as much as Faust.mtimes.
		%>
		%>
		%> @param F the Faust object.
		%> @param ord (optional) the norm order. Respectively 1 or 2 for the 1-norm and 2-norm or 'fro' for the Frobenius norm (by default the 2-norm is computed).
		%>
		%>
		%> @retval norm_Faust the norm (real).
		%>
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.REAL, 1, 2, 50, 100, .5)
		%>	norm(F)
		%>
		%>	ans =
		%>
		%>	   34.2995
		%>	norm(F,2)
		%>
		%>	ans =
		%>
		%>	   34.2995
		%>
		%> @endcode
		%>
		%>
		%>
		%======================================================================
		function norm_Faust=norm(F,varargin)
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
						norm_Faust=mexFaustReal('normfro',F.matrix.objectHandle);
					else
						norm_Faust=mexFaustCplx('normfro',F.matrix.objectHandle);
					end
					return
				end
				if varargin{1} ~= 2 && varargin{1} ~= 1
					error('only L1, L2 or Frobenius norms are supported for the Faust');
				end
				ord = varargin{1};
			end

			if (F.isReal)
				norm_Faust=mexFaustReal('norm',F.matrix.objectHandle, ord);
			else
				norm_Faust=mexFaustCplx('norm',F.matrix.objectHandle, ord);
			end

		end

		%===========================================================================================
		%> @brief Gives the total number of non-zero elements in F's factors.
		%>
		%> The function sums together the number of non-zeros elements of
		%> each factor and returns the result. Note that in fact the sum is
		%>                        computed at Faust creation time and kept in cache.
		%>
		%> @param F the Faust object.
		%>
		%> @retval nz The number of non-zeros.
		%>
		%> <p>@b See @b also Faust.RCG, Faust.density.
		%===========================================================================================
		function nz=nnz(F)
			%% NNZ Number of nonzero elements in a Faust (overloaded Matlab built-in function).
			%
			% nz = nnz(F) is the number of nonzero elements in the Faust F.
			%
			% See also density, RCG.
			if (F.isReal)
				nz=mexFaustReal('nnz',F.matrix.objectHandle);
			else
				nz=mexFaustCplx('nnz',F.matrix.objectHandle);
			end
		end

		%======================================================================
		%> @brief Calculates the density of F, that is, the number of non-zeros
		%> in factors over the total number of elements in dense matrix of F
		%> (which is equal to size(F, 1)*size(F, 2)).
		%>
		%> @b NOTE: this definition of density allows the value to be greater than 1.
		%>
		%> @param F the Faust object.
		%>
		%> @retval dens density of F.
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.REAL, 2, 5, 50, 100, .5)
		%>	dens = density(F)
		%> @endcode
		%>
		%> <p/>@b See @b also Faust.nnz, Faust.RCG
		%======================================================================
		function dens=density(F)
			%% DENSITY Density of the Faust.
			%
			% dens = density(F) when F is a Faust returns the
			% percentage of nonzero elements of F,
			% dens is a number between 0 and 1.
			% In some degenerated case, dens can be greater than 1.
			% If the Faust is empty, return -1.
			%
			% See also RCG, nnz.

			prod_dim=prod(size(F));
			if (prod_dim ~= 0)
				dens=nnz(F)/prod_dim;
			else
				dens = -1;
			end
		end

		%===========================================================================================
		%> @brief Computes the Relative Complexity Gain (inverse of Faust.density).
		%>
		%> RCG is the theoretical gain brought by Faust representation relatively to its dense
		%> matrix equivalent. That gain applies both for storage and multiplication computation
		%> time.
		%>
		%>
		%> @param F	the Faust object.
		%>
		%>
		%> @retval speed_up =  the RCG value (real). If the density is zero it will be Inf. If the density is negative it will be -1.
		%>
		%> <p>@b See @b also Faust.density, Faust.nnz.
		%===========================================================================================
		function speed_up=RCG(F)
			%% RCG Relative Complexity Gain (inverse of the density)
			%
			% speed_up =  RCG(F) when F is Faust, returns the
			% inverse of density of the Faust (i.e the theoretical gain
			% both for storage and multiplication computation time between the Faust and its full storage
			% equivalent full(F)).
			%
			% See also density, nnz.

			dens=density(F);
			if (dens > 0)
				speed_up=1/dens;
			else
				if (dens == 0)
					speed_up=Inf;
				else
					speed_up = -1;
				end
			end
		end


	end
	methods(Static)
		%===========================================================================================
		%> @brief Generates a random Faust.
		%>
		%> @param faust_type must be one of Faust.DENSE, Faust.SPARSE or Faust.MIXTE (the latter is for allowing generation of dense and sparse factors in the same Faust).
		%> @param field	must be Faust.REAL or Faust.COMPLEX.
		%> @param min_num_factors the minimal number of factors generated.
		%> @param max_num_factors the maximal number of factors generated.
		%> @param min_dim_size	the minimal size of column and row dimensions of the Faust generated.
		%> @param max_dim_size	the maximal size of column and row dimensions of the Faust generated.
		%> @param density	the approximate density of factors generated.
		%>
		%>
		%> @retval F the random Faust.
		%>
		%> @b Example
		%> @code
		%>	F = Faust.rand(Faust.MIXTE, Faust.REAL, 2, 5, 50, 100, .5)
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust.
		%===========================================================================================
		function F=rand(faust_type, field, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density)
			if(faust_type ~= Faust.SPARSE && faust_type ~= Faust.DENSE && faust_type ~= Faust.MIXTE)
				error('Faust.rand(): error: faust_type must be among Faust.SPARSE/DENSE/MIXTE.')
			end
			if(field == Faust.COMPLEX)
				F = Faust(mexFaustCplx('rand', faust_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density), false)
			elseif(field == Faust.REAL)
				F = Faust(mexFaustReal('rand', faust_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density), true)
			else
				error('Faust.rand() error: field must be Faust.REAL or Faust.COMPLEX.')
			end
		end
	end
end

