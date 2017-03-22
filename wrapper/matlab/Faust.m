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


classdef Faust
    properties (SetAccess = private, Hidden = true)
        matrix; % Handle to the underlying C++ class instance
	transpose_flag; % boolean to know if the Faust is transpose or not
	isReal;
    end
    methods

        function F = Faust(varargin)
	%% FAUST Constructor - build a Faust from various type of input.
	%
	% Example of use :
	%
	% F = Faust(factors,lambda) 
	% -factor : 1D cell array of matrix (sparse or
	% dense) representing the factor of the Faust
	% -lambda : (optional) multiplicative scalar
	%                 
	% F = Faust(filename) 
	% filename : a filename (mat file) where a Faust is stored with save_Faust
        

	    if (nargin == 1) && ischar(varargin{1})
		filename=varargin{1};
		load(filename);
		if (~exist('faust_factors','var') )
			error('FaustCore : invalid file');
		end
		F=Faust(faust_factors);
	    else
		%check if the factors are real or complex, at least one complex factor means a complex faust 
		factors=varargin{1};
		isRealFlag = 1;
	
		for i=1:length(factors)
			if (~isreal(factors{i}))
				isRealFlag = 0;
			end
		end

		F.matrix = FaustCore(varargin{:});
	   	F.transpose_flag = 0;
	   	F.isReal = isRealFlag;
		
	    end

	   
	   	
		
	 	
	end
	

	
            


        function delete(F)
	%% DELETE Destructor delete the Faust.
	% delete(F)
	%
	% See also Faust
            if (F.isReal)
            	%mexFaustReal('delete', F.matrix.objectHandle);
		delete(F.matrix);
	    end
        end
        

        function C = mtimes(F,A)
	%% MTIMES * Faust Multiplication (overloaded Matlab built-in function).
	%
	% C=mtimes(F,A) is called for syntax 'C=F*A', when F is a Faust matrix and A a full
	% storage matrix, C is also a full matrix storage.
	% 
	% See also mtimes_trans
            if (F.isReal)	
            	C = mexFaustReal('multiply', F.matrix.objectHandle,A,F.transpose_flag);
	    else
		C = mexFaustCplx('multiply', F.matrix.objectHandle,A,F.transpose_flag);
	    end
        end
        
        

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
	   
	isreally_trans=xor(trans,F.transpose_flag);
	if (F.isReal)
		C = mexFaustReal('multiply', F.matrix.objectHandle,A,isreally_trans);
	else
		C = mexFaustCplx('multiply', F.matrix.objectHandle,A,isreally_trans);
	end
            
        end
        

        function A = full(F)
	%% FULL  Convert Faust matrix to full matrix (overloaded Matlab
	% built-in function).
	%
	% A=full(F) converts a Faust matrix F to full storage matrix A.
        if (F.isReal)
        	A=mexFaustReal('full',F.matrix.objectHandle,F.transpose_flag);
	else
		A=mexFaustCplx('full',F.matrix.objectHandle,F.transpose_flag);
	end
            
        end


        function F_trans=transpose(F)
	%% TRANSPOSE .' Non-conjugate transposed Faust (overloaded Matlab built-in function).
	%
	% F_trans = transpose(F) is called for the syntax F.' when F is Faust.
	%
	% WARNING : currently Faust is a real matrix, so the conjugate transposition is the same as the real one
	%
	% See also ctranspose.
                F_trans=ctranspose(F); 
        
	end
        
        function F_ctrans=ctranspose(F)
	%% CTRANSPOSE ' Complex conjugate transposed Faust (overloaded Matlab built-in function).
	%
	% F_trans = ctranspose(F) is called for syntax F' (complex conjugate transpose) when F is a Faust.
	%
	% WARNING : currently Faust is a real matrix, so the conjugate transposition is the same as the real one
	%
	% See also transpose.
        
                F_ctrans=F; % trans and F point share the same C++ underlying object (objectHandle)
                F_ctrans.transpose_flag = xor(1,F.transpose_flag); % inverse the transpose flag
        end

	
	

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

	    % if the Faust is tranposed, inverse the dimension		    	
	    if(F.transpose_flag) 
	    	Size = Size*[0,1;1,0];
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



	

        function factor = get_fact(F,id)
	%% GET_FACT Ith factor of the Faust.
	%
	% A=get_fact(F,id) return the id factor A of the Faust F as a full storage matrix.
	% 
	% Example of use :
	% A=get_fact(F,1) returns the 1st factor of the Faust F.
	% A=get_fact(F,4) returns the 4th factor of the Faust F.
	%
	% See also get_nb_factor.
        
		if (~isa(id,'double'))
			error('get_fact second argument (indice) must either be real positive integers or logicals.');
		end

		if (floor(id) ~= id)
			error('get_fact second argument (indice) must either be real positive integers or logicals.');
		end

		if (F.transpose_flag)
			id = get_nb_factor(F)+1-id;
		end
		if (F.isReal)
			factor = mexFaustReal('get_fact',F.matrix.objectHandle,id);
		else
			factor = mexFaustCplx('get_fact',F.matrix.objectHandle,id);
		end
		if (F.transpose_flag)
			factor = factor';
		end
	end



	function nb_factor = get_nb_factor(F)
	%% GET_NB_FACTOR Number of factor of the Faust.
	%
	% nb_factor = get_nb_factor(F) return the number of factor of the
	% Faust F.
	%
	% See also get_fact.
    		if (F.isReal)
			nb_factor = mexFaustReal('get_nb_factor',F.matrix.objectHandle);
		else
			nb_factor = mexFaustReal('get_nb_factor',F.matrix.objectHandle);
		end
    end

    
	function save(F,filename)
	%% SAVE Save a Faust into a matfile.
	%
	% save(F,filename) save the Faust F into the .mat file specified by
	% filename.
    	
    	
    
		if (~ischar(filename))
			error('second argument must contains a string (a filename)');
		end
		
		nb_fact=get_nb_factor(F);
		
		faust_factors=cell(1,nb_fact);

		for i=1:nb_fact
			faust_factors{i}=get_fact(F,i);
		end
		save(filename,'faust_factors');
		
		
	end


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

    function F = subsasgn(F,S,B)
    %% SUBSASGN (WARNING not implemented) (overloaded Matlab built-in function)
    %
    % This function is no available for Faust class,
    % this function just throw an error
    %
    % F(i,j)=1, F(2:5,3:5)=zeros(4,3) will throw 
    % a Matlab error with this message :
    % 'function not implemented for Faust class'
          error('function not implemented for Faust class');   
        
    end
    
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
    

	function norm_Faust=norm(F,varargin)
	%% NORM Faust norm (overloaded Matlab built-in function).
	%
	% norm(F,2) when F is Faust returns the 2-norm of F
	% norm(F) is the same as norm(F)
	%
	% WARNING : norm(F,typenorm) is only supported when typenorm equals 2
    
	    nb_input = length(varargin);
            if (nb_input > 1)
                error('Too many input arguments');
            end
            
            if nb_input == 1
		if varargin{1} ~= 2
                	error('only 2-norm is supported for the Faust');
		end
            end
	
	    % the transpose flag of the Faust is ignored because norm(A)==norm(A')
	    if (F.isReal)	
	    	norm_Faust=mexFaustReal('norm',F.matrix.objectHandle);
	    else
		norm_Faust=mexFaustCplx('norm',F.matrix.objectHandle);
            end
	    


	end


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
    
end




















