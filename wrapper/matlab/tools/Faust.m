%class FAUST
% representing a given dense matrix by a product of sparse matrix (i.e faust)
% in order to speed-up multiplication by this matrix,
% matlab wrapper class implemented in C++
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://faust.gforge.inria.fr>

classdef Faust
    properties (SetAccess = private, Hidden = true)
        matrix; % Handle to the underlying C++ class instance
	transpose_flag; % boolean to know if the Faust is transpose or not
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = Faust(varargin)
            % Constructor - build a faust from a cell array of matrix and a scalar (optional)
            %                1st input : 1D cell array of matrix (sparse or dense)
            %                2nd input : (optional) multiplicative scalar
	    %              - or from a filename (mat file) where a faust is stored with save_faust
	   this.matrix = FaustCore(varargin{:});
	   this.transpose_flag = 0;	
	end
	

	
            

        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            % destructor delete the faust
            mexFaust('delete', this.objectHandle);
        end
        
        %% Multiplication faust-vector or faust-matrix
        function y = mtimes(this,x)
            % mtimes - overloading of the matlab multiplication (*) function, compatible with matlab matrix and vector
            y = mexFaust('multiply', this.matrix.objectHandle,x,this.transpose_flag);
        end
        
        
        %% Multiplication by a faust or its transpose
        % if trans = 0 multiplication by faust
        % if trans = 1 multiplication by the transpose of a faust
        function y = mtimes_trans(this,x,trans)
	if ~isreal(trans)
		error('invalid argument trans, must be equal to 0 or 1');
	end

	if (trans ~= 1) && (trans ~= 0)
		error('invalid argument trans, must be equal to 0 or 1');
	end
	   
	isreally_trans=xor(trans,this.transpose_flag);
	y = mexFaust('multiply', this.matrix.objectHandle,x,isreally_trans);
            
        end
        
        %% Evaluate the product of a faust_core
        function y = get_product(this)
            % get_product - compute the dense matrix equivalent to the faust (the product of sparse matrix)
            y=mexFaust('get_product',this.matrix.objectHandle,this.transpose_flag);
	
            
        end

        %% Transpose operator
        function trans=transpose(this)
            %transpose - overloading of the matlab transpose operator (.')
                trans=ctranspose(this); % currently faust is a real matrix, so the complex transposition is the same as the real one 

        end
        
        function trans=ctranspose(this)
	%ctranspose - overloading of the matlab transpose operator (')
                trans=this; % trans and this point share the same C++ underlying object (objectHandle)
                trans.transpose_flag = xor(1,this.transpose_flag); % inverse the transpose flag
        end

	
	
        %% Size
        function varargout = size(this,varargin)
            %size - overload of the matlab size function
            
            
	   nb_input = length(varargin);
            
            
            if (nb_input > 1)
                error('Too many input arguments');
            end
            
            if ((nb_input == 1) && (nargout > 1) | (nargout > 2)) 
                error('Too many output arguments');
            end
            
	
	    Size=mexFaust('size',this.matrix.objectHandle);
	    
	    %% if the faust is tranposed, inverse the dimension		    	
	    if(this.transpose_flag) 
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


    	% end : 
    	% serve as the last index in an indexing expression.  In
    	% that context, end = SIZE(X,k) when used as part of the k-th index.
    	% Examples of this use are, X(3:end) and X(1,1:2:end-1) 
    	function end_dim = end(this,k,n)
		% end - overload of the builtin function end
	
		if (n ~= 2)
			error('invalid slicing : faust is a 2D array i.e matrix');
		end

		end_dim=size(this,k);


	end



	

	%% get_fact : return the id factor of the faust as a dense matrix
        function factor = get_fact(this,id)
		% get_fact : return the id factor of the faust as a dense matrix
		if (~isa(id,'double'))
			error('get_fact second argument (indice) must either be real positive integers or logicals.');
		end

		if (floor(id) ~= id)
			error('get_fact second argument (indice) must either be real positive integers or logicals.');
		end

		if (this.transpose_flag)
			id = get_nb_factor(this)+1-id;
		end

		factor = mexFaust('get_fact',this.matrix.objectHandle,id);

		if (this.transpose_flag)
			factor = factor';
		end
	end


	%% get_nb_factor : return the number of factor of the faust
	function nb_factor = get_nb_factor(this)
		% get_nb_factor : return the number of factor of the faust
		nb_factor = mexFaust('get_nb_factor',this.matrix.objectHandle);
	end

	%% save a faust into a matfile
	function save(this,filename)
		% save a faust into a matfile
		if (~ischar(filename))
			error('second argument must contains a string (a filename)');
		end
		
		nb_fact=get_nb_factor(this);
		
		faust_factors=cell(1,nb_fact);

		for i=1:nb_fact
			faust_factors{i}=get_fact(this,i);
		end
		save(filename,'faust_factors');
		
		
	end

	%% subsref : allows operation such as A(i,j) A(:,j)  A(3:4,2:5) but not A(end,end)
	function submatrix=subsref(this,S)
		% overloading of the slicing method only for reading the value of the coeff
		
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
			invalid(' subsref invalid slicing must have 2 index since this is a 2D-array');
		end

		slicing_row=S.subs{1};
		slicing_col=S.subs{2};
        
        	[dim1 dim2]=size(this);
        
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
		submatrix=mtimes_trans(this,identity,transpose_flag);
		
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


	%% norm : compute the 2-norm of a faust
	function norm_faust=norm(this,varargin)

	    nb_input = length(varargin);
            if (nb_input > 1)
                error('Too many input arguments');
            end
            
            if nb_input == 1
		if varargin{1} ~= 2
                	error('only 2-norm is supported for the faust');
		end
            end
	
	    %% the transpose flag of the faust is ignored because norm(A)==norm(A')
	    norm_faust=mexFaust('norm',this.matrix.objectHandle);
            
	    


	end


	%% nnz : Number of nonzero matrix elements.
	function nz=nnz(this)

	    nz=mexFaust('nnz',this.matrix.objectHandle);
            
	end


	%% nnz : density of the Faust
	function dens=density(this)
	     prod_dim=prod(size(this));			
	     if (prod_dim ~= 0)		
	    	dens=nnz(this)/prod_dim;
             else
		dens = -1;
	     end
	end

	%% nnz : Relative Complexity Gain (inverse of the density)
	function speed_up=RCG(this)	
		dens=density(this);		
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


















