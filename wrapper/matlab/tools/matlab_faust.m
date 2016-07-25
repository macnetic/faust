%class MATLAB_FAUST
% representing a given dense matrix by a product of sparse matrix (i.e faust)
% in order to speed-up multiplication by this matrix,
% matlab wrapper class implemented in C++
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://faust.gforge.inria.fr>

classdef matlab_faust 
    properties (SetAccess = public, Hidden = false)
        matrix; % Handle to the underlying C++ class instance
	transpose_flag; % boolean to know if the matlab_faust is transpose or not
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = matlab_faust(varargin)
            % Constructor - build a faust from a cell array of matrix and a scalar (optional)
            %                1st input : 1D cell array of matrix (sparse or dense)
            %                2nd input : (optional) multiplicative scalar
	    %              - or from a filename (mat file) where a faust is stored with save_faust
	   this.matrix = matlab_faust_core(varargin{:});
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
	    if (this.transpose_flag)
            trans='T';
        else
            trans='N';
	    end
	
            y = mexFaust('multiply', this.matrix.objectHandle,x,trans);
        end
        
        
        %% Multiplication by a faust or its transpose
        % if trans = 'N' multiplication by faust
        % if trans = 'T' multiplication the transpose of a faust
        function y = mtimes_trans(this,x,trans)
	    
	    if xor(strcmp(trans,'T'),this.transpose_flag)
            	trans='T';
        else
            	trans='N';
	    end
	
            y = mexFaust('multiply', this.matrix.objectHandle,x,trans);
            
        end
        
        %% Evaluate the product of a faust_core
        function y = get_product(this)
            % get_product - compute the dense matrix equivalent to the faust (the product of sparse matrix)
            if this.transpose_flag
                trans='T';
            else
                trans='N';
            end	
            y=mexFaust('get_product',this.matrix.objectHandle,trans);
	
            
        end

        %% Transpose operator
        function trans=transpose(this)
            %transpose - overloading of the matlab transpose operator (.')
                trans=ctranspose(this); % currently faust is a real matrix, so the complex transposition is the same as the real one 

        end
        
        function trans=ctranspose(this)
	%ctranspose - overloading of the matlab transpose operator (')
                trans=this; % trans and this point share the same C++ underlying object (objectHandle)
                trans.transpose_flag = mod(this.transpose_flag+1,2); % inverse the transpose flag
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
	function save_faust(this,filename)
		% save a faust into a matfile
		if (~ischar(filename))
			error('second argument must contains a string (a filename)');
		end
		
		nb_fact=get_nb_factor(this);
		
		faust_factors=cell(1,nb_fact);

		for i=1:nb_fact
			faust_factors{i}=get_fact(this,i);
		end
		transpose_flag=this.transpose_flag;
		save(filename,'faust_factors','transpose_flag');
		
		
	end
        
    end
    
end












