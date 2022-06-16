%====================================
%> @brief Constructs a Faust corresponding to the product of log2(n) square factors of size n with butterfly supports and random nonzero coefficients.
%>
%>    The random coefficients are drawn i.i.d. according to a standard Gaussian
%>    (real or complex circular according to the type).
%>
%> @param n: order of the butterfly (must be a power of two).
%> @param 'dev', str: 'gpu or 'cpu' to create the Faust on CPU or GPU ('cpu' by default).
%> @param 'field', str	str is either 'real' or 'complex' (the Faust field).
%>                      The default value is 'real'.
%> @param 'class', str 'double' (by default) or 'single' to select the scalar type used for the Faust generated.
%>
%> @retval F a random butterfly support Faust.
%>
%> @b See also: matfaust.fact.butterfly, matfaust.rand_butterfly, matfaust.dft.
%====================================
function F = rand_butterfly(n, varargin)
    import matfaust.dft
	argc = length(varargin);
	dev = 'cpu';
	class = 'double';
    field = 'real';
	if(argc > 0)
		for i=1:2:argc
			if(argc > i)
				% next arg (value corresponding to the key varargin{i})
				tmparg = varargin{i+1};
			end
			switch(varargin{i})
                case 'field'
					if(argc == i || ~ strcmp(tmparg, 'real') && ~ strcmp(tmparg, 'complex'))
						error('field keyword argument is not followed by a valid value: real, complex.')
					else
						field = tmparg;
					end
				case 'dev'
					if(argc == i || ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error('dev keyword argument is not followed by a valid value: cpu, gpu*.')
					else
						dev = tmparg;
					end
				case 'class'
					if(argc == i || ~ strcmp(tmparg, 'double') && ~ strcmp(tmparg, 'single'))
						error('class keyword argument is not followed by a valid value: single, double.')
					else
						class = tmparg;
					end
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ...
                            ~ startsWith(tmparg, 'gpu') && ~ strcmp(tmparg, 'real') && ...
                            ~ strcmp(tmparg, 'complex') && ~ strcmp(tmparg, 'single') && ~ strcmp(tmparg, 'double'))
						error([ tmparg ' unrecognized argument'])
					end
			end
        end
    end
    if strcmp(field, 'complex') && ~ strcmp(class, 'double')
        error('if field is complex class must be double (complex single is not yet supported)')
    end
    DFT = dft(n);
    % ignore the bitreversal permutation
    B = factors(DFT, 1:numfactors(DFT)-1);
    if n == 2
        B = matfaust.Faust(B);
    end
    RB_factors = cell(1, numfactors(B));
    for k=1:numfactors(B)
        rb = factors(B, k);
        if ~ strcmp(field, 'complex')
            % real rb
            rb = real(rb);
            [I, J, ~] = find(rb);
            for i=1:numel(I)
                rb(I(i), J(i)) = randn(1);
            end
        else
            % complex rb
            [I, J, ~] = find(rb);
            for i=1:numel(I)
                rb(I(i), J(i)) = randn(1) + j * randn(1);
            end
        end
        RB_factors{1, k} = rb;
    end
    F = matfaust.Faust(RB_factors);
end
