% experimental block start
%==========================================================================
%> @brief Factorizes the matrix M according to a butterfly support.
%>
%> @param 'dir', str: the direction of factorization 'right'ward or 'left'ward
%>         (more precisely: at each stage of the factorization the most right factor or
%>          the most left factor is split in two).
%>
%> @retval F the Faust which is an approximate of M according to a butterfly support.
%==========================================================================
function F = butterfly(M, varargin)
        import matfaust.Faust
        nargin = length(varargin);
        dir = 'right';
        if(nargin > 0)
                for i=1:nargin
                        switch(varargin{i})
                                case 'dir'
                                        if(nargin < i+1 || ~ any(strcmp(varargin{i+1}, {'right', 'left'})))
                                                error('keyword argument ''dir'' must be followed by ''left'' or ''right''')
                                        else
                                                dir = varargin{i+1};
                                        end
                        end
                end
        end
        if(strcmp(dir, 'right'))
                dir = 1;
        elseif(strcmp(dir, 'left'))
                dir = 0;
        end
        if(isreal(M))
                core_obj = mexButterflyReal(M, dir);
        else
                core_obj = mexButterflyCplx(M, dir);
        end
        F = Faust(core_obj, isreal(M));
end

% experimental block end
