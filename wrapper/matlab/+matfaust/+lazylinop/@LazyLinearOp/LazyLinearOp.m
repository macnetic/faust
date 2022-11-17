% ======================================================================
%> @brief This class implements a lazy linear operator.
%>
%> The evaluation of any defined operation is delayed until proceeding to a multiplication
%> by a dense matrix/vector, a call of LazyLinearOp.toarray.
%>
%> To instantiate a LazyLinearOp look at matfaust.lazylinop.aslazylinearoperator or matfaust.lazylinop.LazyLinearOperator to instantiate from matmat/matvec functions.
%>
%> @warning This code is in a beta status.
% ======================================================================
classdef LazyLinearOp < handle % needed to use references on objects 
    properties (SetAccess = protected, Hidden = true)
        lambdas;
        shape;
        dtype;
        root_obj;
    end
    properties(SetAccess = protected, Constant, Hidden = true)
        MUL = 1;
        T = 2;
        H = 3;
        I = 4; % indeing / slicing
        STRS = {'MUL', 'T', 'H', 'I'}
    end
    methods

        %======================================================================
        %> @brief Constructor. Not meant to be used directly.
        %>
        %> @param lambdas: starting operations.
        %> @param shape: the initial shape of the operator.
        %> @param root_obj: (optional) the initial object the operator is based on.
        %>
        %> <p>@b See @b also matfaust.lazylinop.aslazylinearoperator
        %======================================================================
        function L = LazyLinearOp(lambdas, shape, varargin)
            L.lambdas = lambdas;
            L.shape = shape;

            p = inputParser;
            validDtype = @(dtype) any(strcmp(dtype, {'complex', 'double', 'single', 'undefined'}));
            addOptional(p, 'dtype', 'undefined', validDtype)
            addOptional(p, 'root_obj', 'none')

            parse(p, varargin{:})
            L.dtype = p.Results.dtype;
            L.root_obj = p.Results.root_obj;

        end

        function check_lambdas(L)
            import matfaust.lazylinop.LazyLinearOp
            if ~ iscell(L.lambdas)
                error('lambdas must be a cell array')
            end
            if ~ length(L.lambdas) == 4
                error('lambdas must be 4 cells long')
            end
            for i=1:4
                if ~ isa(L.lambdas{i}, 'function_handle')
                    error([LazyLinearOp.STRS{i} ' lambda is not a proper function handle'])
                end
            end
        end

        function L = index_lambda(L, S)
            import matfaust.lazylinop.LazyLinearOp

            if (~isfield(S,'type')) | (~isfield(S,'subs'))
                error(' subsref invalid structure S missing field type or subs');
            end

            if (~ischar(S.type)) | (~iscell(S.subs))
                error(' subsref invalid structure S, S.type must be a character array, S.subs must be a cell array');
            end

            if ~ strcmp(S.type,'()')
                error(' subsref is only overloaded for () operator');
            end

            if (length(S.subs) ~=2)
                error(' subsref invalid slicing must have 2 index since L is a 2D-array');
            end

            irows = speye(size(L, 1));
            irowsS.type = '()';
            irowsS.subs = {S.subs{1}, ':'};
            icols = speye(size(L, 2));
            icolsS.type = '()';
            icolsS.subs = {':', S.subs{2}};
            L = LazyLinearOp.create_from_op(subsref(irows, irowsS)) * L * LazyLinearOp.create_from_op(subsref(icols, icolsS));
        end

        %======================================================================
        %> @brief size of L.
        %>
        %> @param dim: the dimension (1 for rows, 2 for columns) to get the size of.
        %>
        %> @retval s: by default, an array of two numbers (the number of rows and the number of columns) or the size of only one dimension if dim argument is used.
        %>
        %======================================================================
		function s = size(L, varargin)
			if length(varargin) == 0
				s = L.shape;
			else
				s = L.shape(varargin{1});
			end
		end

        %=============================================================
        %> @brief Returns the LazyLinearOp transpose.
        %=============================================================
        function LT = transpose(L)
            check_meth(L, 'transpose')
            LT = L.lambdas{matfaust.lazylinop.LazyLinearOp.T}();
        end

        %=============================================================
        %> @brief Returns the LazyLinearOp ctranspose.
        %=============================================================
        function LH = ctranspose(L)
            check_meth(L, 'ctranspose')
            LH = L.lambdas{matfaust.lazylinop.LazyLinearOp.H}();
        end

        %=============================================================
        %> @brief Returns the LazyLinearOp conjugate.
        %=============================================================
        function LC = conj(L)
            check_meth(L, 'conj')
            LC = (L').';
        end

        %=============================================================
        %> @brief Returns the LazyLinearOp for indexing.
        %>
        %> @b See @b also: subsref matlab built-in.
        %=============================================================
        function LS = subsref(L, S)

			if (~isfield(S,'type')) | (~isfield(S,'subs'))
				error(' subsref invalid structure S missing field type or subs');
			end

			if (~ischar(S.type)) | (~iscell(S.subs))
				error(' subsref invalid structure S, S.type must be a character array, S.subs must be a cell array');
			end

			if ~strcmp(S.type,'()')
				error(' subsref is only overloaded for () operator');
			end

			if (length(S.subs) ~=2)
				error(' subsref invalid slicing must have 2 index since L is a 2D-array');
			end

            LS = L.lambdas{matfaust.lazylinop.LazyLinearOp.I}(S)
        end


        %=============================================================
        %> @brief Returns the LazyLinearOp for the multiplication self * op
        %> or if op is a full matrix it returns the full matrix (self * op).
        %>
        %> @param op: an object compatible with self for this binary operation.
        %=============================================================
		function LM = mtimes(L, A)
			import matfaust.lazylinop.LazyLinearOp
			if LazyLinearOp.isLazyLinearOp(A)
				if isscalar(L)
					LM = mtimes(A, L);
					return;
				elseif ~ LazyLinearOp.isLazyLinearOp(L)
					LM = ctranspose(mtimes(A.', L.'));
					return;
				end
			end
			check_meth(L, 'mtimes');
			op_is_scalar = all(size(A) == [1, 1]);
			if ~ op_is_scalar && ~ all(size(L, 2) == size(A, 1))
				error('Dimensions must agree')
			end
			if op_is_scalar
				new_size = size(L);
			else
				new_size = [size(L, 1), size(A, 2)];
			end

			function l = mul_index_lambda(L, A, S)
				% L and A must be LazyLinearOp
				import matfaust.lazylinop.LazyLinearOp
				Sr.type = '()';
				Sr.subs = {S.subs{1}, ':'};
				Sc.type = '()';
				Sc.subs = {':', S.subs{2}};
				L.lambdas{L.I}(Sr) * A.lambdas{L.I}(Sc);
			end

			if ~ LazyLinearOp.isLazyLinearOp(A) && ismatrix(A) && isnumeric(A) && any(size(A) ~= [1, 1])
				% A is a dense matrix that is not limited to one element
				LM = L.lambdas{L.MUL}(A);
			else
				if isscalar(A)
					LM_size = size(L);
				else
					if ~ LazyLinearOp.isLazyLinearOp(A)
						A = LazyLinearOp.create_from_op(A);
					end
					LM_size = [size(L, 1), size(A, 2)] 
				end

				lambdas = {@(o) L * (A * o), ... %MUL
					@() A.' * L.', ... % T 
					@() A' * L', ... % H 
					@(S) mul_index_lambda(L, A, S)% I
				};
				LM = LazyLinearOp(lambdas, LM_size);
			end
		end


        %=============================================================
        %> @brief Returns the LazyLinearOp L + op.
        %>
        %> @param op: an object compatible with self for this binary operation.
        %=============================================================
        function LP = plus(L, op)
            import matfaust.lazylinop.LazyLinearOp
            check_meth(L, 'plus')
            if any(size(L) ~= size(op))
                error('dimensions must agree')
            end
            if ~ all(size(op) == [1, 1]) && ~ all(size(L) == size(op))
                error('Dimensions must agree')
            end
            if ~ LazyLinearOp.isLazyLinearOp(op)
                op = LazyLinearOp.create_from_op(op);
            end
            lambdas = {@(o) L * o + op * o, ... %MUL
                @() L.' + op.', ... % T
                @() L' * op', ... % H
                @(S) index_lambda(L, S) + index_lambda(op, S) % I
            };
            LP = LazyLinearOp(lambdas, [size(L, 1), size(op, 2)]);
        end

        %=============================================================
        %> @brief Returns the LazyLinearOp (+ L).
        %=============================================================
        function LUP = uplus(L)
            LUP = L;
        end

        %=============================================================
        %> @brief Returns the LazyLinearOp L - op.
        %>
        %> @param op: an object compatible with self for this binary operation.
        %=============================================================
        function LUP = minus(L, op)
            import matfaust.lazylinop.LazyLinearOp
            check_meth(L, 'minus')
            if any(size(L) ~= size(op))
                error('dimensions must agree')
            end
            if ~ all(size(op) == [1, 1]) && ~ all(size(L) == size(op))
                error('Dimensions must agree')
            end
            if ~ LazyLinearOp.isLazyLinearOp(op)
                op = LazyLinearOp.create_from_op(op);
            end
            lambdas = {@(o) L * o - op * o, ... %MUL
                @() L.' - op.', ... % T
                @() L' * op', ... % H
                @(S) index_lambda(L, S) - index_lambda(op, S) % I
            };
            LUP = LazyLinearOp(lambdas, [size(L, 1), size(op, 2)]);
        end

        %=============================================================
        %> @brief Returns the full matrix resulting from L.
        %=============================================================
        function LF = full(L)
            LF = L * eye(size(L, 2));
        end

        %=============================================================
        %> @brief Returns the LazyLinearOp for the vertical concatenation [L, varargin{:}].
        %>
        %> @b See @b also: vertcat matlab built-in.
        %=============================================================
        function LV = vertcat(L, op)
            import matfaust.lazylinop.LazyLinearOp
            if size(L, 2) ~= size(op, 2)
                error('L and op numbers of columns must be the same.')
            end
            if ~ LazyLinearOp.isLazyLinearOp(op)
                op = LazyLinearOp.create_from_op(op);
            end

            function mul = vstack_mul(L, op, o)
                if ~ LazyLinearOp.isLazyLinearOp(o) && (ismatrix(o) || isscalar(o)) && isnumeric(o)
                    mul = [L * o ; op * o];
                else
                    mul = [L ; op] * o;
                end
            end

            function VI = vstack_id(L, op, S)
                rid = S.subs{1};
                if strcmp(rid, ':')
                    rid = 1:size(L, 1) + size(op, 1);
                end
                if all(rid <= size(L, 1))
                    % the indexing of [L ; op] is limited to L
                    VI = subsref(L, S);
                elseif all(rid > size(L, 1))
                    % the indexing of [L ; op] is limited to op
                    rid = rid - size(L, 1);
                    S.subs{1} = rid;
                    VI = subsref(op, S);
                else
                    % the indexing overlaps L and op rows
                    % verify that rid is ordered increasingly
                    if numel(rid) > 1
                        for i=1:numel(rid)-1
                            j = i + 1;
                            if rid(j) < rid(i)
                                error(['indices must be in increasing order on a',
                                ' LazyLinearOp obtained by concatenation'])
                            end
                        end
                    end
                    SL.type = '()';
                    SL.subs = {rid(rid <= size(L, 1)), S.subs{2}};
                    idL = subsref(L, SL);
                    Sop.type = '()';
                    Sop.subs = {rid(rid > size(L, 1)) - size(L, 1), S.subs{2}};
                    idop = subsref(op, Sop);
                    VI = [idL ; idop];
                end
            end

            lambdas = {@(o) vstack_mul(L, op, o), ... %MUL
                @() horzcat(L.', op.'), ... % T
                @() horzcat(L', op'), ... % H
                @(S) vstack_id(L, op, S) % I
            };

            LV = LazyLinearOp(lambdas, [size(L, 1) + size(op, 1), size(L, 2)]);
        end


        %=============================================================
        %> @brief Returns the LazyLinearOp for the horizontal concatenation [L, varargin{:}].
        %>
        %> @b See @b also: horzcat matlab built-in.
        %=============================================================
        function LH = horzcat(L, op)
            import matfaust.lazylinop.LazyLinearOp
            if size(L, 1) ~= size(op, 1)
                error('L and op numbers of rows must be the same.')
            end
            if ~ LazyLinearOp.isLazyLinearOp(op)
                op = LazyLinearOp.create_from_op(op);
            end

            function mul = hstack_mul(L, op, o)
                % TODO: check dimensions
                mul = L * o(1:size(L, 2), :) + op * o(size(L,2)+1:size(o, 2), :)
            end

            function HI = hstack_id(L, op, S)
                cid = S.subs{2};
                if strcmp(cid, ':')
                    cid = 1:size(L, 2) + size(op, 2);
                end
                if all(cid <= size(L, 2))
                    % the indexing of [L op] is limited to L
                    HI = subsref(L, S);
                elseif all(cid > size(L, 2))
                    % the indexing of [L ; op] is limited to op
                    cid = cid - size(L, 2);
                    S.subs{2} = cid;
                    HI = subsref(op, S);
                else
                    % the indexing overlaps L and op rows
                    % verify that cid is ordered increasingly
                    if numel(cid) > 1
                        for i=1:numel(cid)-1
                            j = i + 1;
                            if cid(j) < cid(i)
                                error(['indices must be in increasing order on a',
                                ' LazyLinearOp obtained by concatenation'])
                            end
                        end
                    end
                    SL.type = '()';
                    SL.subs = {S.subs{1}, cid(cid <= size(L, 2))};
                    idL = subsref(L, SL);
                    Sop.type = '()';
                    Sop.subs = {S.subs{1}, cid(cid > size(L, 2)) - size(L, 1)};
                    idop = subsref(op, Sop);
                    HI = [idL idop];
                end
            end

            lambdas = {@(o) hstack_mul(L, op, o), ... %MUL
                @() horzcat(L.', op.'), ... % T
                @() horzcat(L', op'), ... % H
                @(S) hstack_id(L, op, S) % I
            };

            LH = LazyLinearOp(lambdas, [size(L, 1), size(L, 2) + size(op, 2)]);
        end

        %=============================================================
        %> @brief Returns the LazyLinearOp for the concatenation [L, varargin{:}] or [L; varargin{:}].
        %>
        %> @b See @b also: cat matlab built-in.
        %=============================================================
		function LC = cat(varargin)
			import matfaust.lazylinop.LazyLinearOp
			dim = varargin{1};
			L = varargin{2};
			check_meth(L, 'cat');
			for i=3:nargin
				O = varargin{i};
				if dim == 1
                    L = [L ; O];
				elseif dim == 2
                    L = [L, O];
				end
			end
            LC = L;
		end

        %=============================================================
        %> @brief Returns the LazyLinearOp for real(L).
        %>
        %> @b See @b also: real matlab built-in.
        %=============================================================
        function LR = real(L)
			import matfaust.lazylinop.LazyLinearOp
            function RM = real_mul(L, o)
                if LazyLinearOp.isLazyLinearOp(o) && ismatrix(o) && isnumeric(o) && any(size(o) ~= [1, 1])
                    RM = real(L * real(o)) + real(L * imag(o) * 1j);
                else
                    RM = real(L * o);
                end
            end
            lambdas = {@(o) real_mul(L, o), ... %MUL
                @() real(L.'), ... % T
                @() real(L.'), ... % H
                @(S) real(subsref(L, S)) % I
            };
            LR = LazyLinearOp(lambdas, size(L));
        end

        %=============================================================
        %> @brief Returns the LazyLinearOp for imag(L).
        %>
        %> @b See @b also: imag matlab built-in.
        %=============================================================
        function LI = imag(L)
			import matfaust.lazylinop.LazyLinearOp
            function IM = imag_mul(L, o)
                if LazyLinearOp.isLazyLinearOp(o) && ismatrix(o) && isnumeric(o) && any(size(o) ~= [1, 1])
                    IM = imag(L * real(o)) + imag(L * imag(o) * 1j);
                else
                    IM = imag(L * o)
                end
            end
            lambdas = {@(o) imag_mul(L, o), ... %MUL
                @() imag(L.'), ... % T
                @() imag(L.'), ... % H
                @(S) imag(subsref(L, S)) % I
            };
            LI = LazyLinearOp(lambdas, size(L));
        end

    end
    methods(Static = true)
        %=============================================================
        %> Alias of matfaust.lazylinop.isLazyLinearOp.
        %=============================================================
        function b = isLazyLinearOp(obj)
            b = isa(obj, 'matfaust.lazylinop.LazyLinearOp');
        end

        function lop = create_from_op(obj)
            import matfaust.lazylinop.LazyLinearOp
            lambdas = cell(1, 4);
            lambdasT = cell(1, 4);
            lambdasH = cell(1, 4);
            lambdasC = cell(1, 4);

            lambdas{LazyLinearOp.MUL} = @(op) obj * op;
            lambdasT{LazyLinearOp.MUL} = @(op) obj.' * op;
            lambdasH{LazyLinearOp.MUL} = @(op) obj' * op;
            lambdasC{LazyLinearOp.MUL} = @(op) conj(obj) * op;

            if ~ isreal(obj)
                dtype = 'complex';
            elseif strcmp('double', class(obj))
                dtype = 'double';
            else
                dtype = 'single';
            end

            lop = LazyLinearOp(lambdas, size(obj), 'dtype', dtype, ...
                'root_obj', obj);
            lopT = LazyLinearOp(lambdasT, [size(obj, 2), size(obj, 1)], 'dtype', dtype, ...
                'root_obj', obj);
            lopH = LazyLinearOp(lambdasH, [size(obj, 2), size(obj, 1)], 'dtype', dtype, ...
                'root_obj', obj);
            lopC = LazyLinearOp(lambdasC, size(obj), 'dtype', dtype, ...
                'root_obj', obj);

            lop.lambdas{LazyLinearOp.T} = @() lopT;
            lop.lambdas{LazyLinearOp.H} = @() lopH;
            lop.lambdas{LazyLinearOp.I} = @(S) index_lambda(lop, S);

            lopT.lambdas{LazyLinearOp.T} = @() lop;
            lopT.lambdas{LazyLinearOp.H} = @() lopC;
            lopT.lambdas{LazyLinearOp.I} = @(S) index_lambda(lopT, S);

            lopH.lambdas{LazyLinearOp.T} = @() lopC;
            lopH.lambdas{LazyLinearOp.H} = @() lop;
            lopH.lambdas{LazyLinearOp.I} = @(S) index_lambda(lopH, S);

            lopC.lambdas{LazyLinearOp.T} = @() lopH;
            lopC.lambdas{LazyLinearOp.H} = @() lopT;
            lopC.lambdas{LazyLinearOp.I} = @(S) index_lambda(lopC, S);
        end

        function lop = create_from_funcs(matmat, rmatmat, shape, varargin)
            import matfaust.lazylinop.LazyLinearOp

            p = inputParser;
            validDtype = @(dtype) any(strcmp(dtype, {'complex', 'double', 'single', 'undefined'}));
            addOptional(p, 'dtype', 'undefined', validDtype)

            parse(p, varargin{:})
            dtype = p.Results.dtype;

            MX = @(X) matmat(X);
            %MTX = @(X) rmatmat(X.').';
            MHX = @(X) rmatmat(X);

            lambdas = cell(1, 4);
            lambdasT = cell(1, 4);
            lambdasH = cell(1, 4);
            lambdasC = cell(1, 4);

            lambdas{LazyLinearOp.MUL} = MX;
            lambdasT{LazyLinearOp.MUL} = @(op)  conj(rmatmat(real(op))) - ...
                conj(rmatmat(1j * imag(op)));
            lambdasH{LazyLinearOp.MUL} = MHX;
            lambdasC{LazyLinearOp.MUL} = @(op) conj(matmat(real(op))) - ...
                matmat(1j * imag(op));

            lop = LazyLinearOp(lambdas, shape, 'dtype', dtype);
            lopT = LazyLinearOp(lambdasT, [shape(2), shape(1)], 'dtype', dtype);
            lopH = LazyLinearOp(lambdasH, [shape(2), shape(1)], 'dtype', dtype);
            lopC = LazyLinearOp(lambdasC, shape, 'dtype', dtype);

            lop.lambdas{LazyLinearOp.T} = @() lopT;
            lop.lambdas{LazyLinearOp.H} = @() lopH;
            lop.lambdas{LazyLinearOp.I} = @(S) index_lambda(lop, S);

            lopT.lambdas{LazyLinearOp.T} = @() lop;
            lopT.lambdas{LazyLinearOp.H} = @() lopC;
            lopT.lambdas{LazyLinearOp.I} = @(S) index_lambda(lopT, S);

            lopH.lambdas{LazyLinearOp.T} = @() lopC;
            lopH.lambdas{LazyLinearOp.H} = @() lop;
            lopH.lambdas{LazyLinearOp.I} = @(S) index_lambda(lopH, S);

            lopC.lambdas{LazyLinearOp.T} = @() lopH;
            lopC.lambdas{LazyLinearOp.H} = @() lopT;
            lopC.lambdas{LazyLinearOp.I} = @(S) index_lambda(lopC, S);
        end

    end
    methods(Access = protected)
        function check_meth(L, meth)
            if ~ strcmp(L.root_obj, 'none')  && ~ any(ismember(methods(L.root_obj), meth))
                error(meth+' is not supported by the root object of this LazyLinearOp')
            end
        end
    end
end
