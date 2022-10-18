% ======================================================================
%> @brief This class implements a lazy linear operator.
%>
%> The evaluation of any defined operation is delayed until proceeding to a multiplication
%> by a dense matrix/vector, a call of LazyLinearOp.toarray or an explicit
%> evaluation through LazyLinearOp.eval.
%>
%> To instantiate a LazyLinearOp look at matfaust.lazylinop.asLazyLinearOp.
%>
%> @warning This code is in a beta status.
% ======================================================================
classdef LazyLinearOp
	properties (SetAccess = protected, Hidden = true)
		lambda_stack;
		root_obj;
		shape;
		class; % reserved but not used
	end
	methods
		%======================================================================
		%> @brief Constructor. Not meant to be used directly.
		%>
		%> @param init_lambda: starting operation.
		%> @param shape: the initial shape of the operator.
		%> @param root_obj: the initial object the operator is based on.
		%>
		%>
		%> <p>@b See @b also LazyLinearOp.create, matfaust.lazylinop.asLazyLinearOp.
		%======================================================================
		function L = LazyLinearOp(init_lambda, shape, root_obj)
			L.lambda_stack = init_lambda;
			L.shape = shape;
			L.root_obj = root_obj;
			L.class = 0;
		end

		%======================================================================
		%> @brief size of L.
		%>
		%> @param dim: the dimension (1 for rows, 2 for columns) to get the size of.
		%>
		%> @retval s: by default, an array of two sizes (the number of rows and the number of columns) or the size of only one dimension if dim argument is used.
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
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'transpose');
			Ls = size(L);
			LT = LazyLinearOp(@() transpose(eval(L)), [Ls(2), Ls(1)], L.root_obj);
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp ctranspose.
		%=============================================================
		function LCT = ctranspose(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'ctranspose');
			Ls = size(L);
			LCT = LazyLinearOp(@() ctranspose(eval(L)), [Ls(2), Ls(1)], L.root_obj);
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp conjugate.
		%=============================================================
		function LC = conj(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'conj');
			Ls = size(L);
			LC = LazyLinearOp(@() conj(eval(L)), Ls, L.root_obj);
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp L + op.
		%>
		%> @param op: an object compatible with self for this binary operation.
		%=============================================================
		function Lp = plus(L, op)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'plus');
			if ~ all(size(op) == [1, 1]) && ~ all(size(L) == size(op))
				error('Dimensions must agree')
			end
			Lp = LazyLinearOp(@() eval(L) + LazyLinearOp.eval_if_lazy(op), size(L), L.root_obj);
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp + L.
		%=============================================================
		function LUP = uplus(L)
			LUP = L;
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp L - op.
		%>
		%> @param op: an object compatible with self for this binary operation.
		%=============================================================
		function Lm = minus(L, op)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'minus');
			if ~ all(size(op) == [1, 1]) && ~ all(size(L) == size(op))
				error('Dimensions must agree')
			end
			Lm = LazyLinearOp(@() eval(L) - LazyLinearOp.eval_if_lazy(op), size(L), L.root_obj);
		end

		%=============================================================
		%> @brief Evaluate the LazyLinearOp. All stacked virtual operations are evaluated.
		%>
		%> @retval O: a linear operator object (whose type depends of the LazyLinearOp initialization through matfaust.lazylinop.asLazyLinearOp and the operations made upon this object).
		%>
		%> @b Example
		%> @code
		%> >> import matfaust.lazylinop.asLazyLinearOp
		%> >> import matfaust.rand
		%> >> F = rand(10, 12);
		%> >> lF = asLazyLinearOp(F)
		%>
		%> lF =
		%>
		%>   10x12 LazyLinearOp array with no properties.
		%>
		%> >> twolF = 2 * lF
		%>
		%> twolF =
		%>
		%>   10x12 LazyLinearOp array with no properties.
		%>
		%> >> eval(twolF)
		%>
		%> ans =
		%>
		%> Faust size 10x12, density 2.025, nnz_sum 243, 5 factor(s):
		%> - FACTOR 0 (double) SPARSE, size 10x12, density 0.333333, nnz 40, addr: 0x7f62b6503670
		%> - FACTOR 1 (double) SPARSE, size 12x12, density 0.333333, nnz 48, addr: 0x7f62b623c2c0
		%> - FACTOR 2 (double) SPARSE, size 12x11, density 0.454545, nnz 60, addr: 0x7f62b623c730
		%> - FACTOR 3 (double) SPARSE, size 11x10, density 0.5, nnz 55, addr: 0x7f62b627c3e0
		%> - FACTOR 4 (double) SPARSE, size 10x12, density 0.333333, nnz 40, addr: 0x7f62b623cb30
		%> >> norm(full(eval(twolF)) - full(2*F))
		%>
		%> ans =
		%>
		%>      0
		%>
		%> @endcode
		%=============================================================
		function O = eval(L)
			O = L.lambda_stack();
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp for the multiplication self * op
		%> or if op is a full matrix it returns the full matrix (self * op).
		%>
		%> @param op: an object compatible with self for this binary operation.
		%=============================================================
		function Lm = mtimes(L, op)
			import matfaust.lazylinop.LazyLinearOp
			if isscalar(L) && LazyLinearOp.isLazyLinearOp(op)
				Lm = mtimes(op, L);
				return;
			end
			check_meth(L, 'mtimes');
			op_is_scalar = all(size(op) == [1, 1]);
			if ~ op_is_scalar && ~ all(size(L, 2) == size(op, 1))
				error('Dimensions must agree')
			end
			if op_is_scalar
				new_size = L.shape;
			else
				new_size = [size(L, 1), size(op, 2)];
			end
			if ~ LazyLinearOp.isLazyLinearOp(op) && ismatrix(op) && isnumeric(op) && any(size(op) ~= [1, 1])
				% op is a dense matrix that is not limited to one element
				Lm = eval(L) * op;
			else
				Lm = LazyLinearOp(@() eval(L) * LazyLinearOp.eval_if_lazy(op), new_size, L.root_obj);
			end
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp for mrdivide(L, op)
		%>
		%> @param op: an object compatible with self for this binary operation.
		%>
		%> @b See @b also: mrdivide matlab built-in.
		%=============================================================
		function Lm = mrdivide(L, op)
			Lm = LazyLinearOp(@() mrdivide(eval(L) * LazyLinearOp.eval_if_lazy(op), new_size, L.root_obj));
		end

		%=============================================================
		%> @brief Returns the full matrix resulting from self evaluation.
		%=============================================================
		function D = full(L)
			D = eval(L);
			if is_meth(L, 'full');
				D = full(D);
			else
				error('full is not available on the result of eval')
			end
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp for indexing.
		%>
		%> @b See @b also: subsref matlab built-in.
		%=============================================================
		function SUB = subsref(L, S)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'subsref');

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

			if(numel(S.subs{1}) == 1 && numel(S.subs{2}) == 1 && ~ ischar(S.subs{1}) && ~ ischar(S.subs{2}))
				% accessing a single item
				new_shape = [1, 1];
			else

				end_ids = zeros(1,2);
				start_ids = zeros(1,2);
				slicing = [ true, true ];
				ind_lists = cell(1,2);
				ROW=1;
				COL=2;

				for i=ROW:COL
					ind_list=S.subs{i};
					if ischar(ind_list)
						start_ids(i) = 1;
						end_ids(i) = size(L,i);
						%slicing(i) = true
						ind_list = 1:size(L,i); % needed if other dim is not a slice
					else
						if(any(ind_list < 1))
							error(' Subscript indices must be integers >= 1.')
						elseif(any(ind_list > size(L,i)))
							error(' Index exceeds Faust dimensions.')
						elseif(size(ind_list,2) == 0)
							error(' Cannot create empty Faust')
						end
						% check if indices in range are contiguous and not negative step
						sl_sz = size(ind_list,2);
						if(sl_sz >= 2)
							for j=2:sl_sz
								d = ind_list(j)-ind_list(j-1);
								if(abs(d) > 1 || d < 0)
									slicing(i) = false;
									break
								end
							end
						end
						if(slicing(i))
							start_ids(i) = ind_list(1);
							end_ids(i) = ind_list(end);
						end
					end
					ind_lists{i} = ind_list;
				end

				if(slicing)
					% slicing
					new_shape = [end_ids(ROW) - start_ids(ROW) + 1, end_ids(COL) - start_ids(COL) + 1];
				else
					% indexing
					new_shape = [numel(ind_lists{1}), numel(ind_lists{2})];
				end
			end

			SUB = LazyLinearOp(@() subsref(eval(L), S) , new_shape, L.root_obj);
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp for the horizontal concatenation [L, varargin{:}].
		%>
		%> @b See @b also: horzcat matlab built-in.
		%=============================================================
		function LH = horzcat(varargin)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(varargin{1}, 'horzcat');
			LH = cat(2, varargin{1}, varargin{2:end});
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp for the vertical concatenation [L, varargin{:}].
		%>
		%> @b See @b also: vertcat matlab built-in.
		%=============================================================
		function LV = vertcat(varargin)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(varargin{1}, 'vertcat');
			LV = cat(1, varargin{1}, varargin{2:end});
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
					odim = 2;
					new_size = [size(L, 1) + size(O, 1), size(L, 2)];
				elseif dim == 2
					odim = 1;
					new_size = [size(L, 1), size(L, 2) + size(O, 2)];
				end
				if ~ all(size(L, odim) == size(O, odim))
					error('Dimensions must agree')
				end
				LC = LazyLinearOp(@() cat(dim, eval(L), LazyLinearOp.eval_if_lazy(O)), new_size, L.root_obj);
				L = LC;
			end
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp for real(L).
		%>
		%> @b See @b also: real matlab built-in.
		%=============================================================
		function LR = real(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'real');
			LR = LazyLinearOp(@() real(eval(L)) , L.shape, L.root_obj);
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp for imag(L).
		%>
		%> @b See @b also: imag matlab built-in.
		%=============================================================
		function LI = imag(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'real');
			LI = LazyLinearOp(@() imag(eval(L)) , L.shape, L.root_obj);
		end

	end
	methods(Static, Access = public)

		%=============================================================
		%> @brief Alias of matfaust.lazylinop.asLazyLinearOp.
		%>
		%> @param obj: cf. matfaust.lazylinop.asLazyLinearOp.
		%>
		%> @retval L:  cf. matfaust.lazylinop.asLazyLinearOp.
		%>
		%> @b Example
		%> @code
		%> >> import matfaust.lazylinop.LazyLinearOp
		%> >> M = rand(10, 12);
		%> >> lM = LazyLinearOp.create(M);
		%> >> twolM = lM + lM
		%>
		%> twolM =
		%>
		%>   10x12 LazyLinearOp array with no properties.
		%>
		%> >> F = matfaust.rand(10, 12);
		%> >> lF = LazyLinearOp.create(F);
		%> >> twolF = lF + lF
		%> twolF =
		%>
		%>    10x12 LazyLinearOp array with no properties.
		%> @endcode
		%>
		%> @b See @b also: matfaust.rand.
		%=============================================================
		function L = create(obj)
			import matfaust.lazylinop.LazyLinearOp
			L = LazyLinearOp(@() obj, size(obj), obj);
		end

		%=============================================================
		%> Alias of matfaust.lazylinop.isLazyLinearOp.
		%=============================================================
		function B = isLazyLinearOp(obj)
			B = isa(obj, 'matfaust.lazylinop.LazyLinearOp');
		end
	end
	methods(Static, Access = protected)
		function O = eval_if_lazy(obj)
			import matfaust.lazylinop.LazyLinearOp
			if LazyLinearOp.isLazyLinearOp(obj)
				O = eval(obj);
			else
				O = obj;
			end
		end
	end
	methods(Access = protected)
		function check_meth(obj, meth)
			if ~ is_meth(obj, meth)
				error(meth+' is not supported by the root object of this LazyLinearOp')
			end
		end

		function b = is_meth(obj, meth)
			b = any(ismember(methods(obj), meth));
		end
	end
end

%> @package matfaust.lazylinop @brief The matfaust module for lazy linear operators.
