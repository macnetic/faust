% ======================================================================
%>    This class implements a lazy linear operator.
%>
%>    The evaluation of any defined operation is delayed.
%>
%>    For creation and evaluation look at matfaust.lazylinop.asLazyLinearOp and
%>    LazyLinearOp.eval.
%>
%> @warning This code is in a beta status.
% ======================================================================
classdef LazyLinearOp
	properties (SetAccess = private, Hidden = true)
		lambda_stack
		root_obj
		shape
		class % reserved but not used
	end
	methods
		function L = LazyLinearOp(init_lambda, shape, root_obj)
			L.lambda_stack = init_lambda;
			L.shape = shape;
			L.root_obj = root_obj;
			L.class = 0;
		end

		function check_meth(obj, meth)
			b = any(ismember(methods(obj), meth));
			if ~ b
				error(meth+' is not supported by the rootobject of this LazyLinearOp')
			end
		end

		function s = size(L, varargin)
			if length(varargin) == 0
				s = L.shape;
			else
				s = L.shape(varargin{1});
			end
		end

		function LT = transpose(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'transpose');
			Ls = size(L);
			LT = LazyLinearOp(@() transpose(L.lambda_stack()), [Ls(2), Ls(1)], L.root_obj);
		end

		function LCT = ctranspose(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'ctranspose');
			Ls = size(L);
			LT = LazyLinearOp(@() ctranspose(L.lambda_stack()), [Ls(2), Ls(1)], L.root_obj);
		end

		function LC = conj(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'conj');
			Ls = size(L);
			LC = LazyLinearOp(@() conj(L.lambda_stack()), Ls, L.root_obj);
		end

		function Lp = plus(L, op)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'plus')
			if ~ all(size(op) == [1, 1]) && ~ all(size(L) == size(op))
				error('Dimensions must agree')
			end
			Lp = LazyLinearOp(@() L.lambda_stack() + LazyLinearOp.eval_if_lazy(op), size(L), L.root_obj);
		end

		function LUP = uplus(L) 
			LUP = L;
		end

		function Lm = minus(L, op)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'minus')
			if ~ all(size(op) == [1, 1]) && ~ all(size(L) == size(op))
				error('Dimensions must agree')
			end
			Lm = LazyLinearOp(@() L.lambda_stack() - LazyLinearOp.eval_if_lazy(op), size(L), L.root_obj);
		end

		function O = eval(L)
			O = L.lambda_stack();
		end

		function Lm = mtimes(L, op)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'mtimes')
			op_is_scalar = all(size(op) == [1, 1]);
			if ~ op_is_scalar && ~ all(size(L, 2) == size(op, 1))
				error('Dimensions must agree')
			end
			if op_is_scalar
				new_size = L.shape;
			else
				new_size = [size(L, 1), size(op, 2)];
			end
			Lm = LazyLinearOp(@() L.lambda_stack() * LazyLinearOp.eval_if_lazy(op), new_size, L.root_obj);
		end

		function Lm = mrdivide(L, s)
			if(~ isscalar(s))
				error('Unsupported operand type(s) for /: a LazyLinearOp can only be divided by a scalar.')
			end
			Lm = mtimes(L, 1/s);
		end

		function D = full(L)
			D = eval(L);
		end

		function SUB = subsref(L, S)

			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'subsref')

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

			SUB = LazyLinearOp(@() subsref(L.lambda_stack(), S) , new_shape, L.root_obj);
		end

		function LH = horzcat(L, O)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'horzcat')
			LH = cat(2, L, O)
		end

		function LV = vertcat(L, O)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'vertcat')
			LV = cat(1, L, O)
		end

		function LC = cat(dim, L, O)
			import matfaust.lazylinop.LazyLinearOp
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
			LC = LazyLinearOp(@() cat(dim, L.lambda_stack(), LazyLinearOp.eval_if_lazy(O)), new_size, L.root_obj);
		end

		function LR = real(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'real')
			LR = LazyLinearOp(@() real(L.lambda_stack()) , L.shape, L.root_obj);
		end

		function LI = imag(L)
			import matfaust.lazylinop.LazyLinearOp
			check_meth(L, 'real')
			LI = LazyLinearOp(@() imag(L.lambda_stack()) , L.shape, L.root_obj);
		end

	end
	methods(Static, Access = public)
		function L = create(obj)
			import matfaust.lazylinop.LazyLinearOp
			L = LazyLinearOp(@() obj, size(obj), obj);
		end

		function B = isLazyLinearOp(obj)
			B = isa(obj, 'matfaust.lazylinop.LazyLinearOp');
		end
	end
	methods(Static, Access = private)
		function O = eval_if_lazy(obj)
			import matfaust.lazylinop.LazyLinearOp
			if LazyLinearOp.isLazyLinearOp(obj)
				O = eval(obj);
			else
				O = obj;
			end
		end
	end
end
