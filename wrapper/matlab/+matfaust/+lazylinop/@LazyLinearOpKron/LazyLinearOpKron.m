% ================================
%> @brief This class implements a specialization of LazyLinearOp dedicated to the
%> Kronecker product of two linear operators.
%>
%> @b See @b also: matfaust.lazylinop.kron.
% ================================
classdef LazyLinearOpKron < matfaust.lazylinop.LazyLinearOp
	properties (SetAccess = private, Hidden = true)
		A;
		B;
	end
	methods
		%================================
		%> Constructor for the A x B Kronecker product.
		%================================
		function LK = LazyLinearOpKron(A, B)

			shape = [size(A, 1) * size(B, 1), size(A, 2) * size(B, 2)];
			LK = LK@matfaust.lazylinop.LazyLinearOp(@() A, shape, A);
			LK.A = A;
			LK.B = B;
		end

		%================================
		%> Returns the LazyLinearOpKron conjugate.
		%================================
		function CLK = conj(LK)
			import matfaust.lazylinop.*
			CLK = LazyLinearOpKron(conj(asLazyLinearOp(LK.A)), conj(asLazyLinearOp(LK.B)));
		end

		%================================
		%> Returns the LazyLinearOpKron transpose.
		%================================
		function CLK = transpose(LK)
			import matfaust.lazylinop.*
			CLK = LazyLinearOpKron(transpose(asLazyLinearOp(LK.A)), transpose(asLazyLinearOp((LK.B))));
		end

		%================================
		%> Returns the LazyLinearOpKron adjoint/transconjugate.
		%================================
		function CTLK = ctranspose(LK)
			import matfaust.lazylinop.*
			CTLK = LazyLinearOpKron(ctranspose(asLazyLinearOp(LK.A)), ctranspose(asLazyLinearOp((LK.B))));
		end

		%=============================================================
		%> @brief Returns the LazyLinearOp for the multiplication self * op
		%> or if op is a full matrix it returns the full matrix (self * op).
		%>
		%> @note this specialization is particularly optimized for multiplying the
		%> operator by a vector.
		%>
		%> @param op: an object compatible with self for this binary operation.
		%=============================================================
		function LmK = mtimes(LK, op)
			import matfaust.lazylinop.LazyLinearOp
			if isscalar(LK) && LazyLinearOp.isLazyLinearOp(op)
				LmK = mtimes(op, LK);
				return;
			end
			check_meth(LK, 'mtimes');
			op_is_scalar = all(size(op) == [1, 1]);
			if ~ op_is_scalar && ~ all(size(LK, 2) == size(op, 1))
				error('Dimensions must agree')
			end
			if op_is_scalar
				new_size = LK.shape;
			else
				new_size = [size(LK, 1), size(op, 2)];
			end
			if ~ LazyLinearOp.isLazyLinearOp(op) && ismatrix(op) && isnumeric(op) && ~ issparse(op) && any(size(op) ~= [1, 1]) && any(ismember(methods(op), 'reshape')) && any(ismember(methods(op), 'mtimes')) %&& any(ismember(methods(op), 'subsref'))
				% op is a dense matrix that is not limited to one element
				LmK = zeros(new_size);
				A = LK.A;
				B = LK.B;
				for j=1:new_size(2)
					op_mat = reshape(op(:, j), [size(B, 2), size(A, 2)]);
					LmK(:, j) = reshape(LazyLinearOp.eval_if_lazy(B) * op_mat * transpose(LazyLinearOp.eval_if_lazy(A)), [new_size(1), 1]);
				end
			else
				LmK = LazyLinearOp(@() LK.eval() * LazyLinearOp.eval_if_lazy(op), new_size, LK.root_obj);
			end
		end

		%================================
		%> See LazyLinearOp.eval.
		%================================
		function E = eval(LK)
			A = LK.A;
			B = LK.B;
			if any(ismember(methods(A), 'full'))
				A = full(A);
			end
			if any(ismember(methods(B), 'full'))
				B = full(B);
			end
			E = kron(A, B);
		end

	end
end
