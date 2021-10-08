% ======================================================================
%> @brief Subclass of Faust specialized for orthogonal polynomial basis.
%>
%> This class is used only for the native implementation of the poly functions.
%>
%> @note It is not advisable to use this class directly.
%>
%>
% ======================================================================
classdef FaustPoly < matfaust.Faust
	methods
% ======================================================================
%> @brief Constructor
% ======================================================================
		function F = FaustPoly(varargin)
			F = F@matfaust.Faust(varargin{:});
		end

		%================================================================
		%> Returns the next polynomial basis (one additional dimension compared to this one).
		%===
		%>
		%================================================================
		function M = next(self)
			if(self.isreal)
				if(strcmp(self.dtype, 'float'))
					core_obj = mexPolyRealFloat('nextPolyFaust', self.matrix.objectHandle);
				else
					core_obj = mexPolyReal('nextPolyFaust', self.matrix.objectHandle);
				end
			else
				core_obj = mexPolyCplx('nextPolyFaust', self.matrix.objectHandle);
			end
			M = matfaust.poly.FaustPoly(core_obj, self.isreal, 'cpu', self.dtype);
		end

	end

	methods(Static)
		%================================================================
		%> Non documented (implementation backend).
		%===
		%>
		%================================================================
		function M = poly_(self, coeffs, X)
			if(iscell(X))
				% X is {}: no X passed (see matfaust.poly.poly())
				if(self.isreal)
					if(strcmp(class(self), 'single'))
						core_obj = mexPolyRealFloat('polyFaust', coeffs, self.matrix.objectHandle);
					else
						core_obj = mexPolyReal('polyFaust', coeffs, self.matrix.objectHandle);
					end
				else
					core_obj = mexPolyCplx('polyFaust', coeffs, self.matrix.objectHandle);
				end
				M = matfaust.poly.FaustPoly(core_obj, self.isreal, 'cpu', self.dtype);
			elseif(ismatrix(X))
				if(issparse(X))
					error('X must be a dense matrix')
				end
				if(size(X, 1) ~= size(self,2))
					error('The faust and X dimensions must agree.')
				end
				if(self.isreal)
					if(strcmp(self.dtype, 'float'))
						M = mexPolyRealFloat('mulPolyFaust', coeffs, self.matrix.objectHandle, X);
					else
						M = mexPolyReal('mulPolyFaust', coeffs, self.matrix.objectHandle, X);
					end
				else
					M = mexPolyCplx('mulPolyFaust', coeffs, self.matrix.objectHandle, X);
				end
			end
		end

	end
end
