%====================================
%> This function is a shorthand for matfaust.fact.hierarchical().
%===
%> <p>@b See @b also matfaust.fact.hierarchical
%====================================
function varargout = faust_fact(varargin)
	[F, lambda, p] = matfaust.fact.hierarchical(varargin{:});
	varargout = {F, lambda, p};
end
