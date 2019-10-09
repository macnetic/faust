%====================================
%> This function is a shorthand for matfaust.FaustFactory.fact_hierarchical().
%===
%> <p>@b See @b also FaustFactory.fact_hierarchical
%====================================
function varargout = faust_fact(varargin)
	[F, lambda, p] = matfaust.fact.hierarchical(varargin{:});
	varargout = {F, lambda, p};
end
