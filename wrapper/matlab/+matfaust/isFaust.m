
%================================================================
%> Returns true if obj is a Faust object, false otherwise.
%===
%>
%> @b Example
%> @code
%> import matfaust.*
%> isFaust(1) % returns 0
%> isFaust(FaustFactory.rand(5,10)) % returns 1
%> @endcode
%>
%> <p> @b See @b also Faust.Faust
%================================================================
function bool = isFaust(obj)
	import matfaust.Faust.isFaust
	bool = isFaust(obj);
end
