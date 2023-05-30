%================================
%> @brief Returns true if the gpu_mod plug-in has been loaded correctly, false otherwise.
%>
%> <p>@b See @b also matfaust.is_gpu_mod_working
%>
%===============================
function enabled = is_gpu_mod_enabled()
		enabled = mexFaustReal('is_gpu_mod_enabled');
end
