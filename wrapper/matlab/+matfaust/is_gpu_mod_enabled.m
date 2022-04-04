%================================
%> @brief Returns true if the gpu_mod plug-in has been loaded correctly, false otherwise.
%>
%===============================
function enabled = is_gpu_mod_enabled()
		enabled = mexFaustReal('is_gpu_mod_enabled');
end
