%====================================
%> @brief This function returns true if the gpu_mod is working properly false otherwise.
%>
%> is_gpu_mod_working comes as a complement of matfaust.is_gpu_mod_enabled. The latter ensures that gpu_mod shared library/plugin is properly loaded in memory but doesn't ensure that the GPU is available (for example, the NVIDIA driver might not be installed). The former ensures both that the gpu_mod is loaded and the GPU (device 0) is properly available for computing.
%>
%====================================
function success = is_gpu_mod_working()
	success = false;
	if matfaust.is_gpu_mod_enabled()
		try
			gpuF = matfaust.rand(1, 1, 'dev', 'gpu');
		catch
			sucess = false;
			return;
		end
		success = exist('gpuF');
	end
	% TODO: test another device (dev argument)
end

