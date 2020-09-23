#include "faust_cuda_device.h"
#include "faust_gpu_mod_utils.h"

namespace Faust
{

	void set_device(int32_t dev_id)
	{
		auto gp_funcs = GPUModHandler::get_singleton()->gp_funcs();
		gp_funcs->set_dev(dev_id);
	}

	int32_t get_current_device()
	{
		auto gp_funcs = GPUModHandler::get_singleton()->gp_funcs();
		return gp_funcs->cur_dev();
	}

	std::function<void()> switch_device(int32_t dev_id)
	{
		auto cdev_id = get_current_device();
		set_device(dev_id);
		return [=](){set_device(dev_id);};
	}

	int32_t count_devices()
	{
		auto gp_funcs = GPUModHandler::get_singleton()->gp_funcs();
		return gp_funcs->dev_count();
	}

}
