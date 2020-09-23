#ifndef __FAUST_CUDA_DEVICE__
#define __FAUST_CUDA_DEVICE__
#include <cstdint>
#include <functional>

namespace Faust
{
	/**
	 * \brief Changes the current device.
	 *
	 * \param dev_id the device identifier.
	 */
	void set_device(int32_t dev_id);

	/**
	 * \brief Gets the current device.
	 *
	 * \return the current device identifier.
	 */
	int32_t get_current_device();

	/**
	 * \brief Switches to another device and allows to switch back easily.
	 *
	 * \return a lambda to switch back to the previous device.
	 */
	std::function<void()> switch_device(int32_t dev_id);

	/**
	 * \brief Counts the number of devices available.
	 *
	 * \return the number of devices.
	 */
	int32_t count_devices();
}
#endif
