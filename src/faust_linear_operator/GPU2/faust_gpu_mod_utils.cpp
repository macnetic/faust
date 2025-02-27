#define __GM_LOADER__
#include "gm_interf.h"
#include "faust_gpu_mod_utils.h"


namespace Faust {

	GPUModHandler* GPUModHandler::singleton = nullptr;

	GPUModHandler::GPUModHandler() : gm_handle(nullptr), marr_funcs_double(nullptr), marr_funcs_float(nullptr), marr_funcs_cuComplex(nullptr), marr_funcs_cuDoubleComplex(nullptr), dsm_funcs_double(nullptr), dsm_funcs_float(nullptr), dsm_funcs_cuComplex(nullptr), dsm_funcs_cuDoubleComplex(nullptr), spm_funcs_double(nullptr), spm_funcs_float(nullptr), spm_funcs_cuComplex(nullptr), spm_funcs_cuDoubleComplex(nullptr), gp_funcs_(nullptr), bsr_funcs_double(nullptr), bsr_funcs_float(nullptr), bsr_funcs_cuComplex(nullptr), bsr_funcs_cuDoubleComplex(nullptr)
	{
		// class is only instatiatable with get_singleton (private ctor)
	}

	GPUModHandler::~GPUModHandler()
	{
		if(gm_handle != nullptr)
			gm_close_lib(gm_handle);
		delete_mat_funcs(float); // it looks strange because load_all_mat_funcs is a macro (using a string tokenizer)
		delete_mat_funcs(double); // look def in faust_gpu_mod_utils.h
		delete_mat_funcs(cuComplex);
		delete_mat_funcs(cuDoubleComplex);
		this->singleton = nullptr;
	}

	GPUModHandler* GPUModHandler::get_singleton(const bool silent/*=false*/)
	{
		if(GPUModHandler::singleton == nullptr)
		{
			GPUModHandler::singleton = new GPUModHandler();
			// don't warn user at instantiation time about initializing the lib handler
		}
		else if(GPUModHandler::singleton->gm_handle == nullptr && ! silent)
		{
			cerr << "WARNING: you must call enable_gpu_mod() before using GPUModHandler singleton." << endl;
		}
		return GPUModHandler::singleton;
	}

	void GPUModHandler::load_gm_functions()
	{
		if(gp_funcs_ == nullptr)
		{
			gp_funcs_ = new gm_GenPurposeFunc();
			load_gp_funcs(gm_handle, gp_funcs_);
		}

		if(marr_funcs_double == nullptr)
		{
			load_all_mat_funcs(double); // it looks strange because load_all_mat_funcs is a macro (using a string tokenizer)
			load_all_mat_funcs(float);  // look the definition in faust_gpu_mod_utils.h
			load_all_mat_funcs(cuComplex);
			load_all_mat_funcs(cuDoubleComplex);
		}
	}

	void* GPUModHandler::init_gpu_mod(const string& libpath, bool silent, void* gm_handle)
	{
		if(this->gm_handle == nullptr)
			if(gm_handle == nullptr)
				this->gm_handle = gm_load_lib(libpath.c_str(), silent);
			else
				this->gm_handle = gm_handle;
		else
			if(! silent)
				std::cerr << "Warning: gm_lib is already loaded (can't reload it)." << endl;
		if(this->gm_handle != nullptr)
			load_gm_functions();
		return this->gm_handle;
	}

	void* GPUModHandler::enable_gpu_mod(const string& libpath, bool silent) //TODO: add backend argument: cuda (only available for now), opencl (?)
	{
		gm_handle = init_gpu_mod(libpath, silent, nullptr);
		return gm_handle;
	}

	void* GPUModHandler::enable_gpu_mod(const char* libpath, const bool silent)
	{
		return enable_gpu_mod(string(libpath), silent);
	}

	void GPUModHandler::check_gpu_mod_loaded() const
	{
		if(gm_handle == nullptr)
			throw std::runtime_error("Faust::enable_gpu_mod() must be called before any use of gpu_mod.");
	}

	bool GPUModHandler::is_gpu_mod_loaded() const
	{
		return gm_handle != nullptr && gp_funcs_ != nullptr && gp_funcs_->free_mat != nullptr;
	}

	gm_SparseMatFunc_double* GPUModHandler::spm_funcs(const double &d) const
	{
		return spm_funcs_double;
	}

	gm_SparseMatFunc_float* GPUModHandler::spm_funcs(const float &d) const
	{
		return spm_funcs_float;
	}

	gm_SparseMatFunc_cuComplex* GPUModHandler::spm_funcs(const complex<float> &c) const
	{
		return spm_funcs_cuComplex;
	}

	gm_SparseMatFunc_cuDoubleComplex* GPUModHandler::spm_funcs(const complex<double> &c) const
	{
		return spm_funcs_cuDoubleComplex;
	}

	gm_DenseMatFunc_double* GPUModHandler::dsm_funcs(const double &d) const
	{
		return dsm_funcs_double;
	}

	gm_DenseMatFunc_float* GPUModHandler::dsm_funcs(const float &d) const
	{
		return dsm_funcs_float;
	}

	gm_DenseMatFunc_cuComplex* GPUModHandler::dsm_funcs(const complex<float> &c) const
	{
		return dsm_funcs_cuComplex;
	}

	gm_DenseMatFunc_cuDoubleComplex* GPUModHandler::dsm_funcs(const complex<double> &c) const
	{
		return dsm_funcs_cuDoubleComplex;
	}

	gm_MatArrayFunc_double* GPUModHandler::marr_funcs(const double &d) const
	{
		return marr_funcs_double;
	}

	gm_MatArrayFunc_float* GPUModHandler::marr_funcs(const float &d) const
	{
		return marr_funcs_float;
	}

	gm_MatArrayFunc_cuComplex* GPUModHandler::marr_funcs(const complex<float> &c) const
	{
		return marr_funcs_cuComplex;
	}

	gm_MatArrayFunc_cuDoubleComplex* GPUModHandler::marr_funcs(const complex<double> &c) const
	{
		return marr_funcs_cuDoubleComplex;
	}

	gm_BSRMatFunc_double* GPUModHandler::bsr_funcs(const double &d) const
	{
		return bsr_funcs_double;
	}

	gm_BSRMatFunc_float* GPUModHandler::bsr_funcs(const float &d) const
	{
		return bsr_funcs_float;
	}

	gm_BSRMatFunc_cuComplex* GPUModHandler::bsr_funcs(const complex<float> &c) const
	{
		return bsr_funcs_cuComplex;
	}

	gm_BSRMatFunc_cuDoubleComplex* GPUModHandler::bsr_funcs(const complex<double> &c) const
	{
		return bsr_funcs_cuDoubleComplex;
	}

	gm_GenPurposeFunc* GPUModHandler::gp_funcs() const
	{
		return gp_funcs_;
	}

	void* enable_gpu_mod(const char* libpath, const bool silent)
	{
		return GPUModHandler::get_singleton()->enable_gpu_mod(libpath, silent);
	}

	bool is_gpu_mod_enabled()
	{
		return GPUModHandler::get_singleton(true)->is_gpu_mod_loaded();
	}

	void char2gm_Op(const char& c, gm_Op & op)
	{
		if(c == 'N')
			op = OP_NOTRANSP;
		else if(c == 'T')
			op = OP_TRANSP;
		else if(c == 'H')
			op = OP_CONJTRANSP;
		else
			throw std::runtime_error("invalid character to convert to gm_Op");
	}


}
