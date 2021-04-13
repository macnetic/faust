#ifndef __FAUST_VECT_GPU2__
#define __FAUST_VECT_GPU2__
#define NOMINMAX // avoids VS min/max issue with std::min/max.
#include "faust_MatDense_gpu.h"
#include "faust_MatSparse_gpu.h"

namespace Faust
{
//	template <typename FPP, FDevice DEVICE>
//		class Vect<FPP,DEVICE>;

	template<typename FPP>
		class Vect<FPP,GPU2> : public MatDense<FPP,GPU2>
		{
			friend MatDense<FPP,GPU2>;
			friend MatSparse<FPP,GPU2>;
			friend Transform<FPP,GPU2>; // need to access to get_gpu_mat_ptr
			public:

			Vect();

			Vect(const faust_unsigned_int size,
					const FPP* cpu_data=nullptr,
					const bool no_alloc=false,
					const int32_t dev_id=-1,
					const void* stream=nullptr);

			faust_unsigned_int size() const;

			void resize(const faust_unsigned_int size);
			void operator=(const Vect<FPP,GPU2> &v);
			void operator=(const Vect<FPP,Cpu> &v);
			void operator*=(const Vect<FPP,GPU2> &v);
			void operator/=(const Vect<FPP,GPU2> &v);
			bool operator==(const Vect<FPP,GPU2> &v)const;
			bool operator!=(const Vect<FPP,GPU2> &v)const;
			FPP max_coeff();
			FPP min_coeff();
			FPP dot(const Vect<FPP,GPU2> &v);
			FPP sum() const;
			FPP mean() const;
			void Display() const;
			Vect<FPP, Cpu> tocpu(const void* stream=nullptr) const;
			void tocpu(FPP* cpu_buf, const void* stream=nullptr) const;
			void setValues(const FPP& val);
			FPP mean_relative_error(const Vect<FPP,GPU2>& ref_vec) const;
			// delete parent methods that don't apply to a vector
			void setEyes() = delete;
			void multiplyLeft(MatSparse<FPP,GPU2> const& S, const char transS='N');
			FPP operator[](faust_unsigned_int i);
//			FPP& operator[](faust_unsigned_int i); // impossible to return a reference to a GPU buffer value (therefore it's impossible to modify a vector coeff as it's done with Vect<FPP,Cpu> (e.g. v[i] = 3)
			void set_coeff(faust_unsigned_int i, const FPP& val);
		};
}

#endif
