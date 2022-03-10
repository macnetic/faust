namespace Faust 
{
	template<typename FPP>
		MatBSR<FPP, GPU2>::MatBSR(const MatBSR<FPP, Cpu>& mat,
				const int32_t dev_id/*=-1*/,
				const void* stream/*=nullptr*/) : MatBSR<FPP,GPU2>(
					mat.getNbRow(),
					mat.getNbCol(),
					mat.getNbBlockRow(),
					mat.getNbBlockCol(),
					mat.getNBlocks(),
					mat.get_bdata(),
					mat.get_browptr(),
					mat.get_bcolinds(),
					dev_id,
					stream)

	{

	}

	template<typename FPP>
		MatBSR<FPP, GPU2>::MatBSR() : gpu_mat(nullptr)
		{
		}

	template<typename FPP>
		void MatBSR<FPP, GPU2>::setZeros()
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		size_t MatBSR<FPP, GPU2>::getNBytes() const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		MatType MatBSR<FPP, GPU2>::getType() const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		int32_t MatBSR<FPP, GPU2>::getNbRow() const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		int32_t MatBSR<FPP, GPU2>::getNbCol() const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		MatGeneric<FPP,GPU2>* MatBSR<FPP, GPU2>::clone(const int32_t dev_id/*=*-1*/, const void* stream/*=*nullptr*/) const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		MatGeneric<FPP,GPU2>* MatBSR<FPP, GPU2>::Clone(const bool isOptimize/*=*false*/) const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		void* MatBSR<FPP, GPU2>::get_gpu_mat_ptr() const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		faust_unsigned_int MatBSR<FPP, GPU2>::getNonZeros() const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		void MatBSR<FPP, GPU2>::transpose()
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		void MatBSR<FPP, GPU2>::conjugate()
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		void MatBSR<FPP, GPU2>::adjoint()
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		MatGeneric<FPP,GPU2>* MatBSR<FPP, GPU2>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		MatGeneric<FPP,GPU2>* MatBSR<FPP, GPU2>::get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		Faust::MatGeneric<FPP,GPU2>* MatBSR<FPP, GPU2>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		Faust::MatGeneric<FPP,GPU2>* MatBSR<FPP, GPU2>::get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}


	template<typename FPP>
		void MatBSR<FPP, GPU2>::Display() const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		Real<FPP> MatBSR<FPP, GPU2>::norm() const
		{
			//TODO: implement (maybe by moving into .cpp.in
		}

	template<typename FPP>
		void MatBSR<FPP, GPU2>::set_gpu_mat_ptr(void* gpu_mat)
		{
			this->gpu_mat = gpu_mat;
		}
};
