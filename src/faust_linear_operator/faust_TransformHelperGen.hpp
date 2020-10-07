namespace Faust
{
	template<typename FPP, FDevice DEV>
		TransformHelperGen<FPP,DEV>::TransformHelperGen() : is_transposed(false), is_conjugate(false), is_sliced(false), is_fancy_indexed(false), transform(std::make_shared<Transform<FPP,DEV>>())
	{
	}

	template<typename FPP, FDevice DEV>
		const char TransformHelperGen<FPP,DEV>::isTransposed2char() const
		{
			return this->is_transposed?(this->is_conjugate?'H':'T'):'N';
		}

	template <typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::pack_factors(const faust_unsigned_int id, const PackDir dir)
	{
		if(dir == PACK_RIGHT)
			this->pack_factors(id, size()-1);
		else // dir == PACK_LEFT
			this->pack_factors(0, id);
	}

	template <typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::pack_factors()
	{
		//pack all factors in one
		this->pack_factors(0, this->size()-1);
	}
}
