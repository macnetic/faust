namespace Faust
{
	template<typename FPP>
	MatGeneric<FPP,GPU2>::~MatGeneric()
	{
	}

	template<typename FPP>
	MatGeneric<FPP,GPU2>::MatGeneric() : is_zeros(false), is_identity(false)
	{
	}
}
