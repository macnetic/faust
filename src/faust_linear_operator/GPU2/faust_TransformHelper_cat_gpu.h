/**
 * These functions are for Faust::TransformHelper concatenations.
 */
namespace Faust
{
	template<typename FPP>
		TransformHelper<FPP,GPU2>* vertcat(const std::vector<TransformHelper<FPP,GPU2>*> & THs);

	template<typename FPP>
		TransformHelper<FPP,GPU2>* horzcat(const std::vector<TransformHelper<FPP,GPU2>*> & THs);

}
