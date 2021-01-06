/**
 * These functions are for Faust::TransformHelper concatenations.
 */
namespace Faust
{
	template<typename FPP>
		TransformHelper<FPP,Cpu>* vertcat(const std::vector<TransformHelper<FPP,Cpu>*> & THs);

	template<typename FPP>
		TransformHelper<FPP,Cpu>* horzcat(const std::vector<TransformHelper<FPP,Cpu>*> & THs);
}
