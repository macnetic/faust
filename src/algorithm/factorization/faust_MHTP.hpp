namespace Faust
{
	template<typename FPP>
		MHTPParams<FPP>::MHTPParams() : used(false),
		constant_step_size(false),
		step_size(1e-3),
		palm4msa_period(1000),
		updating_lambda(true),
		sc(StoppingCriterion<FPP>(50)){};

	template<typename FPP>
	std::string MHTPParams<FPP>::to_string() const
	{
		auto sc_str = sc.to_string();
		std::string str = "MHTPParams (START):";
		str += "\r\n";
		str += "StoppingCriterion:";
		str += "\r\n ===";
		str += sc_str;
		str += " === \r\n";
		str += "constant_step_size: ";
		str += std::to_string(constant_step_size);
		str += "\r\n";
		str += "step_size: ";
		str += std::to_string(step_size);
		str += "\r\n";
		str += "palm4msa_period: ";
		str += std::to_string(palm4msa_period);
		str += "\r\n";
		str += "updating_lambda: ";
		str += std::to_string(updating_lambda);
		str += "\r\n";
		str += "MHTPParams END.";
		return str;
	}
}
