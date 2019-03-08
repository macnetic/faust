template<typename FPP,Device DEVICE,typename FPP2>
void Faust::ParamsFGFT<FPP,DEVICE,FPP2>::Display() const
{
	Faust::Params<FPP,DEVICE,FPP2>::Display();
	cout << "init_D isIdentity:" << init_D.estIdentite() << endl;
	cout <<  "init_D info:" << endl;
	init_D.Display();
	cout << "ParamsFGFT init_D norm: " << init_D.norm() << endl;
}

