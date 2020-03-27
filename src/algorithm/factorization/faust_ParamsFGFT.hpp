template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::ParamsFGFT<FPP,DEVICE,FPP2>::Display() const
{
	Faust::Params<FPP,DEVICE,FPP2>::Display();
	cout << "init_Dthis->is_identity:" << init_D.is_id() << endl;
	cout <<  "init_D info:" << endl;
	init_D.Display();
	cout << "ParamsFGFT init_D norm: " << init_D.norm() << endl;
}

