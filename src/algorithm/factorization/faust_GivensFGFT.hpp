template<typename FPP, Device DEVICE>
void Faust::GivensFGFT<FPP,DEVICE>::next_step()
{

}

template<typename FPP, Device DEVICE>
void Faust::GivensFGFT<FPP,DEVICE>::compute_facts()
{

}

template<typename FPP, Device DEVICE>
Faust::GivensFGFT<FPP,DEVICE>::GivensFGFT(Faust::MatDense<FPP,DEVICE>& Lap, unsigned int j) : L(Lap), facts(j), D(Lap.getNbRow(), Lap.getNbCol()), C(Lap.getNbRow(), Lap.getNbRow()), err(j), coord_choices(j)
{

}
