#ifndef FAUST_Transform_ALGEBRA_H
#define FAUST_Transform_ALGEBRA_H

#include "faust_constant.h"
#include "faust_MatSparse.h"

template<typename FPP,Device DEVICE> class MatDense;
template<typename FPP,Device DEVICE> class Vect;
template<typename FPP,Device DEVICE> class Transform;


template<typename FPP>
FPP power_iteration(const Faust::Transform<FPP,Cpu> & A, const int nbr_iter_max,FPP threshold, int & flag);


template<typename FPP>
Faust::Vect<FPP,Cpu> operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu> & v);

template<typename FPP>
Faust::MatDense<FPP,Cpu> operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu> & M);

#include "faust_transform_algebra.hpp"

#endif
