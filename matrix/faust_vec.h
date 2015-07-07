#ifndef FAUST_VEC_H
#define FAUST_VEC_H
 
#include <Eigen/Dense>
#include "faust_constant.h"
#include <vector>
#include <iterator>

#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif

class faust_mat;

class faust_vec
{
 public :
 faust_vec() : vec(), dim(0) {}
 faust_vec(int _dim) : vec(_dim),dim(_dim){} 
	
 faust_real* getData(){return vec.data();}
 const faust_real* getData() const {return vec.data();}
 void setOnes();
 void Display() const;
 int getDim() const {return dim;}
 void resize(const int new_dim);
faust_real norm(){return vec.norm();}
faust_real scalarMultiply(faust_real const scalar){vec = scalar * vec;}
void normalize(){scalarMultiply(1/norm());} 

void operator=(faust_vec const& y);
faust_real operator()(int i)const{return vec(i);}
friend void gemv(const faust_mat & A,const faust_vec & x,faust_vec & y,const faust_real & alpha, const faust_real & beta, char typeA);
private:
  int dim;
  Eigen::Matrix<faust_real, Eigen::Dynamic,1> vec; 
};
#endif