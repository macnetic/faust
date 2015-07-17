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
class faust_spmat;

class faust_vec
{
 public :
 faust_vec() : vec(), dim(0) {}
 faust_vec(int _dim) : vec(_dim),dim(_dim){}
 faust_vec(const faust_vec& v) : vec(v.vec), dim(v.dim){}
	
 faust_real* getData(){return vec.data();}
 const faust_real* getData() const {return vec.data();}
 void setOnes();
 void Display() const;
void print_file(const char* filename)const;
 int getDim() const {return dim;}
 int size() const {return dim;}
 void resize(const int new_dim);
faust_real norm(){return vec.norm();}
faust_real scalarMultiply(faust_real const scalar){vec = scalar * vec;}
void normalize(){scalarMultiply(1/norm());} 


// multiply (*this) =  A * (*this)
void  multiplyLeft(faust_mat const& A){gemv(A, *this, *this, 1.0, 0.0, 'N');}
void  multiplyLeft(faust_spmat const& A);
  
faust_real sum()const{return vec.sum();}
faust_real mean()const{return vec.mean();}


void operator=(faust_vec const& y);

void operator*=(const faust_real alpha);
void operator+=(const faust_real alpha);
void operator-=(const faust_real alpha);

void operator+=(const faust_vec& v);
void operator-=(const faust_vec& v);


faust_real mean_relative_error(const faust_vec& v);

faust_real& operator[](int i){return vec(i);}
const faust_real& operator[](int i)const{return vec(i);}

const faust_real& operator()(int i)const{return vec(i);}

friend void gemv(const faust_mat & A,const faust_vec & x,faust_vec & y,const faust_real & alpha, const faust_real & beta, char typeA);
friend faust_vec solve(const faust_mat & A, const faust_vec & v);
friend void solve(const faust_spmat & A,faust_vec & x, const faust_vec & y);

private:
  int dim;
  Eigen::Matrix<faust_real, Eigen::Dynamic,1> vec; 
};
#endif
