#ifndef FAUST_VEC_H
#define FAUST_VEC_H
 
#include <Eigen/Dense>
#include "faust_constant.h"
#include <vector>
#include <iterator>

#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif

template<typename T> class faust_mat;
template<typename T> class faust_spmat;
template<typename T> class faust_vec;



template<typename T> 
void gemv(const faust_mat<T> & A,const faust_vec<T> & x,faust_vec<T> & y,const T & alpha, const T & beta, char typeA);

template<typename T> 
faust_vec<T> solve(const faust_mat<T> & A, const faust_vec<T> & v);

template<typename T> 
void sp_solve(const faust_mat<T> & A,faust_vec<T> & x, const faust_vec<T> & y);



template<typename T>
class faust_vec
{
 public :
 faust_vec() : dim(0), vec() {}
 faust_vec(const int _dim) : dim(_dim), vec(_dim){}
 faust_vec(const faust_vec<T>& v) : dim(v.dim), vec(v.vec){}
 faust_vec(const faust_unsigned_int dim_, const T* data_);
	
 T* getData(){return vec.data();}
 const T* getData() const {return vec.data();}
 void setOnes();
 void Display() const;
void print_file(const char* filename)const;
 faust_unsigned_int size() const {return dim;}
 void resize(const int new_dim);
T norm(){return vec.norm();}
void scalarMultiply(T const scalar){vec *= scalar;}
void normalize(){scalarMultiply(1/norm());} 


// multiply (*this) =  A * (*this)
void  multiplyLeft(faust_mat<T> const& A){gemv(A, *this, *this, 1.0, 0.0, 'N');}
void  multiplyLeft(faust_spmat<T> const& A);
  
T sum()const{return vec.sum();}
T mean()const{return vec.mean();}


void operator=(faust_vec<T> const& y);

void operator*=(const T alpha);
void operator+=(const T alpha);
void operator-=(const T alpha);

void operator+=(const faust_vec<T>& v);
void operator-=(const faust_vec<T>& v);


T mean_relative_error(const faust_vec<T>& v);

T& operator[](faust_unsigned_int i){return vec(i);}
const T& operator[](faust_unsigned_int i)const{return vec(i);}

const T& operator()(faust_unsigned_int i)const{return vec(i);}


// friend algebra
friend void gemv<>(const faust_mat<T> & A,const faust_vec<T> & x,faust_vec<T> & y,const T & alpha, const T & beta, char typeA);
friend faust_vec<T> solve<>(const faust_mat<T> & A, const faust_vec<T> & v);
friend void sp_solve<>(const faust_mat<T> & A,faust_vec<T> & x, const faust_vec<T> & y);

private:
  faust_unsigned_int dim;
  static const char * class_name;
  Eigen::Matrix<T, Eigen::Dynamic,1> vec;

	#ifdef __COMPILE_TIMERS__
		public: 
			faust_timer t_local_multiplyLeft;
	#endif
  
};


#include "faust_vec.hpp"


#endif
