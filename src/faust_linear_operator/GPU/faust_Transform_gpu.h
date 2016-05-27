#ifndef __FAUST_Transform_GPU_H__
#define __FAUST_Transform_GPU_H__

#include <vector>
#include "faust_MatSparse_gpu.h"







//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{


template<typename FPP,Device DEVICE> class MatDense;
template<typename FPP,Device DEVICE> class MatSparse;
template<typename FPP,Device DEVICE> class Vect;
template<typename FPP,Device DEVICE> class Transform;


    template<typename FPP>
    class Transform<FPP,Gpu>
    {
        public:
            Transform();
            Transform(const std::vector<Faust::MatSparse<FPP,Gpu> >& facts, const FPP lambda_ = (FPP)1.0);
            Transform(const Transform<FPP,Gpu> & A);
            Transform(const std::vector<Faust::MatDense<FPP,Gpu> >&facts);
            void get_facts(std::vector<Faust::MatSparse<FPP,Gpu> >& sparse_facts)const{sparse_facts = data;}
            void get_facts(std::vector<Faust::MatDense<FPP,Gpu> >& facts)const;
            int size()const{return data.size();}

            Faust::MatDense<FPP,Gpu> get_product(BlasHandle<Gpu>,SpBlasHandle<Gpu>)const;
            Faust::MatSparse<FPP,Gpu> get_fact(int id) const;
            int getNbRow() const;
            int getNbCol() const;
            void print_file(const char* filename) const;
            void init_from_file(const char* filename);
            long long int get_total_nnz()const{return totalNonZeros;}
            void clear(){data.resize(0);totalNonZeros=0;}
            void push_back(const Faust::MatSparse<FPP,Gpu>& S);
            void mult( const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y,SpBlasHandle<Gpu> spblas_Handle);
            void push_first(const Faust::MatSparse<FPP,Gpu>& S);
            void pop_back(Faust::MatSparse<FPP,Gpu>& S);
            void pop_first(Faust::MatSparse<FPP,Gpu>& S);
            void pop_first(Faust::MatSparse<FPP,Gpu>& S) const;
            void Display()const;
            void transpose();
            void updateNonZeros();
            //(*this) = (*this) * A
            void multiply(const Transform<FPP,Gpu> & A);
            //(*this) = A * (*this)
            void multiplyLeft(const Transform<FPP,Gpu> & A);
            void scalarMultiply(const FPP scalar);
            FPP spectralNorm(const int nbr_iter_max, FPP threshold, int &flag) const;
            ~Transform(){}


        public:
            void operator=(const Transform<FPP,Gpu>&  f){data=f.data;totalNonZeros=f.totalNonZeros;}
            // add all of the sparse matrices from f.data to this->data
            void operator*=(const FPP  scalar){scalarMultiply(scalar);};
            void operator*=(const Transform<FPP,Gpu>&  f){multiply(f);};
            // add the sparse matrix S to this->data
            void operator*=(const Faust::MatSparse<FPP,Gpu>&  S){push_back(S);totalNonZeros+=S.getNonZeros();}




        private:
            std::vector<Faust::MatSparse<FPP,Gpu> > data;
            long long int totalNonZeros;
            static const char * class_name;



    };

}


#include "faust_Transform_gpu.hpp"

#endif
