/****************************************************************************/
/*                              Description:                                */
/*  C++ Hpp file wrapping the  Class Faust::Transform,                      */
/*  the methods are quite the same but this works with pointer              */
/*  which are easier to handle from Python                                  */
/*                                                                          */     
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/

#include "faust_Transform.h"
#include <iostream>
#include <exception>

    template<typename FPP>
void FaustCoreCpp<FPP>::push_back(FPP* valueMat, unsigned int nbrow,unsigned int nbcol)
{
    Faust::MatDense<FPP,Cpu> dense_mat(valueMat,nbrow,nbcol);
    Faust::MatSparse<FPP,Cpu> sparse_mat(dense_mat);
    //sparse_mat.Display();
    this->transform.push_back(&sparse_mat);




}

template<typename FPP>
void FaustCoreCpp<FPP>::multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x)const
{



    unsigned int nbRowThis,nbColThis;


    nbRowThis = this->getNbRow();
    nbColThis = this->getNbCol();

    if ( (nbrow_y != nbRowThis) | (nbrow_x != nbColThis) | (nbcol_y != nbcol_x) )
    {	
        std::cout<<"nbRowThis "<<nbRowThis<<" must be equal to nb row y  "<<nbrow_y<<std::endl;
        std::cout<<"nbColThis "<<nbColThis<<" must be equal to nb row x  "<<nbrow_x<<std::endl;
        std::cout<<"nbcol_y "<<nbcol_y<<" must be equal to nbcol_x  "<<nbcol_x<<std::endl;
        handleError("FaustCpp"," multiply : invalid dimension");
    }
    if (nbcol_x == 1)
    {
        Faust::Vect<FPP,Cpu> X(nbrow_x,value_x);
        Faust::Vect<FPP,Cpu> Y;


        Y = this->transform.multiply(X);

        memcpy(value_y,Y.getData(),sizeof(FPP)*nbrow_y);
    }else
    {
        Faust::MatDense<FPP,Cpu> X(value_x,nbrow_x,nbcol_x);
        Faust::MatDense<FPP,Cpu> Y;

        Y = this->transform.multiply(X);


        memcpy(value_y,Y.getData(),sizeof(FPP)*nbrow_y*nbcol_y);
    }


}


template<typename FPP>
unsigned long long FaustCoreCpp<FPP>::nnz() const
{
    return this->transform.get_total_nnz();
}


template<typename FPP>
double FaustCoreCpp<FPP>::norm(int ord) const
{
    double precision = 0.001;
    faust_unsigned_int nbr_iter_max = 100;
    int flag; //not used yet
    switch(ord) {
        case 2:
            return this->transform.spectralNorm(nbr_iter_max, precision, flag);
        case 1:
            return this->transform.normL1();
        default:
            handleError("FaustCoreCpp", "norm(int ord) unvalid norm order.");
    }
    return -1;
}

template<typename FPP>
double FaustCoreCpp<FPP>::normFro() const
{
    return this->transform.normFro();
}

template<typename FPP>
double FaustCoreCpp<FPP>::get_nb_factors() const
{
    double nb_fact = this->transform.size();
    return nb_fact;
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::get_fact_nb_rows(unsigned int& i) const
{
    Faust::MatDense<FPP,Cpu> factor_generic = this->transform.get_fact(i);
    unsigned int nb_rows = factor_generic.getNbRow();
    return nb_rows;
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::get_fact_nb_cols(unsigned int& i) const
{
    Faust::MatDense<FPP,Cpu> factor_generic = this->transform.get_fact(i);
    unsigned int nb_cols = factor_generic.getNbCol();
    return nb_cols;
}

template<typename FPP>
void FaustCoreCpp<FPP>::get_fact(unsigned int& i, FPP* fact_ptr) const
{
    Faust::MatDense<FPP,Cpu> dense_factor = this->transform.get_fact(i);
    //TODO: optimize here (we have two copies from C++ object to Py, the first in MatDense::Clone()
    //the second here)
    memcpy(fact_ptr, dense_factor.getData(),
            sizeof(FPP)*dense_factor.getNbCol()*dense_factor.getNbRow());

}

template<typename FPP>
void FaustCoreCpp<FPP>::save_mat_file(const char* filepath) const
{
    //    std::cout << "FaustCoreCpp::save_mat_file()" << std::endl;
    this->transform.save_mat_file(filepath);
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::getNbRow() const
{
    return this->transform.getNbRow();
}

template<typename FPP>
unsigned int FaustCoreCpp<FPP>::getNbCol() const
{
    return this->transform.getNbCol();
}

    template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::transpose()
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform.transpose();
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>();
    core->transform = th;
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::transpose() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

    template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::conjugate()
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform.conjugate();
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>();
    core->transform = th;
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::conjugate() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}

template<typename FPP>
  FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::randFaust(unsigned int t,
            unsigned int min_num_factors, unsigned int max_num_factors,
            unsigned int min_dim_size, unsigned int max_dim_size, float density) {
      Faust::TransformHelper<FPP,Cpu>* th = Faust::TransformHelper<FPP,Cpu>::randFaust(Faust::RandFaustType(t), min_num_factors, max_num_factors, min_dim_size, max_dim_size, density);
      FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>();
      core->transform = th;
      return core;
  }

    template<typename FPP>
FaustCoreCpp<FPP>* FaustCoreCpp<FPP>::adjoint()
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform.adjoint();
    FaustCoreCpp<FPP>* core = new FaustCoreCpp<FPP>();
    core->transform = th;
#ifdef FAUST_VERBOSE
    std::cout << "FaustCoreCpp::adjoint() th=" << th << "core=" << core << std::endl;
#endif
    return core;
}
