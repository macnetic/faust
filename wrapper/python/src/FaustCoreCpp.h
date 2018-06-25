/****************************************************************************/
/*                              Description:                                */
/*  C++ Header file wrapping the  Class Faust::Transform<FPP,Cpu>,          */
/*  the methods are quite the same but this works with pointer              */
/*  which are easier to handle from Python                                  */
/*                                                                          */     
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
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

#ifndef FAUSTCORECPP_H
#define FAUSTCORECPP_H

//#include "faust_transform.h"
#include "faust_MatDense.h"
#include "faust_TransformHelper.h"

template<typename FPP>
class FaustCoreCpp
{

    public :


    //FaustCoreCpp(): transform(){}
    void Display() const { transform.display();}
    void push_back(FPP* valueMat,unsigned int nbrow,unsigned int nbcol);
    unsigned int getNbRow() const;
    unsigned int getNbCol() const;
    void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x/*,bool isTranspose*/)const;
    unsigned long long nnz()const;
    double norm(int ord) const;
    double normFro() const;
    double get_nb_factors() const;
    unsigned int get_fact_nb_rows(unsigned int& i) const;
    unsigned int get_fact_nb_cols(unsigned int& i) const;
    void get_fact(unsigned int& i, FPP* fact_ptr) const;
    void save_mat_file(const char* filepath) const;
    FaustCoreCpp<FPP>* transpose();
    FaustCoreCpp<FPP>* conjugate();
    FaustCoreCpp<FPP>* adjoint();
    static FaustCoreCpp<FPP>* randFaust(unsigned int t,
            unsigned int min_num_factors, unsigned int max_num_factors,
            unsigned int min_dim_size,
            unsigned int max_dim_size, float density);
    private :
    Faust::TransformHelper<FPP,Cpu> transform;
};

#include "FaustCoreCpp.hpp"

#endif
