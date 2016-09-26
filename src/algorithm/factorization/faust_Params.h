/****************************************************************************/
/*                              Description:                                */
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
#ifndef __FAUST_PARAMS_H__
#define __FAUST_PARAMS_H__

#include "faust_constant.h"
#include <vector>
#ifdef __COMPILE_GPU__
   #include "faust_MatDense_gpu.h"
#else
   #include "faust_MatDense.h"
#endif
#include "faust_StoppingCriterion.h"
#include "faust_ConstraintGeneric.h"

/*! \class Faust::Params
* \brief template class representing the parameters for building HierarchicalFact object
*\tparam FPP scalar numeric type, e.g float or double
*/

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    template<typename FPP,Device DEVICE> class MatDense;


    template<typename FPP,Device DEVICE>
    class Params
    {


        public:

        Params(
            const faust_unsigned_int nbRow,
            const faust_unsigned_int nbCol, 
            const unsigned int nbFact_,
            const std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> & cons_,
            const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
            const Faust::StoppingCriterion<FPP>& stop_crit_2facts_ = Faust::StoppingCriterion<FPP>(defaultNiter1),
            const Faust::StoppingCriterion<FPP>& stop_crit_global_  = Faust::StoppingCriterion<FPP>(defaultNiter2),
            const FPP residuum_decrease_speed = defaultDecreaseSpeed,
            const FPP residuum_prcent = defaultResiduumPercent,
            const bool isVerbose_ = defaultVerbosity ,
            const bool isUpdateWayR2L_ = defaultUpdateWayR2L ,
            const bool isFactSideLeft_ = defaultFactSideLeft ,
            const FPP init_lambda_ = defaultLambda ,
            const bool constant_step_size_ = defaultConstantStepSize,
            const FPP step_size_ = defaultStepSize);


        /*!
        *   \brief Faust::Params constructor
        *   \param data : Faust::MatDense<FPP,DEVICE> to hierarchically factorize
        *   \param nbFact_ : Number of factor used for the decomposition
        *	\param	cons_ : Specifies the constraint sets in which each factor should lie.<br>
                            It should be a std::vector<std::vector> of constraint_generic size 2*(nbFact-1),<br>
                            where the jth columns sets the constraints for the jth factorization in two factors:<br>
                                - cons_[1][j] specifies the constraints for the left factor and<br>
                                - cons[2][j] for the right factor.<br>
        *   \param	init_fact_ : specifying the initial factors, could be an empty std::vector (not specifying the factors) <br>
        *	\param stop_crit_2facts : (optional) Faust::StoppingCriterion for each 2 factors factorization step <br>
        *	\param stop_crit_global : (optional) Faust::StoppingCriterion for the factorization each global factorization step <br>
        *	\param isVerbose : (optional) - if true the function outputs the error at each iteration <br>
                                           - if false, hierarchical_fact run silent mode<br>
                                           (default value is false) <br>
        *   \param isUpdateWayR2L_ : (optional)    - if true, the factors are updated from right to left, <br>
                                                    - if false, the factors are updated from left to right.<br>
                                                    (the default value is false) <br>
        *	\param isFactSideLeft_ : (optional) Side to be factorized iteratively:
                                                    - if true, the left <br>
                                                    - if false, the right <br>
        *
        *   The following figure shows the functionalitie of parameters isUpdateWayR2L and isFactSideLeft
        *   <img src="../../doc/html/factSide.jpg" alt="FactSide" width=600px /> <br>
        *   <br>
        *	\param init_lambda_ : (optional) specifies a starting point for the algorithm.<br>
                                The default value is 1.0 <br>
        *	\param constant_step_size : (optional) specifies if the stepsize of the gradient descent is constant.<br>
                                The default value is false. In this case, the stepsize is controlled by the lipschitz modulus of the gradient <br>

        *   \param step_size : (optional) specifies the step size of the gradient descent, USELESS if constant_step_size is false. default value : 1e-16 <br>
        */
        Params(
            const faust_unsigned_int nbRow_,
            const faust_unsigned_int nbCol_,  		
            const unsigned int nbFact_,
            const std::vector<std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> >& cons_,
            const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
            const Faust::StoppingCriterion<FPP>& stop_crit_2facts_ = Faust::StoppingCriterion<FPP>(defaultNiter1),
            const Faust::StoppingCriterion<FPP>& stop_crit_global_  = Faust::StoppingCriterion<FPP>(defaultNiter2),
            const bool isVerbose_ = defaultVerbosity ,
            const bool isUpdateWayR2L_ = defaultUpdateWayR2L ,
            const bool isFactSideLeft_ = defaultFactSideLeft ,
            const FPP init_lambda_ = defaultLambda ,
            const bool constant_step_size_ = defaultConstantStepSize,
            const FPP step_size_ = defaultStepSize);

        void init_from_file(const char* filename);
        Params();

        void check_constraint_validity();
        void check_bool_validity();

        void Display() const;
        ~Params(){}


        public:
        // Required members
        
	// data is now independant from the params class,
        //Faust::MatDense<FPP,DEVICE> data; 
	faust_unsigned_int m_nbRow; // number of row of the matrix to be factorized
	faust_unsigned_int m_nbCol; // number of columns of the matrix to be factorized

        faust_unsigned_int m_nbFact; // number of factors
        std::vector<std::vector<const Faust::ConstraintGeneric<FPP,DEVICE> *> > cons; // vector of constraints
        std::vector<Faust::MatDense<FPP,DEVICE> > init_fact;

        // Optional members (set to default values if not defined)
        Faust::StoppingCriterion<FPP> stop_crit_2facts;
        Faust::StoppingCriterion<FPP> stop_crit_global;
        bool isVerbose;
        bool isUpdateWayR2L;
        bool isFactSideLeft;
        FPP init_lambda;
        bool isConstantStepSize;
        FPP step_size;

        //default value
        static const int defaultNiter1;
        static const int defaultNiter2;
        static const bool defaultVerbosity;
        static const bool defaultFactSideLeft;
        static const bool defaultUpdateWayR2L;
        static const FPP defaultLambda;
        static const bool defaultConstantStepSize;
        static const FPP defaultStepSize;
        static const FPP defaultDecreaseSpeed;
        static const FPP defaultResiduumPercent;

        //const int nb_rows; // number of rows of the first factor
        //const int nb_cols; // number of columns of the last factor

        /*const int nb_it;   // number of iterations
        // if isFaust::StoppingCriterionError then criterion is error else criterion is number of iteration
        bool  isFaust::StoppingCriterionError;
        const faust_real errorThreshold;
        // only used as stopping criterion, if isFaust::StoppingCriterionError, when error is still greater than
        int maxIteration;*/

        // modif AL AL ???
        //static const char* m_className;
        private :
            static const char* m_className;

    };

}

#include "faust_Params.hpp"

#endif
