#ifndef __FAUST_PARAMS_H__
#define __FAUST_PARAMS_H__

#include "faust_constant.h"
#include <vector>
#ifdef __COMPILE_GPU__
   #include "faust_MatDense_gpu.h"
#else
   #include "faust_MatDense.h"
#endif
#include "StoppingCriterion.h"
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
            const Faust::MatDense<FPP,DEVICE>& data_,
            const unsigned int nb_fact_,
            const std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> & cons_,
            const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
            const StoppingCriterion<FPP>& stop_crit_2facts_ = StoppingCriterion<FPP>(defaultNiter1),
            const StoppingCriterion<FPP>& stop_crit_global_  = StoppingCriterion<FPP>(defaultNiter2),
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
        *   \param nb_fact_ : Number of factor used for the decomposition
        *	\param	cons_ : Specifies the constraint sets in which each factor should lie.<br>
                            It should be a std::vector<std::vector> of constraint_generic size 2*(nb_fact-1),<br>
                            where the jth columns sets the constraints for the jth factorization in two factors:<br>
                                - cons_[1][j] specifies the constraints for the left factor and<br>
                                - cons[2][j] for the right factor.<br>
        *   \param	init_fact_ : specifying the initial factors, could be an empty std::vector (not specifying the factors) <br>
        *	\param stop_crit_2facts : (optional) StoppingCriterion for each 2 factors factorization step <br>
        *	\param stop_crit_global : (optional) StoppingCriterion for the factorization each global factorization step <br>
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
            const Faust::MatDense<FPP,DEVICE>& data_,
            const unsigned int nb_fact_,
            const std::vector<std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> >& cons_,
            const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
            const StoppingCriterion<FPP>& stop_crit_2facts_ = StoppingCriterion<FPP>(defaultNiter1),
            const StoppingCriterion<FPP>& stop_crit_global_  = StoppingCriterion<FPP>(defaultNiter2),
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

        ~Params(){}


        public:
        // Required members
        // modif AL AL
        Faust::MatDense<FPP,DEVICE> data;

        faust_unsigned_int nb_fact; // number of factors
        std::vector<std::vector<const Faust::ConstraintGeneric<FPP,DEVICE> *> > cons; // vector of constraints
        std::vector<Faust::MatDense<FPP,DEVICE> > init_fact;

        // Optional members (set to default values if not defined)
        StoppingCriterion<FPP> stop_crit_2facts;
        StoppingCriterion<FPP> stop_crit_global;
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

        void Display() const;

        //const int nb_rows; // number of rows of the first factor
        //const int nb_cols; // number of columns of the last factor

        /*const int nb_it;   // number of iterations
        // if isStoppingCriterionError then criterion is error else criterion is number of iteration
        bool  isStoppingCriterionError;
        const faust_real errorThreshold;
        // only used as stopping criterion, if isStoppingCriterionError, when error is still greater than
        int maxIteration;*/
        static const char* class_name;
        private :

    };

}

#include "faust_Params.hpp"

#endif
