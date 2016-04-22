#ifndef __FAUST_PARAMS_H__
#define __FAUST_PARAMS_H__

#include "faust_constant.h"
#include <vector>
#ifdef __COMPILE_GPU__
   #include "faust_cu_mat.h"
#else
   #include "faust_mat.h"
#endif
#include "stopping_criterion.h"
#include "faust_constraint_generic.h"

#ifdef __COMPILE_GPU__
   template<typename T> class faust_cu_mat;
#else
   template<typename T> class faust_mat;
#endif


/*! \class faust_params
   * \brief template class representing the parameters for building hierarchical_fact object
   *\tparam T scalar numeric type, e.g float or double
   */
template<typename T>
class faust_params
{
#ifdef __COMPILE_GPU__
    typedef faust_cu_mat<T> faust_matrix ;
#else
    typedef faust_mat<T> faust_matrix ;
#endif

    public:

	faust_params(
        const faust_matrix& data_,
        const unsigned int nb_fact_,
        const std::vector<const faust_constraint_generic<T>*> & cons_,
        const std::vector<faust_matrix >& init_fact_,
        const stopping_criterion<T>& stop_crit_2facts_ = stopping_criterion<T>(defaultNiter1),
        const stopping_criterion<T>& stop_crit_global_  = stopping_criterion<T>(defaultNiter2),
        const T residuum_decrease_speed = defaultDecreaseSpeed,
        const T residuum_prcent = defaultResiduumPercent,
        const bool isVerbose_ = defaultVerbosity ,
        const bool isUpdateWayR2L_ = defaultUpdateWayR2L ,
        const bool isFactSideLeft_ = defaultFactSideLeft ,
        const T init_lambda_ = defaultLambda ,
        const bool constant_step_size_ = defaultConstantStepSize,
        const T step_size_ = defaultStepSize);


	/*!
    *   \brief faust_params constructor
    *   \param data : faust_matrix to hierarchically factorize
	*   \param nb_fact_ : Number of factor used for the decomposition
    *	\param	cons_ : Specifies the constraint sets in which each factor should lie.<br>
                        It should be a std::vector<std::vector> of constraint_generic size 2*(nb_fact-1),<br>
                        where the jth columns sets the constraints for the jth factorization in two factors:<br>
                            - cons_[1][j] specifies the constraints for the left factor and<br>
                            - cons[2][j] for the right factor.<br>
    *   \param	init_fact_ : specifying the initial factors, could be an empty std::vector (not specifying the factors) <br>
	*	\param stop_crit_2facts : (optional) stopping_criterion for each 2 factors factorization step <br>
	*	\param stop_crit_global : (optional) stopping_criterion for the factorization each global factorization step <br>
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
    faust_params(
        const faust_matrix& data_,
        const unsigned int nb_fact_,
        const std::vector<std::vector<const faust_constraint_generic<T>*> >& cons_,
        const std::vector<faust_matrix >& init_fact_,
        const stopping_criterion<T>& stop_crit_2facts_ = stopping_criterion<T>(defaultNiter1),
        const stopping_criterion<T>& stop_crit_global_  = stopping_criterion<T>(defaultNiter2),
        const bool isVerbose_ = defaultVerbosity ,
        const bool isUpdateWayR2L_ = defaultUpdateWayR2L ,
        const bool isFactSideLeft_ = defaultFactSideLeft ,
        const T init_lambda_ = defaultLambda ,
        const bool constant_step_size_ = defaultConstantStepSize,
        const T step_size_ = defaultStepSize);

    void init_from_file(const char* filename);
	faust_params();

    void check_constraint_validity();
	void check_bool_validity();

    ~faust_params(){}

  public:
    // Required members
    faust_matrix data;
    faust_unsigned_int nb_fact; // number of factors
    std::vector<std::vector<const faust_constraint_generic<T> *> > cons; // vector of constraints
    std::vector<faust_matrix > init_fact;

    // Optional members (set to default values if not defined)
    stopping_criterion<T> stop_crit_2facts;
    stopping_criterion<T> stop_crit_global;
    bool isVerbose;
    bool isUpdateWayR2L;
    bool isFactSideLeft;
    T init_lambda;
	bool isConstantStepSize;
	T step_size;

	//default value
	static const int defaultNiter1;
	static const int defaultNiter2;
	static const bool defaultVerbosity;
	static const bool defaultFactSideLeft;
	static const bool defaultUpdateWayR2L;
	static const T defaultLambda;
	static const bool defaultConstantStepSize;
	static const T defaultStepSize;
	static const T defaultDecreaseSpeed;
	static const T defaultResiduumPercent;

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

#include "faust_params.hpp"

#endif
