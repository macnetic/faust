#ifndef __PALM4MSA_HPP__
#define __PALM4MSA_HPP__

//#include "palm4MSA.h"

#include <sstream>

#include <iostream>
#include "faust_params.h"
#include "faust_params_palm.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_real.h"
#include "faust_constraint_int.h"
#include "faust_constraint_mat.h"
#include "faust_mat.h"

#include "LinAlgebra.h"
#include "prox.h"
#include "faust_constraint_type.h"
#include "faust_exception.h"

#ifdef __COMPILE_TIMERS__
#include "faust_timer.h"
#endif

//#include"faust_init_from_matio_mat.h"

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#define __SP setprecision(20)<<

using namespace std;

template<typename T>
const char * palm4MSA<T>::class_name="palm4MSA";

template<typename T>
const T palm4MSA<T>::lipschitz_multiplicator=1.001;



template<typename T>
palm4MSA<T>::palm4MSA(const faust_params<T> & params_, const bool isGlobal_) :
    data(params_.data),
    lambda(params_.init_lambda),
    nb_fact(0),
    S(params_.init_fact),
    RorL(vector<faust_mat<T> >(2)),
    ind_fact(0),
    ind_ite(-1),
    verbose(params_.isVerbose),
    isUpdateWayR2L(params_.isUpdateWayR2L),
    isConstantStepSize(params_.isConstantStepSize),
    isGradComputed(false),
    isProjectionComputed(false),
    isLastFact(false),
    isConstraintSet(false),
    isGlobal(isGlobal_),
    isInit(false),
    c(1/params_.step_size)
{
    if(isGlobal)
        stop_crit = stopping_criterion<T>(params_.stop_crit_global);
    else
        stop_crit = stopping_criterion<T>(params_.stop_crit_2facts);

	if (isConstantStepSize)
        isCComputed = true;
	else
        isCComputed = false;

}

template<typename T>
palm4MSA<T>::palm4MSA(const faust_params_palm<T>& params_palm_,
const bool isGlobal_/*=false*/) :
   stop_crit(params_palm_.stop_crit),
   data(params_palm_.data),
   lambda(params_palm_.init_lambda),
   nb_fact(params_palm_.nb_fact),
   S(params_palm_.init_fact),
   RorL(vector<faust_mat<T> >(2)),
   LorR(faust_mat<T>(params_palm_.init_fact[0].getNbRow())),
   const_vec(params_palm_.cons),
   ind_fact(0),
   ind_ite(-1),
   verbose(params_palm_.isVerbose),
   isUpdateWayR2L(params_palm_.isUpdateWayR2L),
   isConstantStepSize(params_palm_.isConstantStepSize),
   isGradComputed(false),
   isProjectionComputed(false),
   isLastFact(false),
   isConstraintSet(false),
   isGlobal(isGlobal_),
   c(1/params_palm_.step_size)
{
   	if (isConstantStepSize)
		isCComputed = true;
	else
		isCComputed = false;

   check_constraint_validity();

}

template<typename T>
void palm4MSA<T>::compute_facts()
{
	while (do_continue())
	{
		next_step();
	}

}

template<typename T>
void palm4MSA<T>::get_facts(faust_core<T> & faust_fact) const
{
	faust_core<T> f(S);
	faust_fact = f;


}


template<typename T>
void palm4MSA<T>::compute_projection()
{
#ifdef __COMPILE_TIMERS__
t_global_compute_projection.start();
t_local_compute_projection.start();
#endif

   
      //faust_mat<T> matrix2project(S[ind_fact]);
      S[ind_fact] -= grad_over_c;
      const_vec[ind_fact]->project(S[ind_fact]);	


    isProjectionComputed = true;
#ifdef __COMPILE_TIMERS__
t_global_compute_projection.stop();
t_local_compute_projection.stop();
#endif
}

template<typename T>
void palm4MSA<T>::compute_grad_over_c()
{
/*static int cmpt = -1;
cmpt++;
char nomFichier[100];*/

#ifdef __COMPILE_TIMERS__
t_global_compute_grad_over_c.start();
t_local_compute_grad_over_c.start();
#endif
    if(!isCComputed)
    {
        throw std::logic_error("c must be set before computing grad/c");
    }

/*! \brief There are 4 ways to compute gradient : <br>
* (0) : lambda*(L'*(lambda*(L*S)*R - X))*R' : complexity = L1*L2*S2 + L1*S2*R2 + L2*L1*R2 + L2*R2*R2; <br>
* (1) : lambda*L'*((lambda*(L*S)*R - X)*R') : complexity = L1*L2*S2 + L1*S2*R2 + L1*R2*S2 + L2*L1*S2; <br>
* (2) : lambda*(L'*(lambda*L*(S*R) - X))*R' : complexity = L2*S2*R2 + L1*L2*R2 + L2*L1*R2 + L2*R2*S2; <br>
* (3) : lambda*L'*((lambda*L*(S*R) - X)*R') : complexity = L2*S2*R2 + L1*L2*R2 + L1*R2*S2 + L2*L1*S2; <br>
*  with L of size L1xL2 <br>
*       S of size L2xS2 <br>
*       R of size S2xR2 <br>
*/
    unsigned long long int L1, L2, R2, S2;
    if (!isUpdateWayR2L)
    {
        L1 = (unsigned long long int) LorR.getNbRow();
        L2 = (unsigned long long int) LorR.getNbCol();
        R2 = (unsigned long long int) RorL[ind_fact].getNbCol();
    }
    else
    {
        L1 = (unsigned long long int) RorL[ind_fact].getNbRow();
        L2 = (unsigned long long int) RorL[ind_fact].getNbCol();
        R2 = (unsigned long long int) LorR.getNbCol();
    }
    S2 = (unsigned long long int) S[ind_fact].getNbCol();
    vector<unsigned long long int > complexity(4,0);
    complexity[0] = L1*L2*S2 + L1*S2*R2 + L2*L1*R2 + L2*R2*R2;
    complexity[1] = L1*L2*S2 + L1*S2*R2 + L1*R2*S2 + L2*L1*S2;
    complexity[2] = L2*S2*R2 + L1*L2*R2 + L2*L1*R2 + L2*R2*S2;
    complexity[3] = L2*S2*R2 + L1*L2*R2 + L1*R2*S2 + L2*L1*S2;

    int idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));


    error = data;
    faust_mat<T> tmp1,tmp3;

    if (idx==0 || idx==1) // computing L*S first, then (L*S)*R
    {
        if (!isUpdateWayR2L)
        {
            // tmp1 = L*S
            multiply(LorR, S[ind_fact], tmp1);

/*sprintf(nomFichier,"LorR_0_%d_host.tmp",cmpt);
LorR.print_file(nomFichier);
sprintf(nomFichier,"RorL%d_0_%d_host.tmp",ind_fact,cmpt);
RorL[ind_fact].print_file(nomFichier);
sprintf(nomFichier,"S%d_0_%d_host.tmp",ind_fact,cmpt);
S[ind_fact].print_file(nomFichier);
sprintf(nomFichier,"tmp1_0_%d_host.tmp",cmpt);
tmp1.print_file(nomFichier);
sprintf(nomFichier,"error_0_%d_host.tmp",cmpt);
error.print_file(nomFichier);
cout << "appel " << cmpt<<" : lambda0 = "<< lambda<<endl;*/
            // error = lambda*tmp1*R - error (= lambda*L*S*R - data )
            gemm<T>(tmp1, RorL[ind_fact], error, lambda, -1.0, 'N', 'N');

/*sprintf(nomFichier,"LorR_1_%d_host.tmp",cmpt);
LorR.print_file(nomFichier);
sprintf(nomFichier,"S%d_1_%d_host.tmp",ind_fact,cmpt);
S[ind_fact].print_file(nomFichier);
sprintf(nomFichier,"tmp1_1_%d_host.tmp",cmpt);
tmp1.print_file(nomFichier);
sprintf(nomFichier,"RorL%d_1_%d_host.tmp",ind_fact,cmpt);
RorL[ind_fact].print_file(nomFichier);
sprintf(nomFichier,"error_1_%d_host.tmp",cmpt);
error.print_file(nomFichier);*/
        }
        else
        {
            // tmp1 = L*S
            multiply(RorL[ind_fact], S[ind_fact], tmp1);
            // error = lambda*tmp1*R - error (= lambda*L*S*R - data )
            gemm<T>(tmp1, LorR, error, lambda, -1.0, 'N', 'N');
        }
    }
    else // computing S*R first, then L*(S*R)
    {
        if (!isUpdateWayR2L)
        {
            // tmp1 = S*R
            multiply(S[ind_fact], RorL[ind_fact], tmp1);

            // error = lambda*L*tmp1 - error (= lambda*L*S*R - data )
            gemm<T>(LorR, tmp1, error, lambda, -1.0, 'N', 'N');

        }
        else
        {
            // tmp1 = S*R
            multiply(S[ind_fact], LorR, tmp1);
            // error = lambda*L*tmp1 - error (= lambda*L*S*R - data )
            gemm<T>(RorL[ind_fact], tmp1, error, lambda, -1.0, 'N', 'N');
        }
    }

    if (idx==0 || idx==2) // computing L'*error first, then (L'*error)*R'
    {
        if (!isUpdateWayR2L)
        {
            // tmp3 = lambda*L'*error (= lambda*L' * (lambda*L*S*R - data) )

            gemm<T>(LorR, error, tmp3, lambda, 0.0, 'T', 'N');

            // grad_over_c = 1/c*tmp3*R' (= 1/c*lambda*L' * (lambda*L*S*R - data) * R' )
            gemm<T>(tmp3, RorL[ind_fact], grad_over_c, 1.0/c, 0.0,'N','T');

        }
        else
        {
            // tmp3 = lambda*L'*error (= lambda*L' * (lambda*L*S*R - data) )
            gemm<T>(RorL[ind_fact], error, tmp3, lambda, 0.0, 'T', 'N');
            // grad_over_c = 1/c*tmp3*R' (= 1/c*lambda*L' * (lambda*L*S*R - data) * R' )
            gemm<T>(tmp3, LorR, grad_over_c, 1.0/c, 0.0,'N','T');
        }
    }
    else // computing error*R' first, then L'*(error*R')
    {
        if (!isUpdateWayR2L)
        {
            // tmp3 = lambda*error*R' (= lambda*(lambda*L*S*R - data) * R' )
            gemm<T>(error, RorL[ind_fact], tmp3, lambda, 0.0, 'N', 'T');
            // grad_over_c = 1/c*L'*tmp3 (= 1/c*L' * lambda*(lambda*L*S*R - data) * R' )
            gemm<T>(LorR, tmp3, grad_over_c,1.0/c, 0.0,'T','N');
        }
        else
        {
            // tmp3 = lambda*error*R' (= lambda * (lambda*L*S*R - data) * R' )
            gemm<T>(error, LorR, tmp3, lambda, 0.0, 'N', 'T');
            // grad_over_c = 1/c*L'*tmp3 (= 1/c*L' * lambda*(lambda*L*S*R - data) * R' )
            gemm<T>(RorL[ind_fact], tmp3, grad_over_c, 1.0/c, 0.0,'T','N');
        }

    }

    isGradComputed = true;

#ifdef __COMPILE_TIMERS__
t_global_compute_grad_over_c.stop();
t_local_compute_grad_over_c.stop();
#endif
}



template<typename T>
void palm4MSA<T>::compute_lambda()
{
#ifdef __COMPILE_TIMERS__
t_global_compute_lambda.start();
t_local_compute_lambda.start();
#endif

   if (!isLastFact)
   {
     handleError(class_name,"compute_lambda : computation of lambda must be done at the end of the iteration through the number of factors");
   }

   // As LorR has also been updated at the end of the last iteration over the facts, LorR matches X_hat, which the product of all factors, including the last one.
   // Xt_Xhat = data'*X_hat
   faust_mat<T> Xt_Xhat;
   gemm<T>(data, LorR, Xt_Xhat, 1.0, 0.0, 'T','N');

   // Xhatt_Xhat = X_hat'*X_hat
   faust_mat<T> Xhatt_Xhat;
   
   gemm<T>(LorR, LorR, Xhatt_Xhat, 1.0, 0.0, 'T','N');
   
   T Xhatt_Xhat_tr = (T) Xhatt_Xhat.trace();
   
   if (Xhatt_Xhat_tr != 0)
   {
		lambda = Xt_Xhat.trace()/Xhatt_Xhat_tr;
		if (isnan(lambda))
		{
			handleError(class_name,"compute_lambda : Xhatt_Xhat_tr is too small or Xt_Xhat.trace is too big so lambda is infinite");
		}
	}else
	{
		handleError(class_name,"compute_lambda : Xhatt_Xhat_tr equal 0 so lambda is infinite");
	}
   //cout<<"lambda : "<<lambda<<endl;
   //cout<<__SP lambda<<endl;

#ifdef __COMPILE_TIMERS__
t_global_compute_lambda.stop();
t_local_compute_lambda.stop();
#endif
}


template<typename T>
void palm4MSA<T>::update_R()
{
#ifdef __COMPILE_TIMERS__
t_global_update_R.start();
t_local_update_R.start();
#endif

   // R[nb_fact-1] est initialise a l'identite lors de la creation de l'objet palm4MSA et n'est pas cense changer
   if (!isUpdateWayR2L)
   {
      RorL.resize(nb_fact);
      RorL[nb_fact-1].resize(const_vec[nb_fact-1]->getCols());
      RorL[nb_fact-1].setEyes();


      for (int i=nb_fact-2 ; i>-1 ; i--)
         //  R[i] = S[i+1] * R[i+1]
         multiply(S[i+1], RorL[i+1], RorL[i]);

   }
   else
   {
      if(!isProjectionComputed)
      {
         throw std::logic_error("Projection must be computed before updating L");
      }
      LorR.multiplyLeft(S[ind_fact]);
   }



#ifdef __COMPILE_TIMERS__
t_global_update_R.stop();
t_local_update_R.stop();
#endif
}


template<typename T>
void palm4MSA<T>::update_L()
{


#ifdef __COMPILE_TIMERS__
t_global_update_L.start();
t_local_update_L.start();
#endif

   if(!isUpdateWayR2L)
   {
      if(!isProjectionComputed)
      {
         throw std::logic_error("Projection must be computed before updating L");
      }

      LorR *= S[ind_fact];

	}
   else
   {
      RorL.resize(nb_fact);
      RorL[0].resize(const_vec[0]->getRows());
      RorL[0].setEyes();
      for (int i=0 ; i<nb_fact-1 ; i++)
         //  R[i] = S[i+1] * R[i+1]
         multiply(RorL[i] , S[i], RorL[i+1]);
   }
#ifdef __COMPILE_TIMERS__
t_global_update_L.stop();
t_local_update_L.stop();
#endif
}

template<typename T>
void palm4MSA<T>::check_constraint_validity()
{
#ifdef __COMPILE_TIMERS__
t_global_check.start();
t_local_check.start();
#endif

   if (nb_fact != S.size())
   {
      handleError(class_name," Wrong initialization: params.nfacts and params.init_facts are in conflict");
   }
#ifdef __COMPILE_TIMERS__
t_global_check.stop();
t_local_check.stop();
#endif
}


#ifdef __PAS_FIXE__
// palm4MSA<T>::compute_c() has been defined as an inline method in palm4MSA.h
#else
template<typename T>
void palm4MSA<T>::compute_c()
{
#ifdef __COMPILE_TIMERS__
	t_global_compute_c.start();
	t_local_compute_c.start();
#endif


   if (!isConstantStepSize)
   {
		faust_int flag1,flag2;

	   int nbr_iter = 10000;
	   T threshold = 1e-16;
	   T nL1=LorR.spectralNorm(nbr_iter,threshold,flag1);
	   T nR1=RorL[ind_fact].spectralNorm(nbr_iter,threshold,flag2);
		c=lipschitz_multiplicator*nR1*nR1*nL1*nL1*lambda*lambda;
   }

   isCComputed = true;


   #ifdef __COMPILE_TIMERS__
	t_global_compute_c.stop();
	t_local_compute_c.stop();
	#endif
}
#endif

template<typename T>
void palm4MSA<T>::init_fact(int nb_facts_)
{
#ifdef __COMPILE_TIMERS__
t_global_init_fact.start();
t_local_init_fact.start();
#endif

  /*if(isInit && isGlobal)
  {
     cerr << "Error in palm4MSA<T>::init_fact : global factorization has already been initialized" <<endl;
     exit(EXIT_FAILURE);
  }
  else if (isGlobal)
  {
     nb_fact = 1;
     S.resize(nb_fact);
     isInit = true;
     return;
  }*/



// if !isGlobal

  if(!isConstraintSet)
  {
     handleError(class_name,"init_fact : constrainst must be set before calling init_fact");
  }
   nb_fact = nb_facts_;
   S.resize(nb_fact);
   if (!isUpdateWayR2L)
   {
      S[0].resize(const_vec[0]->getRows(), const_vec[0]->getCols());
      S[0].setZeros();
      for (int i=1 ; i<nb_fact ; i++)
      {
         S[i].resize(const_vec[i]->getRows(), const_vec[i]->getCols());
         S[i].setEyes();
      }
   }
   else
   {
      for (int i=0 ; i<nb_fact-1 ; i++)
      {
         S[i].resize(const_vec[i]->getRows(), const_vec[i]->getCols());
         S[i].setEyes();
      }
      S[nb_fact-1].resize(const_vec[nb_fact-1]->getRows(), const_vec[nb_fact-1]->getCols());
      S[nb_fact-1].setZeros();
   }


#ifdef __COMPILE_TIMERS__
t_global_init_fact.stop();
t_local_init_fact.stop();
#endif
}

/** \brief
 *
 * \param
 * \param
 */
template<typename T>
void palm4MSA<T>::next_step()
{
#ifdef __COMPILE_TIMERS__
t_global_next_step.start();
t_local_next_step.start();
#endif

    check_constraint_validity();
    // resizing L or R
    if(!isUpdateWayR2L)
    {
        LorR.resize(const_vec[0]->getRows());
        LorR.setEyes();
        update_R();
    }
    else
    {
        LorR.resize(const_vec[nb_fact-1]->getCols());
        LorR.setEyes();
        update_L();
    }


    int* ind_ptr = new int[nb_fact];
    for (int j=0 ; j<nb_fact ; j++)
    if (!isUpdateWayR2L)
        ind_ptr[j] = j;
    else
        ind_ptr[j] = nb_fact-1-j;

    for (int j=0 ; j<nb_fact ; j++)
    {
        if (j == nb_fact-1)
            isLastFact = true;
        else
            isLastFact = false;

        ind_fact = ind_ptr[j];
        if (!isConstantStepSize)
            isCComputed = false;

        isGradComputed = false;
        isProjectionComputed = false;

        if (!isConstantStepSize)
            compute_c();

        compute_grad_over_c();
        compute_projection();

        if(!isUpdateWayR2L)
            update_L();
        else
            update_R();

    }

	
    compute_lambda();


	if (verbose)
    {
        cout << "Iter " << ind_ite << ", RMSE=" << get_RMSE() << endl;
        cout << "Lambda " <<setprecision(20)<< lambda << endl;
    }
    delete[] ind_ptr;
    ind_ptr = NULL;


//cout<<"lambda : "<< lambda<< endl;
#ifdef __COMPILE_TIMERS__
t_global_next_step.stop();
t_local_next_step.stop();
#endif
}

template<typename T>
void palm4MSA<T>::init_fact_from_palm(const palm4MSA& palm2, bool isFactSideLeft)
{
#ifdef __COMPILE_TIMERS__
t_global_init_fact_from_palm.start();
t_local_init_fact_from_palm.start();
#endif


    if (palm2.nb_fact != 2)
    {
        handleError(class_name,"init_fact_from_palm : argument palm2 must contain 2 factors.");
    }

    if(!isConstraintSet)
    {
        handleError(class_name,"init_fact_from_palm : constrainst must be set before calling init_fact_from_palm");
    }

    if(isFactSideLeft)
    {
        S.insert(S.begin(), palm2.S[0]);
        S[1] = palm2.S[1];
    }
    else
    {
        S[S.size()-1] = palm2.S[0];
        S.push_back(palm2.S[1]);
    }
    nb_fact++;

    check_constraint_validity();

#ifdef __COMPILE_TIMERS__
t_global_init_fact_from_palm.stop();
t_local_init_fact_from_palm.stop();
#endif
}


#ifdef __COMPILE_TIMERS__

template<typename T> faust_timer palm4MSA<T>::t_global_compute_projection;
template<typename T> faust_timer palm4MSA<T>::t_global_compute_grad_over_c;
template<typename T> faust_timer palm4MSA<T>::t_global_compute_c;
template<typename T> faust_timer palm4MSA<T>::t_global_compute_lambda;
template<typename T> faust_timer palm4MSA<T>::t_global_update_R;
template<typename T> faust_timer palm4MSA<T>::t_global_update_L;
template<typename T> faust_timer palm4MSA<T>::t_global_check;
template<typename T> faust_timer palm4MSA<T>::t_global_init_fact;
template<typename T> faust_timer palm4MSA<T>::t_global_next_step;
template<typename T> faust_timer palm4MSA<T>::t_global_init_fact_from_palm;


template<typename T> faust_timer palm4MSA<T>::t_prox_const;
template<typename T> faust_timer palm4MSA<T>::t_prox_sp;
template<typename T> faust_timer palm4MSA<T>::t_prox_spcol;
template<typename T> faust_timer palm4MSA<T>::t_prox_splin;
template<typename T> faust_timer palm4MSA<T>::t_prox_normcol;

template<typename T> int palm4MSA<T>::nb_call_prox_const;
template<typename T> int palm4MSA<T>::nb_call_prox_sp;
template<typename T> int palm4MSA<T>::nb_call_prox_spcol;
template<typename T> int palm4MSA<T>::nb_call_prox_splin;
template<typename T> int palm4MSA<T>::nb_call_prox_normcol;




template<typename T>
void palm4MSA<T>::init_local_timers()
{
t_local_compute_projection.reset();
t_local_compute_grad_over_c.reset();
t_local_compute_c.reset();
t_local_compute_lambda.reset();
t_local_update_R.reset();
t_local_update_L.reset();
t_local_check.reset();
t_local_init_fact.reset();
t_local_next_step.reset();
t_local_init_fact_from_palm.reset();
}

template<typename T>
void palm4MSA<T>::print_global_timers()const
{
    cout << "timers in palm4MSA : " << endl;
    cout << "t_global_next_step           = " << t_global_next_step.get_time()           << " s for "<< t_global_next_step.get_nb_call()           << " calls" << endl;
    cout << "t grad + updateL + updateR   = " << t_global_compute_grad_over_c.get_time() + t_global_update_L.get_time() + t_global_update_R.get_time()  << " s for "<< t_global_compute_grad_over_c.get_nb_call()            << " calls of grad" << endl;
    cout << "t_global_compute_c = " << t_global_compute_c.get_time() << " s for "<< t_global_compute_c.get_nb_call() << " calls" << endl;
    cout << "t_global_compute_lambda      = " << t_global_compute_lambda.get_time()      << " s for "<< t_global_compute_lambda.get_nb_call()      << " calls" << endl;
    cout << "t_global_compute_projection  = " << t_global_compute_projection.get_time()  << " s for "<< t_global_compute_projection.get_nb_call()  << " calls" << endl<<endl;
    cout << "t_global_compute_grad_over_c = " << t_global_compute_grad_over_c.get_time() << " s for "<< t_global_compute_grad_over_c.get_nb_call() << " calls" << endl;
    cout << "t_global_update_R            = " << t_global_update_R.get_time()            << " s for "<< t_global_update_R.get_nb_call()            << " calls" << endl;
    cout << "t_global_update_L            = " << t_global_update_L.get_time()            << " s for "<< t_global_update_L.get_nb_call()            << " calls" << endl;
    cout << "t_check                      = " << t_global_check.get_time()               << " s for "<< t_global_check.get_nb_call()               << " calls" << endl;
    cout << "t_global_init_fact           = " << t_global_init_fact.get_time()           << " s for "<< t_global_init_fact.get_nb_call()           << " calls" << endl;
    cout << "t_global_init_fact_from_palm = " << t_global_init_fact_from_palm.get_time() << " s for "<< t_global_init_fact_from_palm.get_nb_call() << " calls" << endl<<endl;
}


template<typename T>
void palm4MSA<T>::print_prox_timers() const
{
/*	cout << "prox timers in palm4MSA : " << endl;
   cout << "total t_prox_const  =  " << t_prox_const.get_time()  << " s for "<< nb_call_prox_const  << " calls" << endl;
   cout << "total t_prox_sp  =  " << t_prox_sp.get_time()  << " s for "<< nb_call_prox_sp  << " calls" << endl;
   cout << "total t_prox_spcol  =  " << t_prox_spcol.get_time()  << " s for "<< nb_call_prox_spcol  << " calls" << endl;
   cout << "total t_prox_splin  =  " << t_prox_splin.get_time()  << " s for "<< nb_call_prox_splin  << " calls" << endl;
   cout << "total t_prox_normcol  =  " << t_prox_normcol.get_time()  << " s for "<< nb_call_prox_normcol  << " calls" << endl;
*/}


template<typename T>
void palm4MSA<T>::print_local_timers()const
{
    cout << "timers in palm4MSA : " << endl;
    cout << "t_local_next_step           = " << t_local_next_step.get_time()           << " s for "<< t_local_next_step.get_nb_call()           << " calls" << endl;
    cout << "t grad + updateL + updateR  = " << t_local_update_L.get_time()+t_local_update_R.get_time()+t_local_compute_grad_over_c.get_time()            << " s for "<< t_local_update_L.get_nb_call()            << " calls of grad" << endl;
    cout << "t local_compute_c  = " << t_local_compute_c.get_time()            << " s for "<< t_local_compute_c.get_nb_call()            << " calls of grad" << endl;
    cout << "t_local_compute_lambda      = " << t_local_compute_lambda.get_time()      << " s for "<< t_local_compute_lambda.get_nb_call()      << " calls" << endl;
    cout << "t_local_compute_projection  = " << t_local_compute_projection.get_time()  << " s for "<< t_local_compute_projection.get_nb_call()  << " calls" << endl<<endl;
    cout << "t_local_compute_grad_over_c = " << t_local_compute_grad_over_c.get_time() << " s for "<< t_local_compute_grad_over_c.get_nb_call() << " calls" << endl;
    cout << "t_local_update_R            = " << t_local_update_R.get_time()            << " s for "<< t_local_update_R.get_nb_call()            << " calls" << endl;
    cout << "t_local_update_L            = " << t_local_update_L.get_time()            << " s for "<< t_local_update_L.get_nb_call()            << " calls" << endl;
    cout << "t_check_                    = " << t_local_check.get_time()               << " s for "<< t_local_check.get_nb_call()               << " calls" << endl;
    cout << "t_local_init_fact           = " << t_local_init_fact.get_time()           << " s for "<< t_local_init_fact.get_nb_call()           << " calls" << endl;
    cout << "t_local_init_fact_from_palm = " << t_local_init_fact_from_palm.get_time() << " s for "<< t_local_init_fact_from_palm.get_nb_call() << " calls" << endl<<endl;
}


#endif


#endif
