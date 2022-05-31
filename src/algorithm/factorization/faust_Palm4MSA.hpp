/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2021):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
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
#ifndef __PALM4MSA_HPP__
#define __PALM4MSA_HPP__

//#include "faust_Palm4MSA.h"

#include <sstream>

#include <iostream>
#include "faust_Params.h"
#include "faust_ParamsPalm.h"
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintInt.h"
#include "faust_ConstraintMat.h"
#include "faust_MatDense.h"

#include "faust_linear_algebra.h"
#include "faust_prod_opt.h"
#include "faust_prox.h"
#include "faust_ConstraintType.h"
#include "faust_exception.h"

#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif

//#include"faust_init_from_matio_mat.h"

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#define __SP setprecision(20)<<

using namespace std;

template<typename FPP,FDevice DEVICE,typename FPP2>
const char * Faust::Palm4MSA<FPP,DEVICE,FPP2>::m_className="Faust::Palm4MSA";

template<typename FPP,FDevice DEVICE,typename FPP2>
const FPP2 Faust::Palm4MSA<FPP,DEVICE,FPP2>::lipschitz_multiplicator=1.001;


template<typename FPP,FDevice DEVICE,typename FPP2>
Faust::Palm4MSA<FPP,DEVICE,FPP2>::Palm4MSA(const Faust::MatDense<FPP,DEVICE>& M, const Faust::Params<FPP,DEVICE,FPP2> & params_, const bool isGlobal_) :
    data(M),
    m_lambda(params_.init_lambda),
    m_nbFact(0),
    S(params_.init_fact),
    RorL(vector<Faust::MatDense<FPP,DEVICE> >(2)),
    m_indFact(0),
    m_indIte(-1),
    verbose(params_.isVerbose),
    isUpdateWayR2L(params_.isUpdateWayR2L),
    isConstantStepSize(params_.isConstantStepSize),
	gradCalcOptMode(params_.gradCalcOptMode),
    isGradComputed(false),
    isProjectionComputed(false),
    isLastFact(false),
    isConstraintSet(false),
    isGlobal(isGlobal_),
    isInit(false),
    c(FPP2(1)/params_.step_size),
	norm2_threshold(params_.norm2_threshold),
    norm2_max_iter(params_.norm2_max_iter),
	is_complex(typeid(data.getData()[0]) == typeid(complex<float>) || typeid(data.getData()[0]) == typeid(complex<double>)),
	TorH(is_complex?'H':'T')
{
   RorL.reserve(params_.m_nbFact);

   if(isGlobal)
      stop_crit = Faust::StoppingCriterion<FPP2>(params_.stop_crit_global);
   else
      stop_crit = Faust::StoppingCriterion<FPP2>(params_.stop_crit_2facts);

	if (isConstantStepSize)
        isCComputed = true;
	else
        isCComputed = false;

   if(verbose)
   {
	   spectral_duration = std::chrono::duration<double>::zero();
	   fgrad_duration = std::chrono::duration<double>::zero();
   }
}

template<typename FPP,FDevice DEVICE,typename FPP2>
Faust::Palm4MSA<FPP,DEVICE,FPP2>::Palm4MSA(const Faust::ParamsPalm<FPP,DEVICE,FPP2>& params_palm_, const bool isGlobal_/*=false*/) :
	stop_crit(params_palm_.stop_crit),
	data(params_palm_.data),
	m_lambda(params_palm_.init_lambda),
	m_nbFact(params_palm_.nbFact),
	S(params_palm_.init_fact),
	RorL(vector<Faust::MatDense<FPP,DEVICE> >(2)),
	LorR(Faust::MatDense<FPP,DEVICE>(params_palm_.init_fact[0].getNbRow())),
	const_vec(params_palm_.cons),
	m_indFact(0),
	m_indIte(-1),
	verbose(params_palm_.isVerbose),
	isUpdateWayR2L(params_palm_.isUpdateWayR2L),
	isConstantStepSize(params_palm_.isConstantStepSize),
	gradCalcOptMode(params_palm_.gradCalcOptMode),
	isGradComputed(false),
	isProjectionComputed(false),
	isLastFact(false),
	isConstraintSet(false),
	isGlobal(isGlobal_),
	c(FPP2(1)/params_palm_.step_size),
    norm2_threshold(params_palm_.norm2_threshold),
    norm2_max_iter(params_palm_.norm2_max_iter),
	is_complex(typeid(data.getData()[0]) == typeid(complex<float>) || typeid(data.getData()[0]) == typeid(complex<double>)
),
	TorH(is_complex?'H':'T')
{
	RorL.reserve(const_vec.size()+1);

   	if (isConstantStepSize)
		isCComputed = true;
	else
		isCComputed = false;

   check_constraint_validity();
   if(verbose)
   {
	   spectral_duration = std::chrono::duration<double>::zero();
	   fgrad_duration = std::chrono::duration<double>::zero();
   }
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_facts()
{
	while (do_continue())
	{
		next_step();
	}
	if(verbose)
	{
		std::cout << "palm4msa spectral time=" << spectral_duration.count() << std::endl;
		std::cout << "palm4msa fgrad time=" << fgrad_duration.count() << std::endl;
		spectral_duration = std::chrono::duration<double>::zero();
		fgrad_duration = std::chrono::duration<double>::zero();
	}
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::get_facts(Faust::Transform<FPP,DEVICE> & faust_fact) const
{
	Faust::Transform<FPP,DEVICE> f(S);
	faust_fact = f;


}


template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_projection()
{
#ifdef __COMPILE_TIMERS__
t_global_compute_projection.start();
t_local_compute_projection.start();
#endif


      //Faust::MatDense<FPP,DEVICE> matrix2project(S[m_indFact]);
      S[m_indFact] -= grad_over_c;
      const_vec[m_indFact]->template project<FPP,DEVICE,FPP2>(S[m_indFact]);


    isProjectionComputed = true;
#ifdef __COMPILE_TIMERS__
t_global_compute_projection.stop();
t_local_compute_projection.stop();
#endif
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_grad_over_c_ext_opt()
{
//#define mul_3_facts multiply_order_opt
#define mul_3_facts multiply_order_opt_all_ends// this one only optimizes the product on factor ends but for three factors it doesn't change anything comparing to multiply_order_opt
	// compute error = m_lambda*L*S*R-data
	error = data;
	std::vector<MatDense<FPP,DEVICE>*> facts;
	std::vector<char> tc_flags;
	if(isUpdateWayR2L)
		facts = { &RorL[m_indFact], &S[m_indFact], &LorR };
	else
		facts = { &LorR, &S[m_indFact], &RorL[m_indFact] };
	mul_3_facts(facts, error, (FPP) m_lambda, (FPP) -1.0);
	// compute m_lambda/c * L'*error*R'
	if(isUpdateWayR2L)
		facts = { &RorL[m_indFact], &error, &LorR };
	else
		facts = {&LorR, &error, &RorL[m_indFact]};
	tc_flags = {TorH, 'N', TorH};
	mul_3_facts(facts, grad_over_c, FPP(m_lambda/c), (FPP)0, tc_flags);
	isGradComputed = true;
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_grad_over_c()
{
	if(verbose)
		fgrad_start = std::chrono::high_resolution_clock::now();
	if(gradCalcOptMode == EXTERNAL_OPT)
		compute_grad_over_c_ext_opt();
	else
		compute_grad_over_c_int_opt();
	if(verbose)
	{
		fgrad_stop = std::chrono::high_resolution_clock::now();
		fgrad_duration += fgrad_stop-fgrad_start;
	}
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_grad_over_c_int_opt()
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
        R2 = (unsigned long long int) RorL[m_indFact].getNbCol();
    }
    else
    {
        L1 = (unsigned long long int) RorL[m_indFact].getNbRow();
        L2 = (unsigned long long int) RorL[m_indFact].getNbCol();
        R2 = (unsigned long long int) LorR.getNbCol();
    }
    S2 = (unsigned long long int) S[m_indFact].getNbCol();
    vector<unsigned long long int > complexity(4,0);
    complexity[0] = L1*L2*S2 + L1*S2*R2 + L2*L1*R2 + L2*R2*R2;
    complexity[1] = L1*L2*S2 + L1*S2*R2 + L1*R2*S2 + L2*L1*S2;
    complexity[2] = L2*S2*R2 + L1*L2*R2 + L2*L1*R2 + L2*R2*S2;
    complexity[3] = L2*S2*R2 + L1*L2*R2 + L1*R2*S2 + L2*L1*S2;

    int idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));


    error = data;
    Faust::MatDense<FPP,DEVICE> tmp1,tmp3;

   if (idx==0 || idx==1 || gradCalcOptMode == DISABLED) // computing L*S first, then (L*S)*R
   {
      if (!isUpdateWayR2L)
      {
         // tmp1 = L*S
         multiply(LorR, S[m_indFact], tmp1);
/*sprintf(nomFichier,"LorR_0_%d_device.tmp",cmpt);
LorR.print_file(nomFichier);
sprintf(nomFichier,"RorL%d_0_%d_host.tmp",m_indFact,cmpt);
RorL[m_indFact].print_file(nomFichier);
sprintf(nomFichier,"S%d_0_%d_host.tmp",m_indFact,cmpt);
S[m_indFact].print_file(nomFichier);
sprintf(nomFichier,"tmp1_0_%d_host.tmp",cmpt);
tmp1.print_file(nomFichier);
sprintf(nomFichier,"error_0_%d_host.tmp",cmpt);
error.print_file(nomFichier);
cout << "appel " << cmpt<<" : lambda0 = "<< m_lambda<<endl;*/
         // error = m_lambda*tmp1*R - error (= m_lambda*L*S*R - data )
         gemm(tmp1, RorL[m_indFact], error, FPP(m_lambda), (FPP)-1.0, 'N', 'N');
/*sprintf(nomFichier,"LorR_1_%d_device.tmp",cmpt);
LorR.print_file(nomFichier);
sprintf(nomFichier,"S%d_1_%d_host.tmp",m_indFact,cmpt);
S[m_indFact].print_file(nomFichier);
sprintf(nomFichier,"tmp1_1_%d_host.tmp",cmpt);
tmp1.print_file(nomFichier);
sprintf(nomFichier,"RorL%d_1_%d_host.tmp",m_indFact,cmpt);
RorL[m_indFact].print_file(nomFichier);
sprintf(nomFichier,"error_1_%d_device.tmp",cmpt);*/
      }
      else
      {
         // tmp1 = L*S
         multiply(RorL[m_indFact], S[m_indFact], tmp1);
         // error = m_lambda*tmp1*R - error (= m_lambda*L*S*R - data )
         gemm(tmp1, LorR, error, FPP(m_lambda),(FPP) -1.0, 'N', 'N');
      }
   }
   else // computing S*R first, then L*(S*R)
   {
      if (!isUpdateWayR2L)
      {
         // tmp1 = S*R
         multiply(S[m_indFact], RorL[m_indFact], tmp1);
/*sprintf(nomFichier,"LorR_0_%d_device.tmp",cmpt);
LorR.print_file(nomFichier);
sprintf(nomFichier,"RorL%d_0_%d_device.tmp",m_indFact,cmpt);
RorL[m_indFact].print_file(nomFichier);
sprintf(nomFichier,"S%d_0_%d_device.tmp",m_indFact,cmpt);
S[m_indFact].print_file(nomFichier);
sprintf(nomFichier,"tmp1_0_%d_device.tmp",cmpt);
tmp1.print_file(nomFichier);
sprintf(nomFichier,"error_0_%d_device.tmp",cmpt);
error.print_file(nomFichier);
cout << "appel " << cmpt<<" : lambda0 = "<< m_lambda<<endl;*/

         // error = m_lambda*L*tmp1 - error (= m_lambda*L*S*R - data )
         gemm(LorR, tmp1, error, FPP(m_lambda),(FPP) -1.0, 'N', 'N');
/*sprintf(nomFichier,"LorR_1_%d_device.tmp",cmpt);
LorR.print_file(nomFichier);
sprintf(nomFichier,"S%d_1_%d_device.tmp",m_indFact,cmpt);
S[m_indFact].print_file(nomFichier);
sprintf(nomFichier,"tmp1_1_%d_device.tmp",cmpt);
tmp1.print_file(nomFichier);
sprintf(nomFichier,"RorL%d_1_%d_device.tmp",m_indFact,cmpt);
RorL[m_indFact].print_file(nomFichier);
sprintf(nomFichier,"error_1_%d_device.tmp",cmpt);*/

      }
      else
      {
         // tmp1 = S*R
         multiply(S[m_indFact], LorR, tmp1);
         // error = m_lambda*L*tmp1 - error (= m_lambda*L*S*R - data )
         gemm(RorL[m_indFact], tmp1, error, FPP(m_lambda),(FPP) -1.0, 'N', 'N');
      }
   }

   if (idx==0 || idx==2 || gradCalcOptMode == DISABLED) // computing L'*error first, then (L'*error)*R'
   {
      if (!isUpdateWayR2L)
      {
         // tmp3 = m_lambda*L'*error (= m_lambda*L' * (m_lambda*L*S*R - data) )
         gemm(LorR, error, tmp3, FPP(m_lambda),(FPP) 0.0, TorH, 'N');
         // grad_over_c = 1/c*tmp3*R' (= 1/c*m_lambda*L' * (m_lambda*L*S*R - data) * R' )
         gemm(tmp3, RorL[m_indFact], grad_over_c, FPP(1.0/c),(FPP) 0.0,'N',TorH);
      }
      else
      {
         // tmp3 = m_lambda*L'*error (= m_lambda*L' * (m_lambda*L*S*R - data) )
         gemm(RorL[m_indFact], error, tmp3, FPP(m_lambda), (FPP) 0.0, TorH, 'N');
         // grad_over_c = 1/c*tmp3*R' (= 1/c*m_lambda*L' * (m_lambda*L*S*R - data) * R' )
         gemm(tmp3, LorR, grad_over_c, FPP(1.0/c), (FPP) (FPP) 0.0,'N',TorH);
      }
   }
   else // computing error*R' first, then L'*(error*R')
   {
      if (!isUpdateWayR2L)
      {
         // tmp3 = m_lambda*error*R' (= m_lambda*(m_lambda*L*S*R - data) * R' )
         gemm(error, RorL[m_indFact], tmp3, FPP(m_lambda), (FPP) 0.0, 'N', TorH);
         // grad_over_c = 1/c*L'*tmp3 (= 1/c*L' * m_lambda*(m_lambda*L*S*R - data) * R' )
         gemm(LorR, tmp3, grad_over_c, FPP(1.0/c), (FPP) 0.0,TorH,'N');
      }
      else
      {
         // tmp3 = m_lambda*error*R' (= m_lambda * (m_lambda*L*S*R - data) * R' )
         gemm(error, LorR, tmp3, FPP(m_lambda), (FPP) 0.0, 'N', TorH);
         // grad_over_c = 1/c*L'*tmp3 (= 1/c*L' * m_lambda*(m_lambda*L*S*R - data) * R' )
         gemm(RorL[m_indFact], tmp3, grad_over_c, FPP(1.0/c), (FPP) 0.0,TorH,'N');
      }

    }

    isGradComputed = true;

#ifdef __COMPILE_TIMERS__
t_global_compute_grad_over_c.stop();
t_local_compute_grad_over_c.stop();
#endif
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_lambda()
{
	compute_lambda(LorR);
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_lambda(Faust::MatDense<FPP,DEVICE>& LorR)
{
#ifdef __COMPILE_TIMERS__
t_global_compute_lambda.start();
t_local_compute_lambda.start();
#endif

   if (!isLastFact)
   {
     handleError(m_className,"compute_lambda : computation of lambda must be done at the end of the iteration through the number of factors");
   }

   // As LorR has also been updated at the end of the last iteration over the facts, LorR matches X_hat, which the product of all factors, including the last one.
   // Xt_Xhat = data'*X_hat
   Faust::MatDense<FPP,DEVICE> Xt_Xhat;
   gemm(data, LorR, Xt_Xhat, (FPP) 1.0, (FPP) 0.0, TorH,'N');
   // Xhatt_Xhat = X_hat'*X_hat
   Faust::MatDense<FPP,DEVICE> Xhatt_Xhat;
   gemm(LorR, LorR, Xhatt_Xhat, (FPP) 1.0, (FPP) 0.0, TorH,'N');

   FPP Xhatt_Xhat_tr = (FPP) Xhatt_Xhat.trace();

   if (Xhatt_Xhat_tr != FPP(0))
   {
		m_lambda = Faust::fabs(Xt_Xhat.trace()/Xhatt_Xhat_tr);
		if (std::isnan(m_lambda))
		{
			handleError(m_className,"compute_lambda : Xhatt_Xhat_tr is too small or Xt_Xhat.trace is too big so lambda is infinite");
		}
	}else
	{
		handleError(m_className,"compute_lambda : Xhatt_Xhat_tr equal 0 so lambda is infinite");
	}
   //cout<<"m_lambda : "<<m_lambda<<endl;
   //cout<<__SP m_lambda<<endl;

#ifdef __COMPILE_TIMERS__
t_global_compute_lambda.stop();
t_local_compute_lambda.stop();
#endif
}


template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::update_R()
{
#ifdef __COMPILE_TIMERS__
t_global_update_R.start();
t_local_update_R.start();
#endif

   // R[nbFact-1] est initialise a l'identite lors de la creation de l'objet Faust::Palm4MSA et n'est pas cense changer
   if (!isUpdateWayR2L)
   {
      RorL.resize(m_nbFact);
      RorL[m_nbFact-1].resize(const_vec[m_nbFact-1]->get_cols());
      RorL[m_nbFact-1].setEyes();


      for (int i=m_nbFact-2 ; i>-1 ; i--)
         //  R[i] = S[i+1] * R[i+1]
         multiply(S[i+1], RorL[i+1], RorL[i]);

   }
   else
   {
      if(!isProjectionComputed)
      {
         throw std::logic_error("Projection must be computed before updating L");
      }
	   // FORMER VERSION //
	   //Faust::MatDense<FPP,DEVICE> LorR_tmp(LorR);
	   //multiply(S[m_indFact],LorR_tmp,LorR);
	   //END FORMER VERSION//
	   multiply(S[m_indFact], LorR, LorR);

      //LorR.multiplyLeft(S[m_indFact]);
   }



#ifdef __COMPILE_TIMERS__
t_global_update_R.stop();
t_local_update_R.stop();
#endif
}


template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::update_L()
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


       // FORMER VERSION //
	   //Faust::MatDense<FPP,DEVICE> LorR_tmp(LorR);
	   //multiply(LorR_tmp,S[m_indFact],LorR);
   	   // END FORMER VERSION //
	   multiply(LorR, S[m_indFact], LorR);
	}
   else
   {
      RorL.resize(m_nbFact);
      RorL[0].resize(const_vec[0]->get_rows());
      RorL[0].setEyes();
      for (int i=0 ; i<m_nbFact-1 ; i++)
         //  R[i] = S[i+1] * R[i+1]
         multiply(RorL[i] , S[i], RorL[i+1]);
   }
#ifdef __COMPILE_TIMERS__
t_global_update_L.stop();
t_local_update_L.stop();
#endif
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::check_constraint_validity()
{
#ifdef __COMPILE_TIMERS__
t_global_check.start();
t_local_check.start();
#endif

   if (m_nbFact != S.size())
   {
      handleError(m_className," Wrong initialization: params.nfacts and params.init_facts are in conflict");
   }
#ifdef __COMPILE_TIMERS__
t_global_check.stop();
t_local_check.stop();
#endif
}


#ifdef __PAS_FIXE__
// Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_c() has been defined as an inline method in Faust::Palm4MSA.h
#else
template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::compute_c()
{
#ifdef __COMPILE_TIMERS__
	t_global_compute_c.start();
	t_local_compute_c.start();
#endif

   if (!isConstantStepSize)
   {
	   int flag1,flag2;
	   if(verbose)
		   spectral_start = std::chrono::high_resolution_clock::now();
	   FPP2 nL1=LorR.spectralNorm(norm2_max_iter,norm2_threshold,flag1);
	   FPP2 nR1=RorL[m_indFact].spectralNorm(norm2_max_iter,norm2_threshold,flag2);
	   if(verbose)
		   spectral_stop = std::chrono::high_resolution_clock::now();
	   c=lipschitz_multiplicator*nR1*nR1*nL1*nL1*m_lambda*m_lambda;
	   if(verbose)
		   spectral_duration += spectral_stop-spectral_start;
   }

   isCComputed = true;


   #ifdef __COMPILE_TIMERS__
	t_global_compute_c.stop();
	t_local_compute_c.stop();
	#endif
}
#endif

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::init_fact(int nbFacts_)
{
#ifdef __COMPILE_TIMERS__
t_global_init_fact.start();
t_local_init_fact.start();
#endif

  /*if(isInit && isGlobal)
  {
     cerr << "Error in Faust::Palm4MSA<FPP,DEVICE,FPP2>::init_fact : global factorization has already been initialized" <<endl;
     exit(EXIT_FAILURE);
  }
  else if (isGlobal)
  {
     m_nbFact = 1;
     S.resize(m_nbFact);
     isInit = true;
     return;
  }*/



// if !isGlobal

  if(!isConstraintSet)
  {
     handleError(m_className,"init_fact : constrainst must be set before calling init_fact");
  }
   m_nbFact = nbFacts_;
   S.resize(m_nbFact);
   if (!isUpdateWayR2L)
   {
      S[0].resize(const_vec[0]->get_rows(), const_vec[0]->get_cols());
      S[0].setZeros();
      for (int i=1 ; i<m_nbFact ; i++)
      {
         S[i].resize(const_vec[i]->get_rows(), const_vec[i]->get_cols());
         S[i].setEyes();
      }
   }
   else
   {
      for (int i=0 ; i<m_nbFact-1 ; i++)
      {
         S[i].resize(const_vec[i]->get_rows(), const_vec[i]->get_cols());
         S[i].setEyes();
      }
      S[m_nbFact-1].resize(const_vec[m_nbFact-1]->get_rows(), const_vec[m_nbFact-1]->get_cols());
      S[m_nbFact-1].setZeros();
   }


#ifdef __COMPILE_TIMERS__
t_global_init_fact.stop();
t_local_init_fact.stop();
#endif
}

template<typename FPP,FDevice DEVICE,typename FPP2>
bool Faust::Palm4MSA<FPP,DEVICE,FPP2>::do_continue()
{
	bool cont;
	FPP2 err = -1;
	if(m_indIte > 1 && stop_crit.isCriterionErr())
	//compute error (only if at least one iteration has been executed
		err = error.norm();
	cont = stop_crit.do_continue(++m_indIte, err);
	if(!cont)
	{
		m_indIte=-1;
		isConstraintSet=false;
	}
	return cont;
}

/** \brief
 *
 * \param
 * \param
 */
template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::next_step()
{
#ifdef __COMPILE_TIMERS__
t_global_next_step.start();
t_local_next_step.start();
#endif

    check_constraint_validity();
    // resizing L or R
    if(!isUpdateWayR2L)
    {
        LorR.resize(const_vec[0]->get_rows());
        LorR.setEyes();
        update_R();
    }
    else
    {
        LorR.resize(const_vec[m_nbFact-1]->get_cols());
        LorR.setEyes();
        update_L();
    }

    int* ind_ptr = new int[m_nbFact];
    for (int j=0 ; j<m_nbFact ; j++)
		if (!isUpdateWayR2L)
			ind_ptr[j] = j;
		else
			ind_ptr[j] = m_nbFact-1-j;

    for (int j=0 ; j<m_nbFact ; j++)
    {
        if (j == m_nbFact-1)
            isLastFact = true;
        else
            isLastFact = false;

        m_indFact = ind_ptr[j];
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
        cout << "Iter " << m_indIte << ", RMSE=" << get_RMSE() << ", RE=" << get_RE()<< endl;
        cout << "m_lambda " <<setprecision(20)<< m_lambda << endl;
    }
    delete[] ind_ptr;
    ind_ptr = NULL;


//cout<<"m_lambda : "<< m_lambda<< endl;
#ifdef __COMPILE_TIMERS__
t_global_next_step.stop();
t_local_next_step.stop();
#endif
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::init_fact_from_palm(const Faust::Palm4MSA<FPP,DEVICE,FPP2>& palm2, bool isFactSideLeft)
{
#ifdef __COMPILE_TIMERS__
t_global_init_fact_from_palm.start();
t_local_init_fact_from_palm.start();
#endif


    if (palm2.m_nbFact != 2)
    {
        handleError(m_className,"init_fact_from_palm : argument palm2 must contain 2 factors.");
    }

    if(!isConstraintSet)
    {
        handleError(m_className,"init_fact_from_palm : constraints must be set before calling init_fact_from_palm");
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
    m_nbFact++;

    check_constraint_validity();

#ifdef __COMPILE_TIMERS__
t_global_init_fact_from_palm.stop();
t_local_init_fact_from_palm.stop();
#endif
}


#ifdef __COMPILE_TIMERS__

template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_compute_projection;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_compute_grad_over_c;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_compute_c;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_compute_lambda;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_update_R;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_update_L;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_check;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_init_fact;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_next_step;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_global_init_fact_from_palm;


template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_prox_const;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_prox_sp;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_prox_spcol;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_prox_splin;
template<typename FPP,FDevice DEVICE,typename FPP2> Faust::Timer Faust::Palm4MSA<FPP,DEVICE,FPP2>::t_prox_normcol;

template<typename FPP,FDevice DEVICE,typename FPP2> int Faust::Palm4MSA<FPP,DEVICE,FPP2>::nb_call_prox_const;
template<typename FPP,FDevice DEVICE,typename FPP2> int Faust::Palm4MSA<FPP,DEVICE,FPP2>::nb_call_prox_sp;
template<typename FPP,FDevice DEVICE,typename FPP2> int Faust::Palm4MSA<FPP,DEVICE,FPP2>::nb_call_prox_spcol;
template<typename FPP,FDevice DEVICE,typename FPP2> int Faust::Palm4MSA<FPP,DEVICE,FPP2>::nb_call_prox_splin;
template<typename FPP,FDevice DEVICE,typename FPP2> int Faust::Palm4MSA<FPP,DEVICE,FPP2>::nb_call_prox_normcol;




template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::init_local_timers()
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

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::print_global_timers()const
{
    cout << "timers in Faust::Palm4MSA : " << endl;
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


template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::print_prox_timers() const
{
/*	cout << "prox timers in Faust::Palm4MSA : " << endl;
   cout << "total t_prox_const  =  " << t_prox_const.get_time()  << " s for "<< nb_call_prox_const  << " calls" << endl;
   cout << "total t_prox_sp  =  " << t_prox_sp.get_time()  << " s for "<< nb_call_prox_sp  << " calls" << endl;
   cout << "total t_prox_spcol  =  " << t_prox_spcol.get_time()  << " s for "<< nb_call_prox_spcol  << " calls" << endl;
   cout << "total t_prox_splin  =  " << t_prox_splin.get_time()  << " s for "<< nb_call_prox_splin  << " calls" << endl;
   cout << "total t_prox_normcol  =  " << t_prox_normcol.get_time()  << " s for "<< nb_call_prox_normcol  << " calls" << endl;
*/}


template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::print_local_timers()const
{
    cout << "timers in Faust::Palm4MSA : " << endl;
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

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::Palm4MSA<FPP,DEVICE,FPP2>::print_state() const
{
	std::cout << std::string(7, '=') << " PALM4MSA state: " << "(iter: "  m_indIte <<  ") ===" << std::endl;
	std::cout << "m_indFact: " << m_indFact << std::endl;
	std::cout << "isUpdateWayR2L: " << isUpdateWayR2L << std::endl;
	std::cout << "m_lambda: " << m_lambda << std::endl;
	std::cout << "S: ";
	for(auto f: S)
		std::cout << f.norm() << " ";
	cout << std::endl;
	std::cout << "RorL: ";
	for(auto f: RorL)
		std::cout << f.norm() << " ";
	cout << std::endl;
	std::cout << "LorR: ";
	std::cout << LorR.norm() << " ";
	std::cout << LorR.to_string(false, true);
	cout << std::endl;
	std::cout << "grad_over_c: " << grad_over_c.norm() << std::endl;
	std::cout << "c: " << c << std::endl;
	std::cout << std::string(25, '=') << std::endl;
}

#endif


#endif
