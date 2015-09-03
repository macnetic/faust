#include "sparsePalm4MSA.h"

#include <iostream>
#include "faust_params.h"
#include "faust_params_palm.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_mat.h"
#include "faust_constraint_real.h"
#include "faust_constraint_int.h"
//#include "faust_mat.h"
#include "LinAlgebra.h"
#include "prox.h"

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

sparsePalm4MSA::sparsePalm4MSA(const faust_params& params_, const bool isGlobal_) :
   data(params_.data),
   lambda(params_.init_lambda),
   nb_fact(0),
   ind_ite(-1),
   lipschitz_multiplicator(1.001),
   verbose(params_.isVerbose),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   //isLambdaComputed(params_.isLambdaComputed),
   isCComputed(false),
   isGradComputed(false),
   isProjectionComputed(false),
   isLastFact(false),
   isConstraintSet(false)

{
   if(isGlobal_)
      stop_crit = stopping_criterion(params_.stop_crit_global);
   else
      stop_crit = stopping_criterion(params_.stop_crit_2facts);
}

sparsePalm4MSA::sparsePalm4MSA(const faust_params_palm& params_palm_) :
   data(params_palm_.data),
   lambda(params_palm_.init_lambda),
   nb_fact(params_palm_.nb_fact),
   sparse_S(),
   dense_S(),
   stop_crit(params_palm_.stop_crit),
   const_vec(params_palm_.cons),
   ind_ite(-1),
   lipschitz_multiplicator(1.001),
   verbose(params_palm_.isVerbose),
   isUpdateWayR2L(params_palm_.isUpdateWayR2L),
   //isLambdaComputed(params_palm_.isLambdaComputed),
   isCComputed(false),
   isGradComputed(false),
   isProjectionComputed(false),
   isLastFact(false),
   isConstraintSet(false)

{
   	if (isUpdateWayR2L)
	{	
		ind_fact = params_palm_.init_fact.size()-1;
		for (int i=0;i<ind_fact;i++)
		{
			dense_S = params_palm_.init_fact[i];
			sparse_S = dense_S;
			L.push_back(sparse_S);
		}
		
		dense_S = params_palm_.init_fact[ind_fact];
		sparse_S = dense_S;
	}
	else
	{	
		ind_fact = 0;
		for (int i=1;i<params_palm_.init_fact.size();i++)
		{
			dense_S = params_palm_.init_fact[i];
			sparse_S = dense_S;
			R.push_back(sparse_S);
		}
		dense_S = params_palm_.init_fact[0];
		sparse_S = dense_S;
	}
	
	
   check_constraint_validity();

}

faust_spmat sparsePalm4MSA::get_residuum(const bool isFactSideLeft) const
{
	if (isFactSideLeft)
	{
		return get_fact(0);
	}else
	{
		return get_fact(nb_fact-1);
	}

}
void sparsePalm4MSA::compute_facts()
{
	while(do_continue())
	{
		next_step();
	}
}

void sparsePalm4MSA::Display() const
{
	std::cout<<"data :  nbRow "<<data.getNbRow()<<" NbCol : "<< data.getNbCol()<<std::endl;
	std::cout<<"nb_facts : "<<nb_fact<<std::endl;
	std::cout<<"Constraints :"<<std::endl;
	std::cout<<"Constraints size:"<<const_vec.size()<<std::endl;
    std::cout<<"UpdatewayR2L:"<<isUpdateWayR2L<<std::endl;
	 std::cout<<"size L and R: "<<L.size()<<" "<<R.size()<<" : sum LSR size "<<L.size()+R.size()+1<<std::endl;
	
	for (int i=0;i<const_vec.size();i++)
		{	
			
			std::cout<<"type : "<<(const_vec[i]->get_constraint_name());
			//const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(current_line[i]);
			//std::cout<<" parameter : "<<const_int->getParameter();
			std::cout<<" DIM1 : "<<const_vec[i]->getRows()<<" DIM2 : "<<const_vec[i]->getCols()<<std::endl;
		}
	std::cout<<std::endl<<std::endl;
}


#ifdef __PAS_FIXE__
// palm4MSA::compute_c() has been defined as an inline method in palm4MSA.h
#else
void sparsePalm4MSA::compute_c()
{
#ifdef __COMPILE_TIMERS__
	t_global_compute_c.start();
	t_local_compute_c.start();
#endif

   //std::cout<<"calcul pas : "<<std::endl;
   // faust_real nL=LorR.spectralNorm();
   //faust_real nR=RorL[ind_fact].spectralNorm();
   //c=lipschitz_multiplicator*nR*nR*nL*nL*lambda*lambda;


   int flag1,flag2;
   
   int nbr_iter = 10000;
   faust_real threshold = 1e-16;
   faust_real nL1=L.spectralNorm(nbr_iter,threshold,flag1);
   faust_real nR1=R.spectralNorm(nbr_iter,threshold,flag2);
    c=lipschitz_multiplicator*nR1*nR1*nL1*nL1*lambda*lambda;

   //std::cout<<" nL : "<<nL <<" nL1 : "<<nL1<<" flag : "<<flag1<<std::endl;
   //std::cout<<" nR : "<<nR <<" nR1 : "<<nR1<<" flag : "<<flag2<<std::endl;
   //std::cout<<" c : "<< c <<" c1 : "<<c1<<std::endl<<std::endl;
   //std::cout<<" c : "<< c <<std::endl;


   isCComputed = true;
   
   #ifdef __COMPILE_TIMERS__
	t_global_compute_c.stop();
	t_local_compute_c.stop();
	#endif	
}
#endif



























void sparsePalm4MSA::next_step()
{	
	#ifdef __COMPILE_TIMERS__
		t_global_next_step.start();
		t_local_next_step.start();
	#endif
	//cout<<"debut sppalm4MSA next_step"<<endl;
	for (int i=0;i<nb_fact;i++)
	{	
		//cout<<"iter i:"<<i<<endl;
		isCComputed = false;
		isGradComputed = false;
		isProjectionComputed = false;
		
		
		compute_c();
		//cout<<"debut compute_grad_over_c"<<endl;
		compute_grad_over_c();
		//cout<<"fin compute_grad_over_c"<<endl;
		compute_projection();
		//cout<<"fin compute_proj"<<endl;
		update_LSR();
		//cout<<"fin update_LSR"<<endl;
		
	}
		compute_lambda();
		last_update();
	#ifdef __COMPILE_TIMERS__
		t_global_next_step.stop();
		t_local_next_step.stop();
	#endif
}

faust_spmat sparsePalm4MSA::get_fact(const int id_fact) const
{	
	int current_id = id_fact;
	if ( (id_fact>=nb_fact) ||  (id_fact < 0) )
	{
		 cerr << "ERROR sparsePalm4MSA::get_fact : id_fact exceeds number of facts" << endl;
      exit(EXIT_FAILURE);
	}
	
	if (current_id<L.size())
	{
		return L.get_fact(current_id);
	}else
	{
		current_id-=L.size();
	}
	
	if (current_id == 0)
	{
		return sparse_S;
	}else
	{
		return R.get_fact(current_id-1);
	}
}


void sparsePalm4MSA::compute_grad_over_c()
{
	#ifdef __COMPILE_TIMERS__
		t_global_compute_grad_over_c.start();
		t_local_compute_grad_over_c.start();
	#endif
	if(!isCComputed) 
   {
      cerr << "c must be set before computing grad/c" << endl;
      exit(EXIT_FAILURE);
   }
	//cout<<"L_size : "<<L.size()<<" : dim "<<L.getNbRow()<<" "<<L.getNbCol()<<endl;
	//cout<<"sparse_S size : "<<sparse_S.getNbRow()<<" "<<sparse_S.getNbCol()<<endl;
	//cout<<"R_size : "<<R.size()<<" : dim "<<R.getNbRow()<<" "<<R.getNbCol()<<endl;
   	faust_core LSR(L);
	faust_mat  LSR_prod,tmp1;
	LSR *= sparse_S;
	LSR *= R;
	LSR *= lambda;
	LSR_prod=LSR.get_product();
	LSR_prod -= data;
	//cout<<"before multiply(L,LSR_prod,tmp1,lambda,'T','R');"<<endl;
	multiply(L,LSR_prod,tmp1,lambda,'T','R');
	//cout<<"after multiply(L,LSR_prod,tmp1,lambda,'T','R');"<<endl;
	multiply(R,tmp1,grad_over_c,1.0/c,'T','L');
	//cout<<"sparse_S size : "<<sparse_S.getNbRow()<<" "<<sparse_S.getNbCol()<<endl;
	//cout<<"dense_S size : "<<dense_S.getNbRow()<<" "<<dense_S.getNbCol()<<endl;
	dense_S = sparse_S;
	
	
	isGradComputed = true;
	#ifdef __COMPILE_TIMERS__
		t_global_compute_grad_over_c.stop();
		t_local_compute_grad_over_c.stop();
	#endif
}




  
void sparsePalm4MSA::compute_projection()
{
#ifdef __COMPILE_TIMERS__
t_global_compute_projection.start();
t_local_compute_projection.start();
#endif

   if (const_vec[ind_fact]->getConstraintType() == CONSTRAINT_NAME_CONST)
   {	
		#ifdef __COMPILE_TIMERS__
			nb_call_prox_const++;
			t_prox_const.start();
			
		#endif
         const faust_constraint_mat* const_mat = dynamic_cast<const faust_constraint_mat*>(const_vec[ind_fact]);
         //S[ind_fact] = const_mat->getParameter();
		 #ifdef __COMPILE_TIMERS__
			t_prox_const.stop();
		#endif
   }
   else
   {
      //faust_mat matrix2project(S[ind_fact]);
	  //cout<<"dense_S size : "<<dense_S.getNbRow()<<" "<<dense_S.getNbCol()<<endl;
	  //cout<<"grad size : "<<grad_over_c.getNbRow()<<" "<<grad_over_c.getNbCol()<<endl;
      dense_S -= grad_over_c;
      switch (const_vec[ind_fact]->getConstraintType())
      {
         case CONSTRAINT_NAME_SP:
         {	
			#ifdef __COMPILE_TIMERS__
				nb_call_prox_sp++;
				t_prox_sp.start();
			#endif
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
			/*
            #if (PROX == 0)
			faust_mat S_back_up=S[ind_fact];
			faust_mat S1,S2;
			prox_sp(S[ind_fact], const_int->getParameter());
			S1=S[ind_fact];
			S[ind_fact]=S_back_up;
			prox_sp_old_old(S[ind_fact], const_int->getParameter());
			S2=S[ind_fact];
			faust_real seuil = 0.000001;
			if (!S1.isEqual(S2,seuil))
			{	
				S_back_up.write_into_file("fail_prox_sp.dat");
				//write_faust_mat_into_matfile(S_back_up,"fail_prox_sp.dat","C");
				cerr<<"erreur prox_sp k= :"<< const_int->getParameter()<<endl;
				exit( EXIT_FAILURE); 
			}
			#endif
			
			
			#if (PROX == 1)
			cout<<"new"<<endl;
			prox_sp_old_old(S[ind_fact], const_int->getParameter());
			#endif
			
			#if (PROX == 2)
			prox_sp(S[ind_fact], const_int->getParameter());
			#endif
			*/

				prox_sp(dense_S, const_int->getParameter());
	

			
			#ifdef __COMPILE_TIMERS__
			t_prox_sp.stop();
			#endif
			
			
         }
         break;

         case CONSTRAINT_NAME_SPCOL:
         {	
			#ifdef __COMPILE_TIMERS__
				nb_call_prox_spcol++;
				t_prox_spcol.start();
				
			#endif
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);

            prox_spcol(dense_S, const_int->getParameter());

			#ifdef __COMPILE_TIMERS__
				t_prox_spcol.stop();
			#endif
         }
         break;

         case CONSTRAINT_NAME_SPLIN:
         {	

			

			#ifdef __COMPILE_TIMERS__
			nb_call_prox_splin++;
			t_prox_splin.start();
			#endif
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
			/*#if (PROX == 0)
				//cout<<"comp"<<endl;		
			faust_mat S_back_up=S[ind_fact];
			faust_mat S1,S2;
			prox_splin(S[ind_fact], const_int->getParameter());
			S1=S[ind_fact];
			S[ind_fact]=S_back_up;
			old_splin(S[ind_fact], const_int->getParameter());
			S2=S[ind_fact];
			faust_real seuil=0.000001;
			if (!S1.isEqual(S2,seuil))
			{	
				S_back_up.write_into_file("fail_prox_splin.dat");
				//write_faust_mat_into_matfile(S_back_up,"fail_prox_splin.mat","C");
				cerr<<"erreur prox_splin k = :"<<const_int->getParameter()<<endl;
				exit( EXIT_FAILURE); 
			}
			#endif
			#if (PROX == 1)
				//cout<<"old"<<endl;
				old_splin(S[ind_fact], const_int->getParameter());
			#endif
			
			#if (PROX == 2)
				//cout<<"new"<<endl;
				prox_splin(S[ind_fact], const_int->getParameter());
			#endif*/

				prox_splin(dense_S, const_int->getParameter());

			
			#ifdef __COMPILE_TIMERS__
				t_prox_splin.stop();
			#endif
			
         }
         break;

         case CONSTRAINT_NAME_NORMCOL:
         {	
			#ifdef __COMPILE_TIMERS__
				nb_call_prox_normcol++;
				t_prox_normcol.start();
			#endif
            const faust_constraint_real* const_real = dynamic_cast<const faust_constraint_real*>(const_vec[ind_fact]);
            prox_normcol(dense_S, const_real->getParameter());
			#ifdef __COMPILE_TIMERS__
				t_prox_normcol.stop();
			#endif
         }
         break;

         case CONSTRAINT_NAME_SPLINCOL:
         {	

            const faust_constraint_real* const_real = dynamic_cast<const faust_constraint_real*>(const_vec[ind_fact]);
            //prox_splincol(S[ind_fact], const_real->getParameter());

         }
         break;

         case CONSTRAINT_NAME_L0PEN:
         {	

            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_l0pen(S[ind_fact], sqrt(2*const_int->getParameter()/c));

         }
         break;

         case CONSTRAINT_NAME_L1PEN:
         {	

            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_l1pen(S[ind_fact], const_int->getParameter()/c);

         }
         break;

         case CONSTRAINT_NAME_WAV:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_wav(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SP_POS:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp_pos(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_BLKDIAG:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SPLIN_TEST:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SUPP:
         {
            const faust_constraint_mat* const_mat = dynamic_cast<const faust_constraint_mat*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_mat->getParameter());
         }
         break;

         case CONSTRAINT_NAME_NORMLIN:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_TOEPLITZ:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_int->getParameter());
         }
         break;

         default:
            cerr << "error in palm4MSA::compute_projection : unknown name of constraint" << endl;
            exit(EXIT_FAILURE);

      }
	  sparse_S = dense_S;
   }
   isProjectionComputed = true;
#ifdef __COMPILE_TIMERS__
t_global_compute_projection.stop();
t_local_compute_projection.stop();
#endif
}
	
	



void sparsePalm4MSA::update_LSR()
{
	#ifdef __COMPILE_TIMERS__
		t_global_update_LSR.start();
		t_local_update_LSR.start();
	#endif
	if(isUpdateWayR2L)
	{
		R.push_first(sparse_S);
		L.pop_back(sparse_S);
		ind_fact--;
	}else
	{
		L.push_back(sparse_S);
		R.pop_first(sparse_S);
		ind_fact++;
	}
	#ifdef __COMPILE_TIMERS__
		t_global_update_LSR.stop();
		t_local_update_LSR.stop();
	#endif
}


void sparsePalm4MSA::compute_lambda()
{
	#ifdef __COMPILE_TIMERS__
		t_global_compute_lambda.start();
		t_local_compute_lambda.start();
	#endif
	faust_mat Xhat_t_X,Xhat_t_Xhat_product;
	faust_core Xhat_t_Xhat;
	
	//at the end of next_step,
	//if updatewayR2L is true, R=Xhat
	//if updatewayR2L is false, L=Xhat
	if (isUpdateWayR2L)
	{	
		Xhat_t_Xhat=R;
		//compute (Xhat') * (X) instead of (X')*Xhat is less complex because
		//Xhat is a faust_core, 
		multiply(R,data,Xhat_t_X,1.0,'T','R');
	}else
	{	
		Xhat_t_Xhat=L;
		multiply(L,data,Xhat_t_X,1.0,'T','R');	
	}
	
	Xhat_t_Xhat.transpose();
	if (isUpdateWayR2L)
	{
		Xhat_t_Xhat.multiply(R);
	}else
	{
		Xhat_t_Xhat.multiply(L);
	}
	Xhat_t_Xhat_product=Xhat_t_Xhat.get_product();
	
	lambda = Xhat_t_X.trace()/Xhat_t_Xhat_product.trace();
	//cout<<"lambda :"<<lambda<<endl;
	
	#ifdef __COMPILE_TIMERS__
		t_global_compute_lambda.stop();
		t_local_compute_lambda.stop();
	#endif	
	
	
}

void sparsePalm4MSA::last_update()
{
	#ifdef __COMPILE_TIMERS__
		t_global_last_update.start();
		t_local_last_update.start();
	#endif
	if (isUpdateWayR2L)
	{
		R.pop_back(sparse_S);
		L=R;
		R.clear();
		ind_fact=nb_fact-1;
	}else
	{
		L.pop_first(sparse_S);
		R=L;
		L.clear();
		ind_fact = 0;
	}
	#ifdef __COMPILE_TIMERS__
		t_global_last_update.stop();
		t_local_last_update.stop();
	#endif	
}






void sparsePalm4MSA::check_constraint_validity()
{
#ifdef __COMPILE_TIMERS__
t_global_check.start();
t_local_check.start();
#endif

   if (nb_fact != (L.size()+1+R.size()) )
   {
      cerr << "Error in sparsePalm4MSA::check_constraint_validity : Wrong initialization: params.nfacts and params.init_facts are in conflict" << endl;
      exit(EXIT_FAILURE);
   }
   if (nb_fact != const_vec.size())
   {
	 cerr << "Error in sparsePalm4MSA::check_constraint_validity : Wrong initialization: nb_fact and the size of const_vec are in conflict" << endl;
      exit(EXIT_FAILURE);   
   }   
#ifdef __COMPILE_TIMERS__
t_global_check.stop();
t_local_check.stop();
#endif
}


void sparsePalm4MSA::init_fact()
{
	#ifdef __COMPILE_TIMERS__
		t_global_init_fact.start();
		t_local_init_fact.start();
	#endif
	
	

	L.clear();
	R.clear();
	if (isUpdateWayR2L)
	{	
		ind_fact = nb_fact-1;
		for (int i=0;i<ind_fact-1;i++)
		{
			sparse_S.resize(const_vec[i]->getRows(),const_vec[i]->getCols());
			sparse_S.setEyes();
			L.push_back(sparse_S);
		}
		
		sparse_S.resize(const_vec[ind_fact]->getRows(),const_vec[ind_fact]->getCols());
		sparse_S.setZeros();
		R.clear();
	}
	else
	{	
		ind_fact = 0;
		for (int i=1;i<nb_fact;i++)
		{
			sparse_S.resize(const_vec[i]->getRows(),const_vec[i]->getCols());
			sparse_S.setEyes();
			R.push_back(sparse_S);
		}
		//std::cout<<"init theorical sparse_S size : "<<const_vec[ind_fact]->getRows()<<" "<<const_vec[ind_fact]->getCols()<<std::endl;
		sparse_S.resize(const_vec[ind_fact]->getRows(),const_vec[ind_fact]->getCols());
		sparse_S.setZeros();
		L.clear();
	}
	//std::cout<<"init sparse_S size : "<<sparse_S.getNbRow()<<" "<<sparse_S.getNbCol()<<std::endl;
	#ifdef __COMPILE_TIMERS__
		t_global_init_fact.stop();
		t_local_init_fact.stop();
	#endif
}


void sparsePalm4MSA::init_fact_from_palm(const sparsePalm4MSA& palm2, bool isFactSideLeft)
{	
#ifdef __COMPILE_TIMERS__
t_global_init_fact_from_palm.start();
t_local_init_fact_from_palm.start();
#endif
	#ifdef PRINT
	cout<<"INIT_FACT_FROM_PALM : "<<endl;
	cout<<"isFactSideLeft : "<<isFactSideLeft<<endl;
	cout<<"palm2.updatewayR2L : "<<palm2.isUpdateWayR2L<<endl;
	cout<<"palm_global.updatewayR2L : "<<isUpdateWayR2L<<endl;
	Display();
	#endif
	
   
	if (palm2.nb_fact != 2)
	{
		cerr<<"sparsePalm4MSA::init_fact_from_palm : argument palm2 must contain 2 factors"<<endl;
		exit(EXIT_FAILURE);
	}
	if(!isConstraintSet)
	{
		cerr << "Error in sparsePalm4MSA::init_fact_from_palm : constrainst must be set before calling init_fact_from_palm" << endl;
		exit(EXIT_FAILURE);
	}
	if(palm2.isUpdateWayR2L != isUpdateWayR2L)
	{
		cerr << "Error in sparsePalm4MSA::init_fact_from_palm : the 2 sparsePalm4MSA hasn't the same isUpdateWayR2L" << endl;
		exit(EXIT_FAILURE);
	}
	faust_spmat residuum;
		
	// left_part is update 
	if (nb_fact <= 2)
	{
		L=palm2.L;
		sparse_S=palm2.sparse_S;
		R=palm2.R;
	}else
	{
		
	
	
		if(isFactSideLeft)
		{	
			
		//palmGlobal : L S []
			if(isUpdateWayR2L)
			{	
				L.pop_first(residuum);
				palm2.L.pop_first(residuum);
				L.push_first(palm2.sparse_S);
				L.push_first(residuum);	
			}else
			{
				R.pop_first(residuum);
				palm2.R.pop_first(residuum);
				sparse_S = palm2.sparse_S;
				R.push_first(residuum);
			}
		
		}else
		{
			if (isUpdateWayR2L)
			{
				L.pop_back(residuum);
				palm2.L.pop_first(residuum);
				sparse_S = palm2.sparse_S;
				L.push_back(residuum);
			}else
			{
				R.pop_back(residuum);
				palm2.R.pop_first(residuum);
				R.push_back(palm2.sparse_S);
				R.push_back(residuum);	
			}
			
			
		}
	}
	#ifdef __COMPILE_TIMERS__
		t_global_init_fact_from_palm.stop();
		t_local_init_fact_from_palm.stop();
	#endif
	
	
}

#ifdef __COMPILE_TIMERS__

faust_timer sparsePalm4MSA::t_global_compute_projection;
faust_timer sparsePalm4MSA::t_global_compute_grad_over_c;
faust_timer sparsePalm4MSA::t_global_compute_c;
faust_timer sparsePalm4MSA::t_global_compute_lambda;
faust_timer sparsePalm4MSA::t_global_update_LSR;
faust_timer sparsePalm4MSA::t_global_last_update;
faust_timer sparsePalm4MSA::t_global_check;
faust_timer sparsePalm4MSA::t_global_init_fact;
faust_timer sparsePalm4MSA::t_global_next_step;
faust_timer sparsePalm4MSA::t_global_init_fact_from_palm;



faust_timer sparsePalm4MSA::t_prox_const;
faust_timer sparsePalm4MSA::t_prox_sp;
faust_timer sparsePalm4MSA::t_prox_spcol;
faust_timer sparsePalm4MSA::t_prox_splin;
faust_timer sparsePalm4MSA::t_prox_normcol;

int sparsePalm4MSA::nb_call_prox_const;
int sparsePalm4MSA::nb_call_prox_sp;
int sparsePalm4MSA::nb_call_prox_spcol;
int sparsePalm4MSA::nb_call_prox_splin;
int sparsePalm4MSA::nb_call_prox_normcol;






void sparsePalm4MSA::init_local_timers()
{
t_local_compute_projection.reset();
t_local_compute_grad_over_c.reset();
t_local_compute_lambda.reset();
t_local_update_LSR.reset();
t_local_last_update.reset();
t_local_check.reset();
t_local_init_fact.reset();
t_local_next_step.reset();
t_local_init_fact_from_palm.reset();
}







void sparsePalm4MSA::print_global_timers()const
{
   cout << "timers in palm4MSA : " << endl;
   cout << "t_global_next_step           = " << t_global_next_step.get_time()           << " s for "<< t_global_next_step.get_nb_call()           << " calls" << endl;
   cout << "t_compute_grad_over_c = " << t_global_compute_grad_over_c.get_time() << " s for "<< t_global_compute_grad_over_c.get_nb_call()            << " calls of grad" << endl;
     cout << "t_global_compute_c = " << t_global_compute_c.get_time() << " s for "<< t_global_compute_c.get_nb_call() << " calls" << endl;
   cout << "t_global_compute_lambda      = " << t_global_compute_lambda.get_time()      << " s for "<< t_global_compute_lambda.get_nb_call()      << " calls" << endl;
   cout << "t_global_compute_projection  = " << t_global_compute_projection.get_time()  << " s for "<< t_global_compute_projection.get_nb_call()  << " calls" << endl<<endl;
    cout << "t_global_update_LSR  = " << t_global_update_LSR.get_time()  << " s for "<< t_global_update_LSR.get_nb_call()  << " calls" << endl<<endl;
	cout << "t_global_last_update  = " << t_global_last_update.get_time()  << " s for "<< t_global_last_update.get_nb_call()  << " calls" << endl<<endl;

}


void sparsePalm4MSA::print_local_timers()const
{
   cout << "timers in palm4MSA : " << endl;
   cout << "t_local_next_step           = " << t_local_next_step.get_time()           << " s for "<< t_local_next_step.get_nb_call()           << " calls" << endl;
   cout << "t_compute_grad_over_c = " << t_local_compute_grad_over_c.get_time() << " s for "<< t_local_compute_grad_over_c.get_nb_call()            << " calls of grad" << endl;
     cout << "t_local_compute_c = " << t_local_compute_c.get_time() << " s for "<< t_local_compute_c.get_nb_call() << " calls" << endl;
   cout << "t_local_compute_lambda      = " << t_local_compute_lambda.get_time()      << " s for "<< t_local_compute_lambda.get_nb_call()      << " calls" << endl;
   cout << "t_local_compute_projection  = " << t_local_compute_projection.get_time()  << " s for "<< t_local_compute_projection.get_nb_call()  << " calls" << endl<<endl;
    cout << "t_local_update_LSR  = " << t_local_update_LSR.get_time()  << " s for "<< t_local_update_LSR.get_nb_call()  << " calls" << endl<<endl;
	cout << "t_local_last_update  = " << t_local_last_update.get_time()  << " s for "<< t_local_last_update.get_nb_call()  << " calls" << endl<<endl;

}

void sparsePalm4MSA::print_prox_timers() const
{
/*	cout << "prox timers in palm4MSA : " << endl;
   cout << "total t_prox_const  =  " << t_prox_const.get_time()  << " s for "<< nb_call_prox_const  << " calls" << endl;
   cout << "total t_prox_sp  =  " << t_prox_sp.get_time()  << " s for "<< nb_call_prox_sp  << " calls" << endl;
   cout << "total t_prox_spcol  =  " << t_prox_spcol.get_time()  << " s for "<< nb_call_prox_spcol  << " calls" << endl;
   cout << "total t_prox_splin  =  " << t_prox_splin.get_time()  << " s for "<< nb_call_prox_splin  << " calls" << endl;
   cout << "total t_prox_normcol  =  " << t_prox_normcol.get_time()  << " s for "<< nb_call_prox_normcol  << " calls" << endl;	
*/}




#endif
