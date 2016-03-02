#ifndef __FAUST_CU_TIMER_H__
#define __FAUST_CU_TIMER_H__

#include <ctime>
#include "faust_constant.h"
#include <cuda_runtime.h>


class faust_cu_timer
{
   public:
      faust_cu_timer();
      faust_cu_timer(const cudaStream_t);
      void start();
      void stop();
      void reset();
      float get_time()const;
      float get_time();
      faust_unsigned_int get_nb_call()const;
      faust_unsigned_int get_nb_call();
	  static const char * class_name; 


   private:
      bool isRunning;
      float result;

      cudaEvent_t debut;
      cudaEvent_t fin;
      cudaStream_t stream;
      faust_unsigned_int nbCall;
     
};


#endif
