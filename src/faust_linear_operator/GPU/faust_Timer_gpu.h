#ifndef __FAUST_TIMER_GPU_H__
#define __FAUST_TIMER_GPU_H__

#include <ctime>
#include "faust_constant_gpu.h"
#include <cuda_runtime.h>

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{


    class Timer_gpu
    {
       public:
          Timer_gpu();
          Timer_gpu(const cudaStream_t);
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
}


#endif
