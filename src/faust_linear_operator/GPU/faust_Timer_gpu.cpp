#include "faust_Timer_gpu.h"
#include "faust_cuda.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"
using namespace std;


Faust::Timer_gpu::Timer_gpu() : isRunning(false), result(0.0f), nbCall(0), stream(0)
{
   faust_cudaEventCreate(&debut);
   faust_cudaEventCreate(&fin);
}

Faust::Timer_gpu::Timer_gpu(const cudaStream_t stream_) : isRunning(false), result(0.0f), nbCall(0), stream(stream_)
{
   faust_cudaEventCreate(&debut);
   faust_cudaEventCreate(&fin);
}


void Faust::Timer_gpu::start()
{
   if(isRunning)
   {
	  handleError(m_className,"Faust::Timer_gpu::start : timer is already started.\n");
   }
   faust_cudaEventRecord(debut, stream);
   isRunning = true;
   nbCall ++;
}

void Faust::Timer_gpu::stop()
{
   if(!isRunning)
   {
	  handleError(m_className,"stop : timer must be started before stopping it\n");
   }
   faust_cudaEventRecord(fin, stream);
   faust_cudaEventSynchronize(fin);
   float ms = 0.0f;
   faust_cudaEventElapsedTime(&ms, debut, fin);
   result += ms/1000.0;

   isRunning = false;
}


void Faust::Timer_gpu::reset()
{
   result=0.0f;
   nbCall = 0;
   if(isRunning)
   {
      faust_cudaEventRecord(debut, stream);
      cerr<<m_className<<"reset : timer has been reset while it was running"<<endl;
      nbCall++;
   }
}

float Faust::Timer_gpu::get_time()
{
   if(isRunning)
   {
      faust_cudaEventRecord(fin, stream);
      faust_cudaEventSynchronize(fin);
      float ms = 0.0f;
      faust_cudaEventElapsedTime(&ms, debut, fin);
      result += ms/1000.0;

      handleError(m_className,"get_time : timer has not been stopped");
   }
   return result;
}

float Faust::Timer_gpu::get_time()const
{
   if(isRunning)
   {
	  handleError(m_className,"get_time : timer has not been stopped");
   }
   return result;
}


faust_unsigned_int Faust::Timer_gpu::get_nb_call()
{
   if(isRunning)
   {
      faust_cudaEventRecord(fin, stream);
      faust_cudaEventSynchronize(fin);
      float ms = 0.0f;
      faust_cudaEventElapsedTime(&ms, debut, fin);
      result += ms/1000.0;

      handleError(m_className,"get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}


faust_unsigned_int Faust::Timer_gpu::get_nb_call()const
{
   if(isRunning)
   {
	  handleError(m_className,"get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}

const char * Faust::Timer_gpu::m_className="Faust::Timer_gpu::";

