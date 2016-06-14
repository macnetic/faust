#include "faust_Timer.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"
using namespace std;


Faust::Timer::Timer() : isRunning(false), result(0.0f), nbCall(0)
{


#if defined(_WIN32)
   QueryPerformanceFrequency(&frequency);
#elif !defined(__linux__)

   handleError(class_name,"Error in Faust::Timer::Timer : OS not supported");

#endif

}

void Faust::Timer::start()
{
   if(isRunning)
   {
	  handleError(class_name,"Faust::Timer::start : timer is already started.\n");
   }
   #if defined(__linux__)
      clock_gettime(CLOCK_MONOTONIC, &debut);
   #elif defined(_WIN32)
      QueryPerformanceCounter(&debut);
   #endif

   isRunning = true;
   nbCall ++;
}

void Faust::Timer::stop()
{
   if(!isRunning)
   {
	  handleError(class_name,"stop : timer must be started before stopping it\n");
   }
   #if defined(__linux__)
      struct timespec fin;
      clock_gettime(CLOCK_MONOTONIC, &fin);
      result += (fin.tv_sec -debut.tv_sec) + (fin.tv_nsec-debut.tv_nsec)/1000000000.0;
   #elif defined(_WIN32)
      LARGE_INTEGER fin;
      QueryPerformanceCounter(&fin);
      result += (fin.QuadPart - debut.QuadPart)*1000.0/frequency.QuadPart;
   #endif
   isRunning = false;
}


void Faust::Timer::reset()
{
   result=0.0f;
   nbCall = 0;
   if(isRunning)
   {
      #if defined(__linux__)
         clock_gettime(CLOCK_MONOTONIC, &debut);
      #elif defined(_WIN32)
         QueryPerformanceCounter(&debut);
      #endif
      cerr<<class_name<<"reset : timer has been reset while it was running"<<endl;

   }
}

float Faust::Timer::get_time()
{
   if(isRunning)
   {

      #if defined(__linux__)
         struct timespec fin;
         clock_gettime(CLOCK_MONOTONIC, &fin);
         result += (fin.tv_sec -debut.tv_sec) + (fin.tv_nsec-debut.tv_nsec)/1000000000.0;
      #elif defined(_WIN32)
         LARGE_INTEGER fin;
         QueryPerformanceCounter(&fin);
         result += (fin.QuadPart - debut.QuadPart)*1000.0/frequency.QuadPart;
      #endif


         handleError(class_name,"get_time : timer has not been stopped");
   }
   return result;
}

float Faust::Timer::get_time()const
{
   if(isRunning)
   {
	  handleError(class_name,"get_time : timer has not been stopped");

   }
   return result;
}


faust_unsigned_int Faust::Timer::get_nb_call()
{
   if(isRunning)
   {

      #if defined(__linux__)
         struct timespec fin;
         clock_gettime(CLOCK_MONOTONIC, &fin);
         result += (fin.tv_sec -debut.tv_sec) + (fin.tv_nsec-debut.tv_nsec)/1000000000.0;
      #elif defined(_WIN32)
         LARGE_INTEGER fin;
         QueryPerformanceCounter(&fin);
         result += (fin.QuadPart - debut.QuadPart)*1000.0/frequency.QuadPart;
      #endif
         handleError(class_name,"get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}


faust_unsigned_int Faust::Timer::get_nb_call()const
{
   if(isRunning)
   {
	  handleError(class_name,"get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}

const char * Faust::Timer::class_name="Faust::Timer::";

