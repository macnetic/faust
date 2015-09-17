#include "faust_timer.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"
using namespace std;


faust_timer::faust_timer() : isRunning(false), result(0.0f), nbCall(0)
{


#if defined(_WIN32)  
   QueryPerformanceFrequency(&frequency); 
#elif !defined(__linux__)
   cerr << "Error in faust_timer::faust_timer : OS not supported" <<endl;
   exit(EXIT_FAILURE);
#endif

}

void faust_timer::start()
{ 
   if(isRunning)
   {
	  ErrorDisplay("faust_timer::start : timer is already started.\n");
   }
   #if defined(__linux__)
      clock_gettime(CLOCK_MONOTONIC, &debut);
   #elif defined(_WIN32)  
      QueryPerformanceCounter(&debut);
   #endif

   isRunning = true;
   nbCall ++;
}

void faust_timer::stop()
{
   if(!isRunning)
   {
	  ErrorDisplay("faust_timer::stop : timer must be started before stopping it\n");
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


void faust_timer::reset()
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
      WarningDisplay("faust_timer::reset : timer has been reset while it was running\n");
         
   }
}

float faust_timer::get_time()
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
          

         WarningDisplay("faust_timer::get_time : timer has not been stopped\n");
   }
   return result;
}

float faust_timer::get_time()const
{
   if(isRunning)
   {
	  ErrorDisplay("faust_timer::get_time : timer has not been stopped\n");		

   }
   return result;
}


long int faust_timer::get_nb_call()
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
         WarningDisplay("faust_timer::get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}


long int faust_timer::get_nb_call()const
{
   if(isRunning)
   {
	  ErrorDisplay("faust_timer::get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}



