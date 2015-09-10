#include "faust_timer.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"
using namespace std;

faust_timer::faust_timer() : isRunning(false), result(0.0f), nbCall(0){}

void faust_timer::start()
{ 
   if(isRunning)
   {
	  ErrorDisplay("faust_timer::start : timer is already started.\n");
   }
   clock_gettime(CLOCK_MONOTONIC, &debut);
   isRunning = true;
   nbCall ++;
}

void faust_timer::stop()
{
   if(!isRunning)
   {
	  ErrorDisplay("faust_timer::stop : timer must be started before stopping it\n");
   }
   struct timespec fin;
   clock_gettime(CLOCK_MONOTONIC, &fin);
   result += (fin.tv_sec -debut.tv_sec) + (fin.tv_nsec-debut.tv_nsec)/1000000000.0;
   isRunning = false;
}


void faust_timer::reset()
{
   result=0.0f;
   nbCall = 0;
   if(isRunning)
   {
	  WarningDisplay("faust_timer::reset : timer has been reset while it was running\n");
      clock_gettime(CLOCK_MONOTONIC, &debut);
   }
}

float faust_timer::get_time()
{
   if(isRunning)
   {
      struct timespec fin;
      clock_gettime(CLOCK_MONOTONIC, &fin);
      result += (fin.tv_sec -debut.tv_sec) + (fin.tv_nsec-debut.tv_nsec)/1000000000.0;
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
      struct timespec fin;
      clock_gettime(CLOCK_MONOTONIC, &fin);
      result += (fin.tv_sec -debut.tv_sec) + (fin.tv_nsec-debut.tv_nsec)/1000000000.0;
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



