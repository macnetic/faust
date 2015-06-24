#include "faust_timer.h"
#include <iostream>
#include <cstdlib>
using namespace std;

faust_timer::faust_timer() : isRunning(false), result(0.0f), nbCall(0){}

void faust_timer::start()
{ 
   if(isRunning)
   {
      cerr << "Error in faust_timer::start : timer is already started." << endl;
      exit(EXIT_FAILURE);
   }
   clock_gettime(CLOCK_MONOTONIC, &debut);
   isRunning = true;
   nbCall ++;
}

void faust_timer::stop()
{
   if(!isRunning)
   {
      cerr << "Error in faust_timer::stop : timer must be started before stopping it" << endl;
      exit(EXIT_FAILURE);
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
      cerr << "Warning in faust_timer::reset : timer has been reset while it was running" << endl;
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
      cerr << "Warning in faust_timer::get_time : timer has not been stopped" << endl;
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
      cerr << "Warning in faust_timer::get_nb_call : timer has not been stopped" << endl;
   }
   return nbCall;
}



