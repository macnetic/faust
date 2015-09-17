#ifndef __FAUST_TIMER_H__
#define __FAUST_TIMER_H__

#include <ctime>

#if defined(_WIN32)  
   #include <windows.h>
#endif

class faust_timer
{
   public:
      faust_timer();
      void start();
      void stop();
      void reset();
      float get_time()const;
      float get_time();
      long int get_nb_call()const;
      long int get_nb_call();


   private:
      bool isRunning;
      float result;
      #if defined(__linux__) 
         struct timespec debut;
      #elif defined(_WIN32)  
         LARGE_INTEGER debut;
         LARGE_INTEGER frequency;  
      #endif
      long int nbCall;
     
};


#endif
