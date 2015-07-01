#ifndef __FAUST_TIMER_H__
#define __FAUST_TIMER_H__

#include <ctime>

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
      struct timespec debut;
      long int nbCall;
     
};


#endif
