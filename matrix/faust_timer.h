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
      float get_time();

   private:
      bool isRunning;
      float result;
      struct timespec debut;
     
};


#endif
