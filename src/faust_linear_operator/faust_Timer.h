#ifndef __FAUST_TIMER_H__
#define __FAUST_TIMER_H__

#include <ctime>
#include "faust_constant.h"

#if defined(_WIN32)
   #include <windows.h>
#endif

/*! \class Faust::Timer
* \brief The class Faust::Timer propose various Transform to evaluate time processing.
*/

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{
    class Timer
    {
        public:
            Timer();
            void start();
            void stop();
            void reset();
            float get_time()const;
            float get_time();
            faust_unsigned_int get_nb_call()const;
            faust_unsigned_int get_nb_call();
            static const char * class_name;

        private:
            //! \brief isRunning : TRUE value indicate that a time is running. <br>
            //! FALSE value indicates no-time running. <br>
            bool isRunning;
            float result;
            #if defined(__linux__)
                struct timespec debut;
            #elif defined(_WIN32)
                LARGE_INTEGER debut;
                LARGE_INTEGER frequency;
            #endif
            faust_unsigned_int nbCall;

    };

}

#endif
