/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#include "faust_Timer.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"
using namespace std;


Faust::Timer::Timer() : isRunning(false), result(0.0f), nbCall(0)
{


#if defined(_WIN32)
   QueryPerformanceFrequency(&frequency);
#elif defined(__MACH__)
   mach_timebase_info_data_t timebase;
   mach_timebase_info(&timebase);
   conversion_factor = (double) timebase.numer / ((double) timebase.denom *1e9);
#elif defined(__linux__)
#else
   handleError(m_className,"Error in Faust::Timer::Timer : OS not supported");
#endif

}

void Faust::Timer::start()
{
   if(isRunning)
   {
	  handleError(m_className,"Faust::Timer::start : timer is already started.\n");
   }
   #if defined(__linux__)
      clock_gettime(CLOCK_MONOTONIC, &debut);
   #elif defined(__MACH__)
	debut=mach_absolute_time();
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
	  handleError(m_className,"stop : timer must be started before stopping it\n");
   }
   #if defined(__linux__)
      struct timespec fin;
      clock_gettime(CLOCK_MONOTONIC, &fin);
      result += (fin.tv_sec -debut.tv_sec) + (fin.tv_nsec-debut.tv_nsec)/1000000000.0;
   #elif defined(__MACH__)
	uint64_t fin;
	fin = mach_absolute_time();
	result+=((double) (fin-debut)) * conversion_factor;
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
      #elif defined(__MACH__)
	debut = mach_absolute_time();
	#elif defined(_WIN32)
         QueryPerformanceCounter(&debut);
      #endif
      cerr<<m_className<<"reset : timer has been reset while it was running"<<endl;

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


         handleError(m_className,"get_time : timer has not been stopped");
   }
   return result;
}

float Faust::Timer::get_time()const
{
   if(isRunning)
   {
	  handleError(m_className,"get_time : timer has not been stopped");

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
         handleError(m_className,"get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}


faust_unsigned_int Faust::Timer::get_nb_call()const
{
   if(isRunning)
   {
	  handleError(m_className,"get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}

const char * Faust::Timer::m_className="Faust::Timer::";

