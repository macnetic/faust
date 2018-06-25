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
#ifndef __FAUST_TIMER_H__
#define __FAUST_TIMER_H__

#include <ctime>
#include "faust_constant.h"

#if defined(_WIN32)
   #include <windows.h>
#elif defined(__MACH__)
   #include <mach/clock.h>
   #include <mach/mach_time.h>
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
            static const char * m_className;

        private:
            //! \brief isRunning : TRUE value indicate that a time is running. <br>
            //! FALSE value indicates no-time running. <br>
            bool isRunning;
            float result;
            #if defined(__linux__)
                struct timespec debut;
	    #elif defined(__MACH__)
		uint64_t debut;
		double conversion_factor;
            #elif defined(_WIN32)
                LARGE_INTEGER debut;
                LARGE_INTEGER frequency;
            #endif
            faust_unsigned_int nbCall;

    };



}

#endif
