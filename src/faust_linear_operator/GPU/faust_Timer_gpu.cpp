/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
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
#include "faust_Timer_gpu.h"
#include "faust_cuda.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"
using namespace std;


Faust::Timer_gpu::Timer_gpu() : isRunning(false), result(0.0f), nbCall(0), stream(0)
{
   faust_cudaEventCreate(&debut);
   faust_cudaEventCreate(&fin);
}

Faust::Timer_gpu::Timer_gpu(const cudaStream_t stream_) : isRunning(false), result(0.0f), nbCall(0), stream(stream_)
{
   faust_cudaEventCreate(&debut);
   faust_cudaEventCreate(&fin);
}


void Faust::Timer_gpu::start()
{
   if(isRunning)
   {
	  handleError(m_className,"Faust::Timer_gpu::start : timer is already started.\n");
   }
   faust_cudaEventRecord(debut, stream);
   isRunning = true;
   nbCall ++;
}

void Faust::Timer_gpu::stop()
{
   if(!isRunning)
   {
	  handleError(m_className,"stop : timer must be started before stopping it\n");
   }
   faust_cudaEventRecord(fin, stream);
   faust_cudaEventSynchronize(fin);
   float ms = 0.0f;
   faust_cudaEventElapsedTime(&ms, debut, fin);
   result += ms/1000.0;

   isRunning = false;
}


void Faust::Timer_gpu::reset()
{
   result=0.0f;
   nbCall = 0;
   if(isRunning)
   {
      faust_cudaEventRecord(debut, stream);
      cerr<<m_className<<"reset : timer has been reset while it was running"<<endl;
      nbCall++;
   }
}

float Faust::Timer_gpu::get_time()
{
   if(isRunning)
   {
      faust_cudaEventRecord(fin, stream);
      faust_cudaEventSynchronize(fin);
      float ms = 0.0f;
      faust_cudaEventElapsedTime(&ms, debut, fin);
      result += ms/1000.0;

      handleError(m_className,"get_time : timer has not been stopped");
   }
   return result;
}

float Faust::Timer_gpu::get_time()const
{
   if(isRunning)
   {
	  handleError(m_className,"get_time : timer has not been stopped");
   }
   return result;
}


faust_unsigned_int Faust::Timer_gpu::get_nb_call()
{
   if(isRunning)
   {
      faust_cudaEventRecord(fin, stream);
      faust_cudaEventSynchronize(fin);
      float ms = 0.0f;
      faust_cudaEventElapsedTime(&ms, debut, fin);
      result += ms/1000.0;

      handleError(m_className,"get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}


faust_unsigned_int Faust::Timer_gpu::get_nb_call()const
{
   if(isRunning)
   {
	  handleError(m_className,"get_nb_call : timer has not been stopped\n");
   }
   return nbCall;
}

const char * Faust::Timer_gpu::m_className="Faust::Timer_gpu::";

