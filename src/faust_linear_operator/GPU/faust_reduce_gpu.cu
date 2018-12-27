/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2018):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
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
#include "faust_reduce_gpu.h"
#include <limits>

/*
#include "faust_cu_vec.h"
#include "faust_cu_matrix.h"

template<>
double faust_cu_reduce<double>(const faust_cu_vec<double>& v)
{
   thrust::device_ptr<double> dev_ptr(v.getData());
   double somme = thrust::reduce(dev_ptr, dev_ptr+v.size());
   return somme;
}

template<>
float faust_cu_reduce<float>(const faust_cu_vec<float>& v)
{
   thrust::device_ptr<float> dev_ptr(v.getData());
   float somme = thrust::reduce(dev_ptr, dev_ptr+v.size());
   return somme;
}


template<>
double faust_cu_reduce<double>(const faust_cu_matrix<double>& M)
{
   thrust::device_ptr<double> dev_ptr(M.getData());
   double somme = thrust::reduce(dev_ptr, dev_ptr+M.getNbRow()*M.getNbCol());
   return somme;
}

template<>
float faust_cu_reduce<float>(const faust_cu_matrix<float>& M)
{
   thrust::device_ptr<float> dev_ptr(M.getData());
   float somme = thrust::reduce(dev_ptr, dev_ptr+M.getNbRow()*M.getNbCol());
   return somme;
}*/


template<>
double faust_cu_sum<double>(const double* data, const int nb_el)
{
   thrust::device_ptr<const double> dev_ptr(data);
   const double somme = thrust::reduce(dev_ptr, dev_ptr+nb_el);
   return somme;
}
template<>
float faust_cu_sum<float>(const float* data, const int nb_el)
{
   thrust::device_ptr<const float> dev_ptr(data);
   const float somme = thrust::reduce(dev_ptr, dev_ptr+nb_el);
   return somme;
}

template<>
double faust_cu_max<double>(const double* data, const int nb_el)
{
   thrust::device_ptr<const double> dev_ptr(data);
   const double maxi = thrust::reduce(dev_ptr, dev_ptr+nb_el, -1.0e300,thrust::maximum<double>());
   return maxi;
}
template<>
float faust_cu_max<float>(const float* data, const int nb_el)
{
   thrust::device_ptr<const float> dev_ptr(data);
   const float maxi = thrust::reduce(dev_ptr, dev_ptr+nb_el, -1.0e300,thrust::maximum<float>());
   return maxi;
}

template<>
double faust_cu_min<double>(const double* data, const int nb_el)
{
   thrust::device_ptr<const double> dev_ptr(data);
   const double mini = thrust::reduce(dev_ptr, dev_ptr+nb_el, 1.0e300,thrust::minimum<double>());
   return mini;
}
template<>
float faust_cu_min<float>(const float* data, const int nb_el)
{
   thrust::device_ptr<const float> dev_ptr(data);
   const float mini = thrust::reduce(dev_ptr, dev_ptr+nb_el, 1.0e300,thrust::minimum<float>());
   return mini;
}

template<>
double faust_cu_norm<double>(const double* data, const int nb_el)
{
   thrust::device_ptr<const double> dev_ptr(data);
   const double frob_norm = std::sqrt(thrust::inner_product(dev_ptr, dev_ptr+nb_el, dev_ptr, 0.0, thrust::plus<double>() ,thrust::multiplies<double>()));
   return frob_norm;
}
template<>
float faust_cu_norm<float>(const float* data, const int nb_el)
{
   thrust::device_ptr<const float> dev_ptr(data);
   const float frob_norm = std::sqrt(thrust::inner_product(dev_ptr, dev_ptr+nb_el, dev_ptr, 0.0f, thrust::plus<float>() ,thrust::multiplies<float>()));
   return frob_norm;
}
