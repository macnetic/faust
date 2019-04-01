/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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

#include <cstring>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "faust_exception.h"

//#include "faust_constant.h"


#ifdef COMPILE_MEX
	#include <mex.h>
#endif
void handleTxt(const char texte [],std::stringstream & ss, va_list ap);


void handleTxt(const char texte [],std::stringstream & ss, va_list ap)
{
	int tailleChaine = 0, i = 0;
	tailleChaine = strlen(texte);
	while (i < tailleChaine)
	{
        if(texte[i] == '%'){
            i++; // on se place au caractère suivant pour savoir ce que c'est
            switch(texte[i])
			{	
				// cas int
				case 'i' :
					ss << (va_arg(ap,int));
				
				// cas faust unsigned_int
                case 'd' :
                    //if (sizeof(faust_unsigned_int) == sizeof(unsigned long int)){
						ss<<(va_arg(ap,unsigned long int));
					/*}else
					{
						handleError("faust_exception : handleTxt : unsupported faust_unsigned_int");
					}*/
					break;
				
				// cas faust_real	
				case 'f' :
					ss<<(va_arg(ap,double));
					break;
				case 's' : 
					ss<<(va_arg(ap,char*));
					break;
				
            }
		}
        else{
			ss<<texte[i];
            //printf("%c",texte[i]);
			
        }
		i++;
	}
	ss<<"\n";
}


void handleWarning(const char texte [], ...)
{
   
    va_list ap;
    va_start(ap, texte);
	std::stringstream ss;
	
	
   
	#ifdef COMPILE_MEX
		ss<<" MATLAB ";
	#else
		ss<<" C++ Warning in ";
	#endif
	
	handleTxt(texte,ss,ap);
	
    va_end(ap);
	
	#ifdef COMPILE_MEX
		mexWarnMsgTxt(ss.str().c_str());
	#else
		//printf(ss.str().c_str());
		std::cerr<<ss.str().c_str();
		
	#endif
}

/*void handleError(const char texte[], ...)
{
    
    
	int tailleChaine = 0, i;
    va_list ap;
    va_start(ap, texte);
	std::stringstream ss;
	i=0;	

    tailleChaine = strlen(texte);
	#ifdef COMPILE_MEX
		ss<<" MATLAB ";
	#else
		ss<<" C++ Error in";
	#endif
	
	handleTxt(texte,ss,ap);
	
    va_end(ap);
	
	#ifdef COMPILE_MEX
		mexErrMsgTxt(ss.str().c_str());
	#else
		std::cerr<<ss.str().c_str();
		exit( EXIT_FAILURE);
	#endif
}*/




