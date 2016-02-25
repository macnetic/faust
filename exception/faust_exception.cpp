
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




