#include <cstring>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "faust_exception.h"


#ifdef COMPILE_MEX
	#include <mex.h>
#endif
void HandleTxt(const char texte [],std::stringstream & ss, va_list ap);

void HandleTxt(const char texte [],std::stringstream & ss, va_list ap)
{
	int tailleChaine = 0, i = 0;
	tailleChaine = strlen(texte);
	while (i < tailleChaine)
	{
        if(texte[i] == '%'){
            i++; // on se place au caractère suivant pour savoir ce que c'est
            switch(texte[i])
			{
				
                case 'd' :
                    //printf("%d", va_arg(ap, int));	
                    ss<<(va_arg(ap, int));
					break;
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


void WarningDisplay(const char texte [], ...)
{
   
    va_list ap;
    va_start(ap, texte);
	std::stringstream ss;
	
	
   
	#ifdef COMPILE_MEX
		ss<<" MATLAB ";
	#else
		ss<<" C++ Warning in ";
	#endif
	
	HandleTxt(texte,ss,ap);
	
    va_end(ap);
	
	#ifdef COMPILE_MEX
		mexWarnMsgTxt(ss.str().c_str());
	#else
		//printf(ss.str().c_str());
		std::cerr<<ss.str().c_str();
		
	#endif
}

void ErrorDisplay(const char texte[], ...)
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
	
	HandleTxt(texte,ss,ap);
	
    va_end(ap);
	
	#ifdef COMPILE_MEX
		mexErrMsgTxt(ss.str().c_str());
	#else
		std::cerr<<ss.str().c_str();
		exit( EXIT_FAILURE);
	#endif
}