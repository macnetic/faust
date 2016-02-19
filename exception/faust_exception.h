#ifndef __FAUST_EXCEPTION_H__
#define __FAUST_EXCEPTION_H__

#include <stdarg.h>
#include <cstdlib>
#include <stdexcept>


void handleError(const char* classe_name , const  char* txt);
//void handleError(const  char* txt);
void handleWarning(const char texte [], ...);
//void faust_exception_throw(const char* texte);


#endif
