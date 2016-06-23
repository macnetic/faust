#ifndef __FAUST_EXCEPTION_H__
#define __FAUST_EXCEPTION_H__

#include <stdarg.h>
#include <cstdlib>
#include <stdexcept>
#include <sstream>



void handleWarning(const char texte [], ...);


/// macro pour la gestion des errors (__LINE__ et __FILE__)
#define handleError(m_className,message) do {                                \
		std::stringstream complete_message;		\
        complete_message<<"Error in file " <<__FILE__<<std::endl<<" line "<<__LINE__<< " , class  " <<m_className<<" : "<<std::endl<<message; \
	    throw std::logic_error(complete_message.str());} while (0)

#endif
