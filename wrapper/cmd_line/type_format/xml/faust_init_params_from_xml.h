#ifndef INIT_PARAMS_FROM_XML_H
#define INIT_PARAMS_FROM_XML_H

#include "faust_Params.h"
#include "faust_ParamsPalm.h"
#include <vector>
#include <libxml/parser.h>





	template<typename FPP,Device DEVICE>
	void init_params_from_xml(const char * filename,Faust::Params<FPP,DEVICE> & params);
	template<typename FPP,Device DEVICE>
	void init_palm_params_from_xml(const char * filename,Faust::ParamsPalm<FPP,DEVICE> & params);



	template<typename FPP,Device DEVICE>
	void add_constraint(std::vector<const Faust::ConstraintGeneric<FPP,DEVICE> *> & consS,char* type, char * parameter, char* dim1,char* dim2);




#include "faust_init_params_from_xml.hpp"

#endif
