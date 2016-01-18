#ifndef INIT_PARAMS_FROM_XML_H
#define INIT_PARAMS_FROM_XML_H

#include "faust_params.h"
#include <vector>
#include <libxml/parser.h>

template<typename T>
void init_params_from_xml(const char * filename,faust_params<T> & params);



template<typename T>
void add_constraint(std::vector<const faust_constraint_generic*> & consS,char* type, char * parameter, char* dim1,char* dim2);
#include "faust_init_params_from_xml.hpp"

#endif