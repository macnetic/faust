#include <stdlib.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <iostream>
#include <vector>



std::vector<xmlChar*> get_content(xmlChar * expression, xmlXPathContextPtr ctxt);
xmlXPathContextPtr get_context(const char * filename);
