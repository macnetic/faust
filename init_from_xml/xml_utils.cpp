#include "xml_utils.h"
#include <iostream>

std::vector<xmlChar*> get_content(xmlChar * expression, xmlXPathContextPtr ctxt)
{
	std::vector<xmlChar*> vector_content;

	xmlXPathObjectPtr xpathRes = xmlXPathEvalExpression(expression, ctxt);

	if (xpathRes != NULL)
	{
		if (xpathRes->nodesetval != NULL)	
		{
			for (int i=0;i<xpathRes->nodesetval->nodeNr;i++)
			{
				vector_content.push_back(xpathRes->nodesetval->nodeTab[i]->content);
			}
		}
	}

	xmlXPathFreeObject(xpathRes);
	return vector_content;
	
}