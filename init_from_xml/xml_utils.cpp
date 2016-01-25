#include "xml_utils.h"
#include <iostream>
#include "faust_exception.h"

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

xmlXPathContextPtr get_context(const char * filename)
{
	xmlDocPtr doc;
	 
	// Ouverture du document
	xmlKeepBlanksDefault(0); // Ignore les noeuds texte composant la mise en forme

	 doc = xmlParseFile(filename);
	// Initialisation de l'environnement XPath
	xmlXPathInit();
	// Création du contexte
	xmlXPathContextPtr ctxt = xmlXPathNewContext(doc); // doc est un xmlDocPtr représentant notre catalogue
	if (ctxt == NULL) {
		handleError("xml_utils",
			"get_context :Error while the creation of the xpath context ");
	}
}