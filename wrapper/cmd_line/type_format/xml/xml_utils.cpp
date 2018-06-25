/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
	return ctxt;
}
