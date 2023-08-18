##############################################################################
##                              Description:                                ##
##  Cmake file to configure the test that will be send to CDash plateform   ##
##									    ##                ##                                                                          ##							
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.inria.fr>                         ##
##                                                                          ##
##                              License:                                    ##
##  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
##                      Luc Le Magoarou, Remi Gribonval                     ##
##                      INRIA Rennes, FRANCE                                ##
##                      http://www.inria.fr/                                ##
##                                                                          ##
##  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
##  General Public License.                                                 ##
##  This program is free software: you can redistribute it and/or modify    ##
##  it under the terms of the GNU Affero General Public License as          ##
##  published by the Free Software Foundation.                              ##
##                                                                          ##
##  This program is distributed in the hope that it will be useful, but     ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
##  See the GNU Affero General Public License for more details.             ##
##                                                                          ##
##  You should have received a copy of the GNU Affero General Public        ##
##  License along with this program.                                        ##
##  If not, see <http://www.gnu.org/licenses/>.                             ##
##                                                                          ##
##                             Contacts:                                    ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##      Hakim Hadj-Djilani: hakim.hadj-djilani@inria.fr                     ##
##      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
##      Adrien Leman    : adrien.leman@inria.fr                             ##
##      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##############################################################################


set(CTEST_PROJECT_NAME "faust")
set(CTEST_NIGHTLY_START_TIME "18:30:00 CET")
#set(CTEST_DROP_METHOD "http")
#set(CTEST_DROP_SITE "cdash.inria.fr")
#set(CTEST_DROP_LOCATION "/CDash/submit.php?project=faust")
set(CTEST_DROP_SITE "cdash-ci.inria.fr")
set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_LOCATION "/submit.php?project=faust")
set(CTEST_DROP_SITE_CDASH TRUE)


set(CTEST_FULL_OUTPUT TRUE)
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE TRUE)
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE TRUE)
set(CTEST_CURL_OPTIONS "CURLOPT_SSL_VERIFYPEER_OFF")

if(DEFINED ENV{CI_COMMIT_SHA}) # it doesn't work without DEFINED for env. var.
	set(CTEST_BUILD_NAME ${CTEST_BUILD_NAME}-$ENV{CI_COMMIT_SHA})
endif()
