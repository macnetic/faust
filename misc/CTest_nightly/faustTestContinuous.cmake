# DO NOT MODIFIED !!!
# This script launch the nightly faust project from the "run_nightly_XXX.sh"
# All the parameters, config optional, path of library, environment variables are defined in "run_nightly_XXX.sh"


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4)

IF(WIN32)
	SET(CTEST_CMAKE_GENERATOR "Visual Studio 9 2008")
	SET (CTEST_SOURCE_DIRECTORY "$ENV{USERPROFILE}/faust_test")
	SET (CTEST_BINARY_DIRECTORY "$ENV{USERPROFILE}/faust_test_build")
	SET (CTEST_INSTALL_DIR      "$ENV{USERPROFILE}/faust_output")
ELSEIF(UNIX)
	SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
	SET (CTEST_SOURCE_DIRECTORY "$ENV{PATH_DIR_RUN_NIGHTLY}/faust_test")
	SET (CTEST_BINARY_DIRECTORY "$ENV{PATH_DIR_RUN_NIGHTLY}/faust_test_build")
	SET (CTEST_INSTALL_DIR      "$ENV{PATH_DIR_RUN_NIGHTLY}/faust_output")	
ENDIF(WIN32)


FIND_PROGRAM(CTEST_SVN_COMMAND NAMES svn)
SET (CTEST_SVN_CHECKOUT  "${CTEST_SVN_COMMAND} checkout --username testcdash --password testcdash https://scm.gforge.inria.fr/svn/faust/trunk/devcpp \"${CTEST_SOURCE_DIRECTORY}\"")
SET (CTEST_CHECKOUT_COMMAND "${CTEST_SVN_CHECKOUT}")


site_name(CTEST_SITE)
set(CTEST_BUILD_NAME "${CMAKE_SYSTEM}_${CMAKE_HOST_SYSTEM_PROCESSOR}")

# met les variables d'environnement a jour pour trouver les differentes librairies
#SET(ENV{PATH} "$ENV{PATH}:/usr/local/MATLAB/R2014b/bin")
#SET(ENV{OPENBLASDIR} "/home/ci/local/OPENBLAS/")
#SET(ENV{HDF5_ROOT_DIR} "/home/ci/local/HDF5")
SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

CTEST_START("Continuous")
CTEST_UPDATE()
#IF(res EQUAL -1) 
#	RETURN()
#ENDIF(res EQUAL -1)
###  run build and test for continuous integration model only if there was a new update
#IF(res EQUAL 0)
#	RETURN()
#ENDIF(res EQUAL 0)


# Option of the ctest build, set from run_nightly_XXX.sh 
set (CTEST_CONFIGURE_OPTIONS_TMP "$ENV{CTEST_CONFIGURE_OPTIONS}" )
set (CTEST_CONFIGURE_OPTIONS_OK "-DCMAKE_INSTALL_PREFIX=${CTEST_INSTALL_DIR};${CTEST_CONFIGURE_OPTIONS_TMP}")
message(STATUS "CTEST_CONFIGURE_OPTIONS_OK=${CTEST_CONFIGURE_OPTIONS_OK}")
#message(STATUS " CTEST_CONFIGURE\(OPTIONS \"${CTEST_CONFIGURE_OPTIONS_OK}\"\) ")
CTEST_CONFIGURE(OPTIONS "${CTEST_CONFIGURE_OPTIONS_OK}")
#CTEST_CONFIGURE(OPTIONS "-DBUILD_DEBUG:BOOL=OFF;-DDASH_TESTING:BOOL=ON;-DBUILD_COVERAGE:BOOL=ON;-DFAUST_GEN_DOC:BOOL=OFF;-DCMAKE_INSTALL_PREFIX=${CTEST_INSTALL_DIR};-DFAUST_USE_GPU:BOOL=ON;-DFAUST_USE_OPENBLAS=ON")

#CTEST_BUILD()
CTEST_BUILD( TARGET install)

IF(UNIX)
	SET(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:${CTEST_INSTALL_DIR}/lib")
ENDIF(UNIX)

CTEST_TEST()
if (WITH_COVERAGE AND CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif (WITH_COVERAGE AND CTEST_COVERAGE_COMMAND)
if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
CTEST_SUBMIT()
