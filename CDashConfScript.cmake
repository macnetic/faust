cmake_minimum_required(VERSION 3.0.2)


set (CTEST_SOURCE_DIRECTORY "${CTEST_SCRIPT_DIRECTORY}")

set(CTEST_BUILD_NAME "${CMAKE_SYSTEM}_${CMAKE_HOST_SYSTEM_PROCESSOR}")


if(WIN32)
	set(CTEST_SITE "FaustWin")
	set (CTEST_CMAKE_GENERATOR "MinGW Makefiles")
	# set above seems to be ignored
	set(CONF_OPTIONS "${CONF_OPTIONS} -G 'MinGW Makefiles' -DMATIO_INC_DIR=$ENV{MATIO_INC_DIR} -DMATIO_LIB_FILE=$ENV{MATIO_LIB_FILE}")
elseif(APPLE AND UNIX)
	set (CTEST_CMAKE_GENERATOR "Unix Makefiles")
	set(CTEST_SITE "FaustMacOS")
elseif(UNIX)
	set (CTEST_CMAKE_GENERATOR "Unix Makefiles")
	set(CTEST_SITE "FaustLinux")
else()
	message(FATAL_ERROR "Unknown system.")
endif()

if($ENV{NOCPPTESTS} MATCHES "ON")
	set(CONF_OPTIONS "${CONF_OPTIONS} -DNOCPPTESTS=ON")
endif()

if($ENV{BUILD_WRAPPER_PYTHON} MATCHES "ON")
	set(CTEST_SITE "${CTEST_SITE}Python")
	#set(BUILD_WRAPPER_PYTHON ON CACHE BOOL "" FORCE) #ignored by configure
	set(CONF_OPTIONS "${CONF_OPTIONS} -DBUILD_WRAPPER_PYTHON=ON -DBUILD_WRAPPER_MATLAB=OFF")
endif()
if($ENV{BUILD_WRAPPER_MATLAB} MATCHES "ON")
	set(CTEST_SITE "${CTEST_SITE}Matlab")
	#set(BUILD_WRAPPER_MATLAB ON CACHE BOOL "" FORCE)
	set(CONF_OPTIONS "${CONF_OPTIONS} -DBUILD_WRAPPER_MATLAB=ON -DBUILD_WRAPPER_PYTHON=OFF")
endif()

if($ENV{BUILD_MULTITHREAD} MATCHES "ON")
	set(CONF_OPTIONS "${CONF_OPTIONS} -DBUILD_MULTITHREAD=ON")
endif()

if($ENV{FAUST_TORCH} MATCHES "ON")
	set(CONF_OPTIONS "${CONF_OPTIONS} -DFAUST_TORCH=ON")
endif()

if($ENV{NOPY2} MATCHES "ON")
	set(CONF_OPTIONS "${CONF_OPTIONS} -DNOPY2=ON")
endif()

if(DEFINED ENV{CMAKE_PREFIX_PATH})
	set(CONF_OPTIONS "${CONF_OPTIONS} -DCMAKE_PREFIX_PATH=$ENV{CMAKE_PREFIX_PATH}")
endif()

if(DEFINED ENV{OpenMP_gomp_LIBRARY})
	set(CONF_OPTIONS "-DOpenMP_gomp_LIBRARY=$ENV{OpenMP_gomp_LIBRARY} ${CONF_OPTIONS}")
endif()

if(DEFINED ENV{OpenMP_INC_DIR})
	set(CONF_OPTIONS "-DOpenMP_INC_DIR=$ENV{OpenMP_INC_DIR} ${CONF_OPTIONS}")
endif()

set(CTEST_BINARY_DIRECTORY "build_${CTEST_SITE}")
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

# https://docs.gitlab.com/ee/ci/variables/
message(STATUS "The git branch is:" $ENV{CI_COMMIT_REF_NAME})
message(STATUS "The git commit is:" $ENV{CI_COMMIT_SHA})

if($ENV{SLOW_TESTS} MATCHES "OFF")
	set(CONF_OPTIONS "${CONF_OPTIONS} -DSLOW_TESTS=OFF")
else()
	set(CONF_OPTIONS "${CONF_OPTIONS} -DSLOW_TESTS=ON")
endif()

if(UNIX)
	set(CONF_OPTIONS "${CONF_OPTIONS} -DCMAKE_INSTALL_PREFIX=$ENV{HOME}")
else() #WIN32
	set(CONF_OPTIONS "${CONF_OPTIONS} -DCMAKE_INSTALL_PREFIX=C:/Users/$ENV{USERNAME}")
endif()

set(CONF_OPTIONS "${CONF_OPTIONS} -DEXPERIMENTAL_PKG=ON")

#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY}) # no need to empty build dir. because
# we use CTEST_START_WITH_EMPTY_BINARY_DIRECTORY above and gitlab-runner starts with a new one

CTEST_START("Experimental") # because we don't update the code (gitlab-runner does it for us)
message(STATUS "The site name is: " ${CTEST_SITE})
#CTEST_START("Nightly")
#CTEST_START("Continuous")
#CTEST_UPDATE() # consistently with experimental mode

set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} ${CONF_OPTIONS} ${CTEST_SOURCE_DIRECTORY}") # cmake is anyway the default configure command
message(STATUS "CONFIGURE COMMAND: ${CTEST_CONFIGURE_COMMAND}")
CTEST_CONFIGURE() #OPTIONS ${CONF_OPTIONS} doesn't work (even with a list()) so we set the ctest_configure_command above
# no OPTIONS (arg)
CTEST_BUILD(TARGET install) #need to install for python tests (quickstart.py)
#CTEST_BUILD()

#IF(UNIX)
#	set(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:${CTEST_INSTALL_DIR}/lib")
#ENDIF(UNIX)
# no shared libraries to look at


CTEST_TEST()
if (WITH_COVERAGE AND CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif (WITH_COVERAGE AND CTEST_COVERAGE_COMMAND)
if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
CTEST_SUBMIT()
