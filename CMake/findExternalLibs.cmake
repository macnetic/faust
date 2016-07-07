###### find external libraries ######
if(WIN32)
	set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} "lib" )
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} ".a")
endif(WIN32)


# functions tools is looking for library and include files in the environment variables and/or externals Path using the cmake find command.
include(CMake/add_library_path.cmake)
include(CMake/check_external_libraries.cmake)


#set(OPENBLAS_DEFAULT_DIR=/opt/OpenBLAS CACHE INTERNAL "")
#set(EIGEN_DEFAULT_DIR=/usr/include/eigen3 CACHE INTERNAL "")

#set(OPENBLAS_EXTERNALS_DIR=${PROJECT_SOURCE_DIR}/externals/OpenBLAS CACHE INTERNAL "")
#set(EIGEN_EXTERNALS_DIR=${PROJECT_SOURCE_DIR}/externals/eigen CACHE INTERNAL "") 

#message(STATUS "DEFAULT EIGEN_DEFAULT_DIR=${EIGEN_DEFAULT_DIR}")
#message(STATUS "DEFAULT EIGEN_EXTERNALS_DIR=${EIGEN_EXTERNALS_DIR}")
# /sw /opt/local
#add_library_path(LIBRARY_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{HDF5_ROOT_DIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/usr/local/lib" "/opt" "/opt/local" "/usr/lib/x86_64-linux-gnu/" "${PROJECT_SOURCE_DIR}/externals" )
#add_include_path(INCLUDE_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{MATIODIRINC}" "$ENV{OPENBLASDIR}" "$ENV{EIGENDIR}" "/usr" "/usr/local" "/usr/include/libxml2" "/opt"  "/opt/local" "${PROJECT_SOURCE_DIR}/externals" )

#add_library_path(LIBRARY_PATH_LIST_TMP "${OPENBLAS_DEFAULT_DIR}" "${OPENBLAS_EXTERNALS_DIR}")
#add_include_path(INCLUDE_PATH_LIST_TMP "${OPENBLAS_DEFAULT_DIR}" "${OPENBLAS_EXTERNALS_DIR}" "${EIGEN_DEFAULT_DIR}" "${EIGEN_EXTERNALS_DIR}")


message(STATUS "******* Check externals library ***********")

# Default path library (where the library is automatically install)
add_library_path(LIBRARY_PATH_LIST_TMP_DEFAULT 	"/opt/OpenBLAS"
												"${PROJECT_SOURCE_DIR}/externals/unix/OpenBLAS"
												"/usr/local" #pour matio et xml2
												"/usr/local/lib" #pour matio et xml2 
												"/usr/local/lib64" #pour matio et xml2
												"${PROJECT_SOURCE_DIR}/externals/unix/matio/src/.libs" #pour matio
												"${PROJECT_SOURCE_DIR}/externals/unix/libxml2/.libs" # pour libxml
)

add_include_path(INCLUDE_PATH_LIST_TMP_DEFAULT 	"/opt/OpenBLAS" 
												"${PROJECT_SOURCE_DIR}/externals/unix/OpenBLAS"
												"/usr/include/eigen3"
												"${PROJECT_SOURCE_DIR}/externals/unix/eigen"
												"/usr/local" #pour matio et xml2 
												"/usr/include/libxml2"
												"${PROJECT_SOURCE_DIR}/externals/unix/matio/src" #pour matio
												"${PROJECT_SOURCE_DIR}/externals/unix/libxml2/include" # pour libxml
)

#set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_DEFAULT}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
#set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_DEFAULT})# CACHE PATH "List of include paths used as PATH parameter in find_path")

#message(STATUS "DEFAULT LIBRARY_PATH_LIST=${LIBRARY_PATH_LIST}")
#message(STATUS "DEFAULT INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")

######## EIGEN ##############
include(CMake/findEIGENLib.cmake)
######## OPENBLAS ##############
if (BUILD_OPENBLAS)
	include(CMake/findOPENBLASLib.cmake)
endif(BUILD_OPENBLAS)
######## MATLAB INCLUDE AND LIBRARY ##################
if (BUILD_MATLAB_MEX_FILES)
	include(CMake/findMatlab.cmake)
endif(BUILD_MATLAB_MEX_FILES)
######## MATIO INCLUDE AND LIBRARY ##################
if (BUILD_READ_MAT_FILE)
	include(CMake/findMATIOLib.cmake)
endif(BUILD_READ_MAT_FILE)
######## XML INCLUDE AND LIBRARY ##################
if (BUILD_READ_XML_FILE)
	include(CMake/findXML2Lib.cmake)
endif(BUILD_READ_XML_FILE)



#add_library_path(LIBRARY_PATH_LIST_TMP3 "$ENV{CUDADIR}" "$ENV{HDF5_ROOT_DIR}" "/usr/lib/x86_64-linux-gnu/")
#add_include_path(INCLUDE_PATH_LIST_TMP3 "$ENV{CUDADIR}" "/usr/include/libxml2")

add_library_path(LIBRARY_PATH_LIST_TMP3 "$ENV{CUDADIR}" )
add_include_path(INCLUDE_PATH_LIST_TMP3 "$ENV{CUDADIR}" )


set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP3}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP3}) # CACHE PATH "List of include paths used as PATH parameter in find_path")


# trop gros (100 mb le tar...) ne sert pas à grand chose dans le projet. 
#check_external_libraries(hdf5 HDF5_LIB_FILE 0)

#LDFLAGS = -L$(CUDADIR)/lib64 -L$(MATIODIR)/lib -lpthread -lm -lcublas -lcudart -lcusparse -lstdc++ -lgfortran -lz -lmatio -lhdf5




if (BUILD_USE_GPU)
	check_external_libraries(cublas CUBLAS_LIB_FILE 1)
	check_external_libraries(cudart CUDART_LIB_FILE 1)
	check_external_libraries(cusparse CUSPARSE_LIB_FILE 1)

	check_external_includes("cublas_v2.h" CUBLAS_V2_INC_DIR 1)
	check_external_includes("cusparse.h" CUSPARSE_INC_DIR 1)
	check_external_includes("cuda.h" CUDA_INC_DIR 1)
	check_external_includes("cuda_runtime_api.h" CUDA_RUNTIME_API_INC_DIR 1)
endif(BUILD_USE_GPU)



