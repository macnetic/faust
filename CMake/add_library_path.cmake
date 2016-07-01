function(add_library_path LIB_TMP_PATH)
	if(${ARGC} LESS 2)
		message(FATAL_ERROR "ARGC=${ARGC} : Wrong number of arguments in add_library_path macro")
	endif()
	string ( REGEX REPLACE "\\\\" "/" ARGN2 "${ARGN}")
	#set(list_var "${ARGN}")
	#foreach(loop_var IN LISTS list_var)
	foreach(loop_var IN LISTS ARGN2)
		if(EXISTS "${loop_var}")
			if(EXISTS "${loop_var}/lib64")
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var}/lib64 )
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} PARENT_SCOPE )
			elseif(EXISTS "${loop_var}/lib")
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var}/lib )
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} PARENT_SCOPE )
			else()
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var}  )
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} PARENT_SCOPE )
			endif()
		else()
			#message(STATUS "${loop_var} directory has been ignored, as it does not exist")
		endif()
	endforeach()
endfunction()


function(add_include_path INC_TMP_PATH)
	if(${ARGC} LESS 2)
		message(FATAL_ERROR "ARGC=${ARGC} : Wrong number of arguments in add_include_path macro")
	endif()
	string ( REGEX REPLACE "\\\\" "/" ARGN2 "${ARGN}")
	#set(list_var "${ARGN}")
	#foreach(loop_var IN LISTS list_var)
	foreach(loop_var IN LISTS ARGN2)
		if(EXISTS "${loop_var}")
			if(EXISTS "${loop_var}/include")
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} "${loop_var}/include")
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} PARENT_SCOPE )
			elseif(EXISTS "${loop_var}/inc")
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} "${loop_var}/inc" )
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} PARENT_SCOPE )
			else()
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} ${loop_var} )
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} PARENT_SCOPE )
			endif()
		else()
			#message(STATUS "${loop_var} directory has been ignored, as it does not exist")
		endif()
	endforeach()
endfunction()





