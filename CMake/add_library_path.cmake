macro(add_library_path LIB_TMP_PATH)
	if(${ARGC} LESS 2)
		message(FATAL_ERROR "ARGC=${ARGC} : Wrong number of arguments in add_library_path macro")
	endif()
	set(list_var "${ARGN}")
	foreach(loop_var IN LISTS list_var)
		if(EXISTS "${loop_var}")
			if(EXISTS "${loop_var}/lib64")
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var}/lib64 )
			elseif(EXISTS "${loop_var}/lib")
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var}/lib )
			else()
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var} )
			endif()
		else()
			message(STATUS "${loop_var} directory has been ignored, as it does not exist")
		endif()
	endforeach()
endmacro()






