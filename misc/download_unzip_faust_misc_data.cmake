if(NOT EXISTS faust_misc_data.zip)
	message(STATUS "Downloading FAµST misc data archive (needed for tests).")
	file(DOWNLOAD https://zenodo.org/records/10613337/files/faust_misc_data-2024-02-03.zip faust_misc_data.zip SHOW_PROGRESS)
endif()
if(NOT EXISTS data)
	if(NOT FAUST_MISC_DIR)
		set(FAUST_MISC_DIR ./misc)
	endif()
	file(ARCHIVE_EXTRACT INPUT faust_misc_data.zip DESTINATION ${FAUST_MISC_DIR})
endif()

