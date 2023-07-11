if(NOT EXISTS faust_misc_data.zip)
	message(STATUS "Downloading FAµST misc data archive (needed for tests).")
	file(DOWNLOAD https://gitlab.inria.fr/faustgrp/gforge_files/-/raw/master/faust_misc_data-2023-07-11.zip faust_misc_data.zip SHOW_PROGRESS)
endif()
if(NOT EXISTS data)
	file(ARCHIVE_EXTRACT INPUT faust_misc_data.zip DESTINATION ${FAUST_MISC_DIR})
endif()

