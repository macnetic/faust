file(GLOB MATLAB_SCRIPTS "*.m")
file(GLOB MATLAB_SCRIPTS "*.m")

# configure_file(${FAUST_MEXINTERFACETOOLS_SRC_DIR}/matlab_faust.m.in ${FAUST_TESTMATLABSCRIPT_BIN_DIR}/matlab_faust.m COPYONLY)





if(FAUST_USE_MEX)
	foreach(matlab_script localize_src.m figure_localize_src.m)
      configure_file(${FAUST_DEMO_SRC_DIR}/${matlab_script} ${FAUST_DEMO_BIN_DIR}/${matlab_script} COPYONLY)
   endforeach()
	configure_file(${FAUST_DEMO_SRC_DIR}/set_path.m.in ${FAUST_DEMO_BIN_DIR}/set_path.m @ONLY)

endif()


