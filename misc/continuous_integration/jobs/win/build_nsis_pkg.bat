:: "Needed env. variables: FAUST_VERSION, EXPERIMENTAL_PKG, BUILD_CONFIG, DURL, DFILE, WIN_PY_VER
if NOT EXIST build (mkdir build) else (rmdir /S /Q build & mkdir build)
cd build
:: build first only matfaust
cmake -G "Visual Studio 16 2019" -DBUILD_WRAPPER_MATLAB=ON -DBUILD_WRAPPER_PYTHON=OFF -DSLOW_TESTS=OFF -DCPACK_PACKAGE_VERSION=%FAUST_VERSION% -DBUILD_DOCUMENTATION=ON -DEXCLUDE_FAUST_LIB_INSTALL=ON -DCMAKE_INSTALL_PREFIX=win_pkg_build -DBUILD_TESTING=OFF -DMATIO_LIB_FILE=C:/faust_libs/libmatio_standalone.lib -DAPI_DOC_BASE_URL="file:///C:/Program Files/Faust/doc/" -DREMOTE_DATA_URL="%DURL%" -DREMOTE_DATA_FILE="%DFILE%" -DEXPERIMENTAL_PKG=%EXPERIMENTAL_PKG% -DUSE_GPU_MOD=ON -DCMAKE_PREFIX_PATH=../gpu_mod -DPYTHON_ENCODING=windows-1252 -DBUILD_MULTITHREAD=ON -DBUILD_FLOAT_PYX=ON -DBUILD_FLOAT_MEX=ON -DCMAKE_BUILD_TYPE=%BUILD_CONFIG% -DVCOMPLIB_PATH=C:\faust_libs\vcomp140.dll ..
cmake --build . --config %BUILD_CONFIG%
:: now build pyfaust too
:: specifying a consistent python version with pkg_win_purepy_rev
cmake -DBUILD_WRAPPER_PYTHON=ON -DMATIO_LIB_FILE=C:/faust_libs/libmatio_standalone.lib ..
cd wrapper\python
:: it costs no time of building if pkg_win_purepy_* already compiled the .pyd lib
py -%WIN_PY_VER% setup.py build_ext --inplace
cd ..\..
makensis faust.nsi

