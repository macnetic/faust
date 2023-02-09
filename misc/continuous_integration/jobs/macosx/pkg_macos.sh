#!/usr/bin/env bash

# Needed environment variables: FAUST_VERSION, MACOS_PY_VER, DURL, DFILE, EXPERIMENTAL_PKG, MACOS_PKG_STORE_PATH (optional)

export PYTHON_PATH=$(which python$MACOS_PY_VER)
if [[ ! -d 'build' ]]
then
	mkdir build
fi
cd build
cmake -DOpenMP_gomp_LIBRARY=/opt/local/lib/libomp/libgomp.dylib -DBUILD_WRAPPER_PYTHON=ON -DBUILD_DOCUMENTATION=ON -DCMAKE_INSTALL_PREFIX=/opt/local/faust -DCPACK_PACKAGE_VERSION=$FAUST_VERSION -DCMAKE_BUILD_TYPE=Release -DEXCLUDE_FAUST_LIB_INSTALL=ON -DBUILD_TESTING=OFF -DREMOTE_DATA_URL="$DURL" -DREMOTE_DATA_FILE="$DFILE" -DBUILD_MULTITHREAD=ON -DNOPY2=ON -DCMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-8.0 -DBUILD_FLOAT_PYX=ON -DBUILD_FLOAT_MEX=ON -DEXPERIMENTAL_PKG=$EXPERIMENTAL_PKG ..
make LANG=en_GB.UTF-8
# to compile matlab wrappers use matlab libiomp5 but backup the clang lib first
#- sudo mv /opt/local/lib/libomp/libomp.dylib /opt/local/lib/libomp/libomp.dylib.bak
sudo cp /opt/local/lib/libomp/libiomp5_matlab.dylib /opt/local/lib/libomp/libomp.dylib
cmake -DBUILD_WRAPPER_MATLAB=ON ..
make LANG=en_GB.UTF-8
# restore clang libomp
#- sudo cp /opt/local/lib/libomp/libomp.dylib /opt/local/lib/libomp/libiomp5.dylib
sudo cp /opt/local/lib/libomp/libomp_macports.dylib /opt/local/lib/libomp/libomp.dylib
# ensure the linking of libomp to the python wrapper shared lib is right (not overridden by matlab libomp)
otool -L wrapper/python/_FaustCorePy.cpython-*-darwin.so
for f in wrapper/python/_FaustCorePy.cpython-*-darwin.so; do sudo install_name_tool -change @rpath/libiomp5.dylib /opt/local/lib/libomp/libomp.dylib $f;done
otool -L wrapper/python/_FaustCorePy.cpython-*-darwin.so
sudo rm -Rf /opt/local/faust # install dir
sudo make install LANG=en_GB.UTF-8
# ensure libomp path also in install path
#TODO: why doing it two times?
for f in wrapper/python/_FaustCorePy.cpython-*-darwin.so;do sudo install_name_tool -change @rpath/libiomp5.dylib /opt/local/lib/libomp/libomp.dylib $f;done
otool -L wrapper/python/_FaustCorePy.cpython-*-darwin.so
#'sudo hdiutil create -volname Faust-$FAUST_VERSION-MatlabR2016a-Py2.7 -srcfolder /opt/local/faust -ov -format UDRW faust-$FAUST_VERSION'
sudo pkgbuild --identifier fr.inria.faust --version $FAUST_VERSION --root /opt/local/faust --install-location /opt/local/faust --scripts . ./faust-$FAUST_VERSION.pkg
for FILE in $(find /usr/local/lib ! -iname "libmatio*" -maxdepth 1 -mindepth 1); do filter_list+="--filter $(basename $FILE) "; done;
sudo pkgbuild --identifier fr.inria.faust.matio --version 1.5.12 --root /usr/local/lib $filter_list --install-location /usr/local/lib ./matio-bin-1.5.12.pkg
sudo pkgbuild --identifier fr.inria.faust.openmp --version 367070 --root /opt/local/lib/libomp --install-location /opt/local/lib/libomp ./libomp-367070.pkg
productbuild --synthesize --package ./matio-bin-1.5.12.pkg --package libomp-367070.pkg --package faust-$FAUST_VERSION.pkg ./distribution.plist
sed -e 's/\(.*pkg-ref id=.fr.inria.faust".*\)/\1<title>FAµST '$FAUST_VERSION'<\/title><license file="licenses.html"\/><readme file="installer_readme.html"\/>/' distribution.plist > tmp.plist; mv tmp.plist distribution.plist
productbuild --distribution ./distribution.plist --package-path matio-bin-1.5.12.pkg --package-path libomp-367070.pkg --package-path faust-$FAUST_VERSION.pkg --resources doc ./faust-matio-omp-$FAUST_VERSION.pkg
mv -f ./faust-matio-omp-$FAUST_VERSION.pkg ./faust-$FAUST_VERSION.pkg
#- if [[ -d $MACOS_PKG_STORE_PATH ]]; then sudo cp faust-$FAUST_VERSION.dmg faust-$FAUST_VERSION.pkg $MACOS_PKG_STORE_PATH; fi
if [[ -d "$MACOS_PKG_STORE_PATH" ]]; then sudo cp faust-$FAUST_VERSION.pkg $MACOS_PKG_STORE_PATH; fi

