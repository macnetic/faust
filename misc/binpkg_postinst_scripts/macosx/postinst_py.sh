#!/bin/bash

echo "$0 -- FAµST python wrapper post-install script START" | tee -a /tmp/log_faust_install

uid=$(id -u)
[[ ! "$uid" = 0 ]] && echo "Error: need to be root." >&2 | tee -a /tmp/log_faust_install >&2 && exit 1 # not needed, we're supposed to be root at installation time

[[ -n "$1" && -d "$1" ]] && FAUST_PY_WRAPPER_PATH="$1" || FAUST_PY_WRAPPER_PATH="@CMAKE_INSTALL_PYTHON_PREFIX@"

[[ -z "$FAUST_PY_WRAPPER_PATH" ]] && echo "USAGE: $0 <path>" | tee -a /tmp/log_faust_install >&2 && exit 2 # shouldn't happen with cmake
[[ ! -d "$FAUST_PY_WRAPPER_PATH" ]] && echo "ERROR: directory $FAUST_PY_WRAPPER_PATH doesn't exist" | tee -a /tmp/log_faust_install >&2 && exit 3 # shouldn't happen with cmake

function link_py_files(){
	PYTHON=$1
	DEBUG=1
	[[ -n "$DEBUG" ]] && echo PYTHON=$PYTHON | tee -a /tmp/log_faust_install
	[[ -d "$PYTHON" ]] && return
	PY_MAJOR_MINOR_VER=$($PYTHON --version 2>&1 | sed -e 's/.*[[:blank:]]\{1,\}\([[:digit:]]\.[[:digit:]]\)\.[[:digit:]]\{1,\}/\1/')
	PY_SITE_PACKAGES_PATH=$($PYTHON -c 'import os,site,re;print([path for path in site.getsitepackages() if re.match(".*"+os.path.sep+"site-packages$", path) and os.path.exists(path) ][0])')
	[[ ! -d "${PY_SITE_PACKAGES_PATH}" ]] && return
	[[ -n "$DEBUG" ]] && echo PY_SITE_PACKAGES_PATH=${PY_SITE_PACKAGES_PATH} | tee -a /tmp/log_faust_install
	PYFILES=${FAUST_PY_WRAPPER_PATH}/pyfaust
	[[ "$PY_MAJOR_MINOR_VER" = 3* ]] && PYFILES+=" ${FAUST_PY_WRAPPER_PATH}/_FaustCorePy.cpython-3*so" || \
					    PYFILES+=" ${FAUST_PY_WRAPPER_PATH}/_FaustCorePy.so"
	for PYFILE in $PYFILES
	do
		[[ -n "$DEBUG" ]] && echo ln -sf "$PYFILE" "${PY_SITE_PACKAGES_PATH}/" | tee -a /tmp/log_faust_install
		[[ -r "$PYFILE" ]] && ln -sf "$PYFILE" "${PY_SITE_PACKAGES_PATH}/"
	done && echo "Linked py wrapper version $PY_MAJOR_MINOR_VER into ${PY_SITE_PACKAGES_PATH}." | tee -a /tmp/log_faust_install
}


# python on macOS is always in /usr/local or /opt/local with macports
for DIR in /usr/local /opt/local
do
	for FILE in $(find $DIR -type l -name "python*") $(find $DIR -type f -name "python*")
	do
		echo "PYTHON found: $file" | tee -a /tmp/log_faust_install
		[[ -x "$FILE" ]] && link_py_files "$FILE"
	done
done

echo "$0 -- FAµST python wrapper post-install script END" | tee -a /tmp/log_faust_install
echo

