#!/bin/bash

uid=$(id -u)
[[ ! "$uid" = 0 ]] && echo "Error: need to be root." >&2 && exit 1 # not needed, we're supposed to be root at installation time

[[ -n "$1" && -d "$1" ]] && FAUST_PY_WRAPPER_PATH="$1" || FAUST_PY_WRAPPER_PATH="@CMAKE_INSTALL_PYTHON_PREFIX@"

[[ -z "$FAUST_PY_WRAPPER_PATH" ]] && echo "USAGE: $0 <path>" && exit 2
[[ ! -d "$FAUST_PY_WRAPPER_PATH" ]] && echo "ERROR: directory $FAUST_PY_WRAPPER_PATH doesn't exist" && exit 3

function link_py_files(){
	PY_MAJOR_VER=$1
	PY_MAJOR_MINOR_VER=$(python$PY_MAJOR_VER --version 2>&1 | sed -e 's/.*[[:blank:]]\{1,\}\([[:digit:]]\.[[:digit:]]\)\.[[:digit:]]\{1,\}/\1/')
	PY_SITE_PACKAGES_PATH=$(python$PY_MAJOR_VER -c 'import os,site,re;print([path for path in site.getsitepackages() if re.match(".*"+os.path.sep+"site-packages$", path) and os.path.exists(path) ][0])')
	PYFILES=${FAUST_PY_WRAPPER_PATH}/FaustPy.py
	[[ "$PY_MAJOR_MINOR_VER" = 3* ]] && PYFILES+=" ${FAUST_PY_WRAPPER_PATH}/FaustCorePy.cpython-3*so" || \
					    PYFILES+=" ${FAUST_PY_WRAPPER_PATH}/FaustCorePy.so"
	for PYFILE in $PYFILES
	do
		[[ -n "$DEBUG" ]] && echo ln -sf "$PYFILE" "${PY_SITE_PACKAGES_PATH}/"
		[[ -r "$PYFILE" ]] && ln -sf "$PYFILE" "${PY_SITE_PACKAGES_PATH}/"
	done && echo "Linked py wrapper version $PY_MAJOR_VER into ${PY_SITE_PACKAGES_PATH}."

}

for V in 2 3
do
	which python$V 2>&1 > /dev/null && link_py_files $V
done

exit 0 # to avoid exit failure result when only one version of py is found (or it'll make the macOS installer fail)

