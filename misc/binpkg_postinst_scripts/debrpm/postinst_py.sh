#!/bin/bash

uid=$(id -u)
[[ ! "$uid" = 0 ]] && echo "Error: need to be root." >&2 && exit 1 # not needed, we're supposed to be root at installation time

[[ -n "$1" && -d "$1" ]] && FAUST_PY_WRAPPER_PATH="$1" || FAUST_PY_WRAPPER_PATH="@CMAKE_INSTALL_PYTHON_PREFIX@"

[[ -z "$FAUST_PY_WRAPPER_PATH" ]] && echo "USAGE: $0 <path>" && exit 2
[[ ! -d "$FAUST_PY_WRAPPER_PATH" ]] && echo "ERROR: directory $FAUST_PY_WRAPPER_PATH doesn't exist" && exit 3

SUPPORTED_PY3=@PY3_MINOR_VER@ # needed because python3 has ABI changes between minor versions (the pyfaust shared lib. is compiled specifically for python major and minor version)

function link_py_files(){
	PY_MAJOR_VER=$1
	PY_MAJOR_MINOR_VER=$(python$PY_MAJOR_VER --version 2>&1 | sed -e 's/.*[[:blank:]]\{1,\}\([[:digit:]]\.[[:digit:]]\)\.[[:digit:]]\{1,\}/\1/')
	# on fedora python package path are suffixed by site-packages, on ubuntu it's rather dist-packages (TODO: on centos ? on debian ?)
	PY_SITE_PACKAGES_PATH=$(python$PY_MAJOR_VER -c 'import os,site,re;print([path for path in site.getsitepackages() if (re.match(".*"+os.path.sep+"site-packages$", path) or re.match(".*"+os.path.sep+"dist-packages$", path)) and os.path.exists(path) ][0])')
	[[ ! -d ${PY_SITE_PACKAGES_PATH} ]] && echo -e "\033[1mWARNING\033[0m: couldn't link the Faust py wrapper for python$PY_MAJOR_VER to your system. Either python$PY_MAJOR_VER is not installed on your system and that's pretty normal or it failed for another reason and you'll have to add Faust to your PYTHONPATH manually for using it (see the documentation)." >&2 && return 1
	PYFILES=${FAUST_PY_WRAPPER_PATH}/pyfaust
	[[ -n "$DEBUG" ]] && echo PY_MAJOR_VER=$PY_MAJOR_MINOR_VER
	if [[ "$PY_MAJOR_VER" = 3* ]]
	then
		[[ ! "$PY_MAJOR_MINOR_VER" = 3.$SUPPORTED_PY3 ]] && echo -e "\033[1mWARNING\033[0m: your python3 version ($PY_MAJOR_MINOR_VER) is not supported by Faust (only 3.$SUPPORTED_PY3 is supported). Please rely on Faust for python2." && return 2 || PYFILES+=" ${FAUST_PY_WRAPPER_PATH}/_FaustCorePy.cpython-3${SUPPORTED_PY3}m-x86_64-linux-gnu.so"
	else #py2
		PYFILES+=" ${FAUST_PY_WRAPPER_PATH}/_FaustCorePy.so"
	fi
	[[ -n "$DEBUG" ]] && echo PYFILES=$PYFILES
	local PYLINKS_FAIL_COUNT
	typeset -i PYLINKS_FAIL_COUNT=0
	for PYFILE in $PYFILES
	do
		PYLINK="${PY_SITE_PACKAGES_PATH}/$(basename "$PYFILE")"
		[[ -n "$DEBUG" ]] && echo ln -sf "$PYFILE" "${PY_SITE_PACKAGES_PATH}/"
		[[ -n "$DEBUG" ]] && echo PYLINK=$PYLINK
		[[ -r "${PYLINK}" ]] && rm -f "$PYLINK"
		[[ -r "$PYFILE" ]] && ln -sf "$PYFILE" "${PY_SITE_PACKAGES_PATH}/"
		[[ ! -r "$PYLINK" ]] && PYLINKS_FAIL_COUNT+=1
		[[ -n "$DEBUG" ]] && echo  PYLINKS_FAIL_COUNT=$PYLINKS_FAIL_COUNT
	done
	[[ $PYLINKS_FAIL_COUNT = 0 ]] && echo -e "\033[1mNOTICE\033[0m: Linked py wrapper for python$PY_MAJOR_VER into ${PY_SITE_PACKAGES_PATH}. Installation's ok for python$PY_MAJOR_MINOR_VER."
}

for V in 2 3
do
	which python$V 2>&1 > /dev/null && link_py_files $V
done


