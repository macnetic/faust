#!/bin/bash
#set -e # won't continue any further if an error occurs

# This script pupose is to build a conda pkg from wheel online pip package
# for linux, macosx or windows (by running the script on the same system) and upload it

SYSTEM_ERROR="system must be either: linux, macosx or win"
PYVER_ERROR="python version must be in the format major_version.minor_version (e.g. 3.9)"
META_YAML_PATH_ERROR="The meta.yaml.in filepath must exist and the filename must be meta.yaml.in"
PYPI_FILES_URL='https://pypi.org/project/pyfaust/#files'
TOKEN_FILE="$HOME/conda_pyfaust_token"

function usage
{
	echo "USAGE: $0 <system> <python_version> <path_to_meta.yaml.in>"
	echo $SYSTEM_ERROR
	echo $PYVER_ERROR
	echo "Note: the package must be available on pypi.org to be built and uploaded to conda repository."
}

function find_pkg
{
	local P OS
	if [[ "$SYSTEM" = macosx ]]
	then
		OS=osx
		P=$(dirname $CONDA_PATH_ENV_PREFIX)
	else
		OS=$SYSTEM
		P=$CONDA_PATH_ENV_PREFIX
	fi
	find $P -iname "*.tar.bz2" | grep conda-bld | grep pyfaust | grep $OS | grep $PYVER_NODOT | grep $PYFAUST_VERSION
}

###### COMMAND-LINE PARSING
if [ $# -lt 3 ]
then
	usage
	exit 1
fi


if [[ "$1" =~ linux|macosx|win ]]
then
	SYSTEM="$1"
else
	echo "Error: $SYSTEM_ERROR"
	echo
	usage
	exit 1
fi

if [[ "$2" =~ [[:digit:].]+ ]]
then
	PYVER="$2"
	PYVER_NODOT=$(echo $PYVER | tr -d '.')
else
	echo "Error: $PYVER_ERROR"
	echo
	usage
	exit 1
fi

if [[ ! "$3" = *meta.yaml.in || ! -r "$3" ]]
then
	echo "Error: $META_YAML_PATH_ERROR"
	echo
	usage
	exit 2
else
	META_YAML_PATH="$3"
fi

[[ -n "$DEBUG" ]] && echo system=$SYSTEM
[[ -n "$DEBUG" ]] && echo "python_version=$PYVER ($PYVER_NODOT)"
[[ -n "$DEBUG" ]] && echo meta.yaml.in path=$META_YAML_PATH

##### FIRST STEP: find the python wheel package URL on pypi.org

WHL_URL=$(curl  $PYPI_FILES_URL | sed -ne 's/.*"\([^"]\{1,\}\.whl\)".*/\1/p' | grep $SYSTEM | grep $PYVER_NODOT)

[[ -z "$WHL_URL" ]] && echo "No wheel package found on pypi.org for the system $SYSTEM_ERROR and the python version $PYVER" && exit 3

[[ -n "$DEBUG" ]] && echo WHEEL_URL=${WHL_URL}

PYFAUST_VERSION=$(echo $WHL_URL | sed -e 's/.*pyfaust-\([[:digit:].]\{1,\}\)-.*/\1/')
[[ -n "$DEBUG" ]] && echo PYFAUST_VERSION=$PYFAUST_VERSION

#### SECOND STEP: generate the meta.yaml file for conda pkg building using the meta.yaml.in model file
mkdir -p pyfaust
sed -e "s%@PYFAUST_VERSION@%$PYFAUST_VERSION%;s%@WHL_URL@%$WHL_URL%" $META_YAML_PATH > pyfaust/meta.yaml

#### THIRD STEP: build the package using conda

if ! which conda
then
	echo "Error: conda isn't in the PATH or is not installed."
	exit 4
fi

# set the environment
CONDA_ENV=build_upload_pyfaust_venv_$PYVER_NODOT
if ! conda env list | grep $CONDA_ENV
then
	[[ -n "$DEBUG" ]] && echo Conda env name to create: $CONDA_ENV
	conda create -y -n $CONDA_ENV python==$PYVER && conda run -n $CONDA_ENV conda install -y conda-build anaconda-client
	[[ ! "$?" = 0 ]] && echo "Error: failed to setup the conda environment" && exit 5
	#conda init bash
fi
# conda-activate is only for interactive shell, so use run
#conda activate $CONDA_ENV || echo "Error: failed to activate the conda env $ENV" && exit 6
CONDA_PATH_ENV_PREFIX=$(conda info | grep "envs directories" | head -1 | cut -d: -f2)
PKG=$(find_pkg)
if [[ -z "$PKG" || ! -r "$PKG" ]]
then
	conda run -n $CONDA_ENV conda build pyfaust
fi
[[ ! "$?" = 0 ]] && echo "Error: failed to build pyfaust conda package" && exit 7
# find the package
PKG=$(find_pkg)
#echo pyfaust | conda run -n $CONDA_ENV anaconda login
[[ -z "$PKG" || ! -r "$PKG" ]] && echo "Error: no built package was found in the virtual env directory." && exit 8
[[ ! -r $TOKEN_FILE ]] && echo "Error: the anaconda token file $TOKEN wasn't found." && exit 9
conda run -n $CONDA_ENV anaconda -t $TOKEN_FILE upload -u pyfaust $PKG
# removing the package: conda run -n $CONDA_ENV anaconda -t ~/conda_pyfaust_token remove
