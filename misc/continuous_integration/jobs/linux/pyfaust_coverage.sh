#!/usr/bin/env bash

# This script outputs a test coverage report for pyfaust

#set -euo pipefail

[[ -z "$NUX_PY_VER" ]] && echo "I need NUX_PY_VER variable to be set" && exit 1

cmake -P $(dirname $0)/../../../../misc/download_unzip_faust_misc_data.cmake # tests need data


# use conda rather than venv because python built on linux VM has not been
# built with --enable-loadable-sqlite-extensions (and sqlite-devel)
# but it is a need of coverave python package
CONDA_ENV=test_coverage_pyfaust
conda create -y -n $CONDA_ENV python==$NUX_PY_VER
#conda activate test_coverage_pyfaust

shopt -s expand_aliases

# ease execution of command in $CONDA_ENV
alias CE="conda run -n $CONDA_ENV"

# install pyfaust with pip
if [[ $# -ge 1 ]]
then
    [[ ! "$1" =~ ^.*/?pyfaust.*whl$ ]] && echo "First arg. must be a filepath of a pyfaust whl package" >&2 && exit 1
    CE pip install $1
else
    CE pip install $(find ./ -name "*pyfaust*whl")
fi
CE conda install -c conda-forge -y coverage
CE coverage erase # just in case
PYFAUST_DIR=$(dirname $(CE python -c "import pyfaust as pf; print(pf.__file__)"))
CE coverage run --source $PYFAUST_DIR misc/test/src/Python/test_FaustPy.py
#coverage run -a --source $PYFAUST_DIR $PYFAUST_DIR/tests/run.py # only real and cpu, all tests are ran below
# take doctest into account
for MOD in factparams proj poly tools fact demo;do PY=$(CE python -c "import doctest;import pyfaust as pf;from os.path import dirname;fn = dirname(pf.__file__)+'/${MOD}.py';print(fn)"); echo -e '\nif __name__ == "__main__":\n    import doctest;\n    doctest.testmod()' >> $PY;done
for MOD in factparams proj poly tools fact demo;do CE coverage run -a --source $PYFAUST_DIR $PYFAUST_DIR/$MOD.py;done
echo > test_mod.py
for T in real complex float32
do
    echo "import pyfaust.tests; pyfaust.tests.run_tests('cpu', '"$T"')" >> test_mod.py
    if echo $CI_RUNNER_TAGS | grep -w cuda
    then
        echo "import pyfaust.tests; pyfaust.tests.run_tests('gpu', '"$T"')" >> test_mod.py
    fi
done
CE coverage run -a --source $PYFAUST_DIR test_mod.py
CE coverage report $(find $PYFAUST_DIR -name "*.py" | grep -v -w tests) | tee /tmp/pyfaust_cov_report
COV=$(sed -e '/^$/d' /tmp/pyfaust_cov_report | tail -1  | sed -e 's/.*\s\+//')
[[ ! $COV =~ ^[[:digit:]]+%$ ]] && echo "Failed to get total coverage" && exit 2
echo Coverage: $COV
CE coverage html $(find $PYFAUST_DIR -name "*.py" | grep -v -w tests)
#conda deactivate
