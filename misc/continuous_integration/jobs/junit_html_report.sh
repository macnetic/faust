#!/usr/bin/env bash

[[ $# -lt 1 ]] && echo "I need a report name" && exit 1
REPORT_NAME="$1"
python3 -m venv junit_venv # use default python3 in PATH (the same has been used in ctest building -- cf. CDashConfScript)
source junit_venv/bin/activate
export BUILD_DIR=$(ls -d build  | tail -1)
pip install junit2html
junit2html ${BUILD_DIR}/junit_output.xml ${BUILD_DIR}/${REPORT_NAME}.html
deactivate
rm -Rf junit_venv
echo The test report will be published here: https://gitlab.inria.fr/faustgrp/faust/-/jobs/${CI_JOB_ID}/artifacts/external_file/${BUILD_DIR}/${REPORT_NAME}.html
