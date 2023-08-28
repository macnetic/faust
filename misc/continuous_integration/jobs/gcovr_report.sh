#!/usr/bin/env bash

python3 -m venv gcovr_venv # use default python3 in PATH (the same has been used in ctest building -- cf. CDashConfScript)
source gcovr_venv/bin/activate
export BUILD_DIR=$(ls -d build* | tail -1)
pip install gcovr
[[ ! -d ctest_coverage ]] && mkdir ctest_coverage
gcovr -r ./src/ ${BUILD_DIR} --html-details ctest_coverage/ctest_coverage.html
echo Coverage: $(grep -A 3 Lines: ctest_coverage/ctest_coverage.html | tail -1 | sed -e "s/[^[:digit:]]\+//;s/[^%]\+$//")
deactivate
rm -Rf gcovr_venv
echo The test coverage report will be published here: https://gitlab.inria.fr/faustgrp/faust/-/jobs/${CI_JOB_ID}/artifacts/external_file/ctest_coverage/ctest_coverage.html
