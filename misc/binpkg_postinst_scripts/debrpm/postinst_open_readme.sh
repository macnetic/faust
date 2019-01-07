#!/bin/bash

# TODO: ideally, it should open the doc but for the moment it indicates the path to it
[[ -r "@CMAKE_INSTALL_PREFIX@/doc/html/index.html" ]] && echo "The documentation is here: @CMAKE_INSTALL_PREFIX@/doc/html/index.html"
