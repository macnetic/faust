#!/bin/bash

echo "$0 -- FAµST doc launcher post-install script START" | tee -a /tmp/log_faust_install
[[ -r "@CMAKE_INSTALL_PREFIX@/doc/html/index.html" ]] && open @CMAKE_INSTALL_PREFIX@/doc/html/index.html
echo "$0 -- FAµST doc launcher post-install script END" | tee -a /tmp/log_faust_install
