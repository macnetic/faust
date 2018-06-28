#!/bin/bash

# for debian and ubuntu packages libmatio-dev are only creating libmatio.so file
# but faust build relies on a symblink toward this lib which is named libmatio.so.4
# so create it!

which dpkg-query >/dev/null && LIBMATIO_SO_PATH=$(dpkg-query -L libmatio-dev 2>/dev/null | grep libmatio.so)
ln -sf $LIBMATIO_SO_PATH ${LIBMATIO_SO_PATH}.4



