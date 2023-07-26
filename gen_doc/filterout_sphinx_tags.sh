#!/usr/bin/env bash

set -eu


# Filterout sphinx tags from doxygen html rendered html files.

[[ $# -lt 1 || ! -d "$1" ]] && echo "Error: the script waits a directory (containing .py to filter) as first argument." >&2 && exit 1

for F in $(find $1 -name "*.py")
do
	#sed -i 's/\(:py:func:\|:py:class:\)`\(\.\)\?\([^`]\{1,\}\)`/\3/g' $F # certain sed versions on macOS don't handle \?, \| and -i flag
	sed -e 's/\(:py:func:\)`\(\.\)\{0,1\}\([^`]\{1,\}\)`/\3/g;s/\(:py:class:\)`\(\.\)\{0,1\}\([^`]\{1,\}\)`/\3/g' < $F > /tmp/$(basename $F)
	mv /tmp/$(basename $F) $F
done
