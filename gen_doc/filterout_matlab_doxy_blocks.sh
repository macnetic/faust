#!/bin/bash

# this script removes the doxygen block of code doc from a matlab script
# it is very useful because otherwise it will be used as inline documentation in matlab (with all specific syntax of doxygen, nobody wants that).
# even if an inner function code doc is existing
# A doxygen code block starts with % followed in the same line by at least five `=' characters.
# It terminates with the same form of sequence.

[[ ! -r "$1" ]] && echo "$1 is not readable/existing as a file." >&2 && exit 0 # not an error, just skipping
#sed -i '/%=\{5,\}/,/%=\{5,\}/d' $1 # macos default sed can't do this (it gives an error:sed: 1: command c expects \ followed by text
sed -e '/%=\{5,\}/,/%=\{5,\}/d' $1 > ${1}_tmp
mv "${1}_tmp" "$1"
