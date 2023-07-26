#!/bin/bash

# This script is intented to sphinx doc that uses math dollar extension and needs to replace \f$ of doxygen by $

[[ $# -lt 1 || ! -d "$1" ]] && echo "Error: the script waits a directory (containing .py to filter) as first argument." >&2 && exit 1

for F in $(find $1 -name "*.py")
do
        sed -i 's/\\f\$/$/g;s/\\f\[/$$/g;s/\\f\]/$$/g;s/\\see/See also:/;s/<br\/>//g;s/<b>/**/g;s/<\/b>/**/g;s/&nbsp;//g' $F
done
