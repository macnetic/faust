#!/bin/bash

# This script is for sphinx doc that:
# - uses math dollar extension and needs to replace \f$ of doxygen by $
# - many other specialized syntaxes.

[[ $# -lt 1 || ! -d "$1" ]] && echo "Error: the script waits a directory (containing .py to filter) as first argument." >&2 && exit 1

for F in $(find $1 -name "*.py")
do
        sed -i 's/\\f\$/$/g;s/\\f\[/$$/g;s/\\f\]/$$/g;s/\\see/.. seealso::/;s/<br\/>//g;s/<b>/**/g;s/<\/b>/**/g;s/&nbsp;//g' $F
        sed -i 's/<a href="\([^"]\+\)">\([^<]\+\)<\/a>/`\2 <\1>`_/g' $F # full link on one line
        sed -i 's/NOTE:\|Note:/.. note::/g' $F
        sed -i 's/Note:/.. note::/g' $F
        sed -i 's/[wW]arning:\|WARNING:/.. warning::/g' $F
done
