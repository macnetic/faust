#!/usr/bin/env bash

# Filter html files to prevent >> (matlab) or >>> (python) to be selected in order to ease the copy of example for users
# It is inspired from a sphinx feature

[[ $# -lt 1 || ! -d "$1" ]] && echo "Error: the script waits a directory (containing .html to filter) as first argument." >&2 && exit 1

for F in $(find $1 -name "*.html")
do
	#sed -i 's/\&gt;\&gt;\(\&gt;\)\?/<span style="user-select: none;">&<\/span>/g' $F
	sed -i 's/<div class="line">/<div class="line" style="user-select: none;">/g' $F
	#sed -i 's/\(\&gt;\&gt;\&gt;\|\&gt;\&gt;\)\(.*\)<\/div>/\1<span style="user-select: auto;">\2<\/span><\/div>/g' $F
	sed -i 's/<div class="line" style="user-select: none;">\(<span class="stringliteral">\)\?\(\&gt;\&gt;\&gt;\|\s\?\&gt;\&gt;\)\(.*\)/<div class="line">\1<span style="user-select: none;">\2<\/span>\3/g' $F
done
