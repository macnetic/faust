#!/usr/bin/env bash
#
# This script outputs a CHANGELOG.html file that groups commits per tag in time
# order (from the soonest to the oldest) and links them to the gitlab repo.

set -e

TAGS=($(git tag -l --sort=-*authordate))

[[ -n "$DEBUG" ]] && echo "TAGS=${TAGS[*]}"

REPO_URL=https://gitlab.inria.fr/faustgrp/faust/-/commit/

echo > CHANGELOG.html
for I in `seq 0 $((${#TAGS[*]}-2))`
do
	TAG=${TAGS[$I]}
	echo "<h3>Version $TAG:</h3>" >> CHANGELOG.html
	git log ${TAGS[$I]} ^${TAGS[$(($I+1))]} --oneline | sed -e 's/</\&lt;/g;s/>/\&gt;/g;s/$/<br\/>/;s%\([^[:blank:]]\+\)%<a href="'$REPO_URL'\1">\1</a>%' >> CHANGELOG.html
	echo "<hr/>" >> CHANGELOG.html
done
