#!/bin/bash

function usage
{
	echo "$0 <faust_version> [<filepath>]"
}

function sha256sum_or_not_found
{
	[ -r "$1" ] && sha256sum $1 || echo "file not found"
}

function generate_html
{
	echo "<b>Latest FAµST version: $VERSION ($(date +%D))"
	echo "<table cellpadding="5">"
	echo "<tr><th>System</th><th>Link</th><th>SHA256</th></tr>"
	for CONF in "Linux RPM::faust-$VERSION-x86_64.rpm" "Linux DEB::faust-$VERSION-x86_64.deb" "Windows 10 EXE::faust-$VERSION-amd64.exe" "Mac OS X PKG::faust-$VERSION.pkg"
	do
		SYS=${CONF%%::*}
		FILE=${CONF##*::}
		echo "<tr><td>$SYS</td><td><a href=\"https://faustgrp.gitlabpages.inria.fr/faust/last-doc/packages/$FILE\">$FILE</a></td><td>$(sha256sum_or_not_found $FILE)</td></tr>"
	done
	echo "</table>"

}

[[ $# < 1 ]] && usage && exit 1

VERSION=$1

[[ ! VERSION =~ ^([[:digit:]].?)+$ ]] || (echo "$VERSION is not a valid version tag.")

if [[ $# -gt 3 ]]
then
	FILEPATH=$2
else
	FILEPATH="download.html"
fi

generate_html > $FILEPATH
