#!/bin/sh

# DEPRECATION: this script was intended to upload FAµST packages on another
# gitlab project but using package registry is a way better

if [[ -z "$UP_TOKEN" ]]
then
	echo UP_TOKEN variable is not set in the environment
	exit 1
fi

# on each new tag the branch packages (of faust-packages project) is deleted and created again to avoid storing too many packages
# UP_TOKEN must be configured in CI/CD variables or as environment variable
curl --request DELETE --header "PRIVATE-TOKEN: $UP_TOKEN"  "https://gitlab.inria.fr/api/v4/projects/29950/repository/branches/packages"
curl --request POST --header "PRIVATE-TOKEN: $UP_TOKEN" "https://gitlab.inria.fr/api/v4/projects/29950/repository/branches?branch=packages&ref=master"

for PKG in build/faust-${CI_COMMIT_TAG}.pkg build/faust-$CI_COMMIT_TAG-amd64.exe build/faust-$CI_COMMIT_TAG-x86_64.*
do
	echo
	curl --request POST --form "branch=packages"  --form "commit_message=Add package $(basename $PKG)"  --form "start_branch=packages"  --form "actions[][action]=create"  --form "actions[][file_path]=$(basename $PKG)"  --form "actions[][content]=@$PKG"  --header "PRIVATE-TOKEN: $UP_TOKEN" "https://gitlab.inria.fr/api/v4/projects/29950/repository/commits" || break
done
