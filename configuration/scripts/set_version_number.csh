#!/bin/csh -f

if ( $#argv < 1 ) then
  echo "$0 requires one argument, none passed"
  exit -1
endif
if ( $#argv > 1 ) then
  echo "$0 requires one argument, passed = $argv"
  exit -1
endif

set versno = $1
set cdate = `date +%Y-%m-%d`
#echo "$0 versno = $versno"

cp -f doc/source/conf.py doc/source/conf.py.bu

sed -i 's|^.*version.*=.*$|version = u'"'"${versno}"'"' | g' doc/source/conf.py 
sed -i 's|^.*release.*=.*$|version = u'"'"${versno}"'"' | g' doc/source/conf.py 

cp -f .zenodo.json .zenodo.json.bu

sed -i 's|^\(.*CICE-Consortium/CICE:\).*$|\1 CICE Version '${versno}'", | g' .zenodo.json
sed -i 's|^\(.*"version":\).*$|\1 "'${versno}'", | g' .zenodo.json
sed -i 's|^\(.*"publication_date":\).*$|\1 "'${cdate}'", | g' .zenodo.json

echo "CICE ${versno}"  >! cicecore/version.txt

echo "$0 completed successfully"

exit 0
