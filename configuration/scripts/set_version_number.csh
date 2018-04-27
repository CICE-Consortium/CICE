#! /bin/csh -f

if ( $#argv < 1 ) then
  echo "$0 requires one argument, none passed"
  exit -1
endif
if ( $#argv > 1 ) then
  echo "$0 requires one argument, passed = $argv"
  exit -1
endif

set versno = $1
#echo "$0 versno = $versno"

cp -f doc/source/conf.py doc/source/conf.py.bu

sed -i 's|^.*version.*=.*$|version = u'"'"${versno}"'"' | g' doc/source/conf.py 
sed -i 's|^.*release.*=.*$|version = u'"'"${versno}"'"' | g' doc/source/conf.py 

echo "CICE ${versno}"  >! cicecore/version.txt

echo "$0 completed successfully"

exit 0
