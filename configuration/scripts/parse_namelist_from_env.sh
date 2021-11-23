#!/bin/bash -f

if [[ "$#" -ne 1 ]]; then
  echo "$0 ERROR: requires 1 argument, the namelist to modify"
  exit -1
fi

scriptname=`basename "$0"`
filename=$1

#echo "$0 $1" 
echo "running $scriptname"

sed -i -e 's|ICE_SANDBOX|'"${ICE_SANDBOX}"'|g' $filename
sed -i -e 's|ICE_MACHINE_INPUTDATA|'"${ICE_MACHINE_INPUTDATA}"'|g' $filename

exit 0
