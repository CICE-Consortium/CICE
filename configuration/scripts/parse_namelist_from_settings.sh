#!/bin/bash -f

if [[ "$#" -ne 2 ]]; then
  echo "$0 ERROR: requires 2 arguments, the namelist and settings"
  exit -1
fi

filename=$1
filemods=$2

#echo "$0 $1 $2" 
echo "running parse_namelist_from_settings.sh"

sed -i.sedbak -e 's|ICE_SANDBOX|'"${ICE_SANDBOX}"'|g' $filename
sed -i.sedbak -e 's|ICE_MACHINE_INPUTDATA|'"${ICE_MACHINE_INPUTDATA}"'|g' $filename

exit 0
