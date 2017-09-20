#!/bin/bash -f

if [[ "$#" -ne 2 ]]; then
  echo "$0 ERROR: requires 2 arguments, the filename and filemods"
  exit -1
fi

filename=$1
filemods=$2

#echo "$0 $1 $2" 
echo "running parse_namelist.sh"

while read -r line
do
  if [[ "$line" =~ ^\s*$|^\s*#.*|^\s*!.* ]]; then
#    echo "skip $line"
     :
  else
    #vname=`echo $line | sed "s|\s*\(\S*\)\s*=\s*\(\S*\).*$|\1|g"`
    #value=`echo $line | sed "s|\s*\(\S*\)\s*=\s*\(\S*\).*$|\2|g"`
    vname=`echo $line | sed "s|^[[:space:]]*\([^[:space:]]*\)[[:space:]]*=[[:space:]]*\([^[:space:]]*\).*$|\1|g"`
    value=`echo $line | sed "s|^[[:space:]]*\([^[:space:]]*\)[[:space:]]*=[[:space:]]*\([^[:space:]]*\).*$|\2|g"`
#    echo "$line $vname $value"

    #sed -i 's|\(^\s*'"$vname"'\s*\=\s*\)\(.*$\)|\1'"$value"'|g' $filename
    sed -i.sedbak -e 's|\(^[[:space:]]*'"$vname"'[[:space:]]*=[[:space:]]*\)\(.*$\)|\1'"$value"'|g' $filename
    if [[ -e "${filename}.sedbak" ]]; then
      rm ${filename}.sedbak
    fi
  fi

done < "$filemods"

exit 0
