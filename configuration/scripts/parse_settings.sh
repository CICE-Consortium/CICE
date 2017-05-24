#!/bin/bash -f

if [[ "$#" -ne 2 ]]; then
  echo "$0 requires 2 arguments, the filename and filemods"
  exit 0
fi

filename=$1
filemods=$2

echo "$0 $1 $2" 

while read -r line
do
  if [[ "$line" =~ ^\s*$|^\s*#.*|^\s*!.* ]]; then
    echo "skip $line"
  else
    vname=`echo $line | sed "s|\(^\s*set\S*\)\s\{1,100\}\(\S*\)\s\{1,100\}\(\S*\).*$|\2|g"`
    value=`echo $line | sed "s|\(^\s*set\S*\)\s\{1,100\}\(\S*\)\s\{1,100\}\(\S*\).*$|\3|g"`
#    echo "$line $vname $value"

    sed -i 's|\(^\s*set.* '"$vname"' \)[^#]*\(#*.*$\)|\1 '"$value"'  \2|g' $filename
  fi

done < "$filemods"
