#!/bin/bash -f

if [[ "$#" -ne 2 ]]; then
  echo "$0 ERROR: requires 2 arguments, the filename and filemods"
  exit -1
fi

filename=$1
filemods=$2

#echo "$0 $1 $2" 
echo "running parse_settings.sh"

while read -r line
do
  if [[ "$line" =~ ^\s*$|^\s*#.*|^\s*!.* ]]; then
#    echo "skip $line"
     :
  else
    #vname=`echo $line | sed "s|\(^\s*set\S*\)\s\{1,100\}\(\S*\)\s\{1,100\}\(\S*\).*$|\2|g"`
    #value=`echo $line | sed "s|\(^\s*set\S*\)\s\{1,100\}\(\S*\)\s\{1,100\}\(\S*\).*$|\3|g"`
    vname=`echo $line | sed "s|\(^[[:space:]]*set[^[:space:]]*\)[[:space:]][[:space:]]*\([^[:space:]]*\)[[:space:]][[:space:]]*\([^[:space:]]*\).*$|\2|g"`
    value=`echo $line | sed "s|\(^[[:space:]]*set[^[:space:]]*\)[[:space:]][[:space:]]*\([^[:space:]]*\)[[:space:]][[:space:]]*\([^[:space:]]*\).*$|\3|g"`
#    echo "$line $vname $value"

    #sed -i 's|\(^\s*set.* '"$vname"' \)[^#]*\(#*.*$\)|\1 '"$value"'  \2|g' $filename
    sed -i.sedbak -e 's|\(^[[:space:]]*set.* '"$vname"' \)[^#]*\(#*.*$\)|\1 '"$value"'  \2|g' $filename
    if [[ -e "${filename}.sedbak" ]]; then
      rm ${filename}.sedbak
    fi

  fi

done < "$filemods"

exit 0
