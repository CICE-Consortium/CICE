#!/bin/bash -f

if [[ "$#" -ne 2 ]]; then
  echo "$0 ERROR: requires 2 arguments, the filename and filemods"
  exit -1
fi

filename=$1
filemods=$2

#echo "$0 $1 $2" 
echo "running parse_namelist.sh"
foundstring="FoundSTRING"

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
    cp ${filename} ${filename}.check
    sed -i -e 's|\(^[[:space:]]*'"$vname"'[[:space:]]*=[[:space:]]*\)\(.*$\)|\1'"$foundstring"'|g' ${filename}.check
    grep -q ${foundstring} ${filename}.check
    if [ $? -eq 0 ]; then
      sed -i.sedbak -e 's|\(^[[:space:]]*'"$vname"'[[:space:]]*=[[:space:]]*\)\(.*$\)|\1'"$value"'|g' ${filename}
      if [[ -e "${filename}.sedbak" ]]; then
        rm ${filename}.sedbak
      fi
    else
      echo "$0 ERROR: parsing error for ${vname}"
      exit -99
    fi
    rm ${filename}.check

  fi

done < "$filemods"

exit 0
