#!/bin/bash -f

if [[ "$#" -ne 2 ]]; then
  echo "$0 ERROR: requires 2 arguments, the filename and filemods"
  exit -1
fi

scriptname=`basename "$0"`
filename=$1
filemods=$2

#echo "$0 $1 $2" 
echo "running parse_namelist.sh"
foundstring="FoundSTRING"
vnamearray=()
valuearray=()

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

    found=${foundstring}
    for i in "${!vnamearray[@]}"; do
      if [[ "${found}" == "${foundstring}" ]]; then
        vn=${vnamearray[$i]}
        vv=${valuearray[$i]}
#        echo "names/values $i ${vname} ${vn} ${value} ${vv}"
        if [[ "$vname" == "$vn" ]]; then
          found=$i
          if [[ "$value" != "${vv}" ]]; then
#            echo "names/values $i ${vname} ${vn} ${value} ${vv}"
            echo "${scriptname} WARNING: re-overriding $vname from ${vv} to ${value}"
          fi
        fi
      fi
    done

    #sed -i 's|\(^\s*'"$vname"'\s*\=\s*\)\(.*$\)|\1'"$value"'|g' $filename
    cp ${filename} ${filename}.check
    sed -i.sedbak -e 's|\(^[[:space:]]*'"$vname"'[[:space:]]*=[[:space:]]*\)\(.*$\)|\1'"$foundstring"'|g' ${filename}.check
    grep -q ${foundstring} ${filename}.check
    if [ $? -eq 0 ]; then
      sed -i.sedbak -e 's|\(^[[:space:]]*'"$vname"'[[:space:]]*=[[:space:]]*\)\(.*$\)|\1'"$value"'|g' ${filename}
      if [[ "${found}" == "${foundstring}" ]]; then
        vnamearray+=($vname)
        valuearray+=($value)
      else
        valuearray[$found]=${value}
      fi
      if [[ -e "${filename}.sedbak" ]]; then
        rm ${filename}.sedbak
      fi
    else
      echo "${scriptname} ERROR: parsing error for ${vname}"
      exit -99
    fi
    rm ${filename}.check ${filename}.check.sedbak

  fi

done < "$filemods"

exit 0
