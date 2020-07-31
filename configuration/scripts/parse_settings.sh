#!/bin/bash -f

if [[ "$#" -ne 2 ]]; then
  echo "$0 ERROR: requires 2 arguments, the filename and filemods"
  exit -1
fi

scriptname=`basename "$0"`
filename=$1
filemods=$2

#echo "$0 $1 $2" 
echo "running parse_settings.sh"
foundstring="FoundSTRING"
vnamearray=()
valuearray=()

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

    #sed -i 's|\(^\s*set.* '"$vname"' \)[^#]*\(#*.*$\)|\1 '"$value"'  \2|g' $filename
    sed -i.sedbak -e 's|\(^[[:space:]]*set.* '"$vname"' \)[^#]*\(#*.*$\)|\1 '"$value"'  \2|g' $filename

    if [[ "${found}" == "${foundstring}" ]]; then
      vnamearray+=($vname)
      valuearray+=($value)
    else
      valuearray[$found]=${value}
    fi

    if [[ -e "${filename}.sedbak" ]]; then
      rm ${filename}.sedbak
    fi

  fi

done < "$filemods"

exit 0
