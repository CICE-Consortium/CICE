
#--- cice.codecov.csh ---

#if ( ${use_curl} == 1 ) then
#  bash -c "bash <(curl -s https://codecov.io/bash) -n '${report_name}' -y ./codecov.yml "
#else
#  bash -c "bash <(wget -O - https://codecov.io/bash) -n '${report_name}' -y ./codecov.yml "
#endif

if ( ${use_curl} == 1 ) then
  curl https://codecov.io/bash -o codecov.bash
else
  wget https://codecov.io/bash -O codecov.bash
endif
chmod +x codecov.bash
sed -i.sedbak 's|mktemp /tmp/|mktemp ./|g' codecov.bash
bash -c "bash ./codecov.bash -n '${report_name}' -y ./codecov.yml"

sleep 10
rm -r -f ./*/codecov_output
