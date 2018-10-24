#!/bin/csh -f

if (! -e results.log) then
  echo " "
  echo "${0}: ERROR results.log does not exist, try running results.csh"
  echo " "
  exit -1
endif

set wikirepo = "https://github.com/CICE-Consortium/Test-Results.wiki.git"
set wikiname = Test-Results.wiki

rm -r -f ${wikiname}
git clone ${wikirepo} ${wikiname}
if ($status != 0) then
  echo " "
  echo "${0}: ERROR git clone failed"
  echo " "
  exit -1
endif

set repo = `grep "#repo = " results.log | cut -c 9-`
set bran = `grep "#bran = " results.log | cut -c 9-`
set hash = `grep "#hash = " results.log | cut -c 9-`
set shhash   = `grep "#hshs = " results.log | cut -c 9-`
set hashuser = `grep "#hshu = " results.log | cut -c 9-`
set hashdate = `grep "#hshd = " results.log | cut -c 9-`
set cdat = `grep "#date = " results.log | cut -c 9-`
set ctim = `grep "#time = " results.log | cut -c 9-`
set user = `grep "#user = " results.log | cut -c 9-`
set mach = `grep "#mach = " results.log | cut -c 9-`
set vers = `grep "#vers = " results.log | cut -c 9-`
set totl = `grep "#totl = " results.log | cut -c 9-`
set pass = `grep "#pass = " results.log | cut -c 9-`
set fail = `grep "#fail = " results.log | cut -c 9-`
set cases = `grep -v "#" results.log | grep ${mach}_ | cut -d " " -f 2 | sort -u`
set compilers = `grep -v "#" results.log | grep ${mach}_ | cut -d "_" -f 2 | sort -u`

#echo "debug ${repo}"
#echo "debug ${bran}"
#echo "debug ${hash}"
#echo "debug ${shhash}"
#echo "debug ${hashuser}"
#echo "debug ${hashdate}"
#echo "debug ${cdat}"
#echo "debug ${ctim}"
#echo "debug ${user}"
#echo "debug ${mach}"
#echo "debug ${vers}"
#echo "debug ${totl}"
#echo "debug ${pass}"
#echo "debug ${fail}"
#echo "debug ${cases}"

set xcdat = `echo $cdat | sed 's|-||g' | cut -c 3-`
set xctim = `echo $ctim | sed 's|:||g'`
set shrepo = `echo $repo | tr '[A-Z]' '[a-z]'`

set tsubdir = cice_master
set hfile = "cice_by_hash"
set mfile = "cice_by_mach"
set vfile = "cice_by_vers"
set bfile = "cice_by_bran"
if ("${shrepo}" !~ "*cice-consortium*") then
  set tsubdir = cice_dev
  set hfile = {$hfile}_forks
  set mfile = {$mfile}_forks
  set vfile = {$vfile}_forks
  set bfile = {$bfile}_forks
endif

set noglob
set green  = "\![#00C000](https://placehold.it/15/00C000/000000?text=+)"
set red    = "\![#F00000](https://placehold.it/15/F00000/000000?text=+)"
set orange = "\![#FFA500](https://placehold.it/15/FFA500/000000?text=+)"
set yellow = "\![#FFE600](https://placehold.it/15/FFE600/000000?text=+)"
set gray   = "\![#AAAAAA](https://placehold.it/15/AAAAAA/000000?text=+)"
unset noglob

#==============================================================
# Create results table
#==============================================================

foreach compiler ( ${compilers} )

  set ofile = "${shhash}.${mach}.${compiler}.${xcdat}.${xctim}"
  set outfile = "${wikiname}/${tsubdir}/${ofile}.md"
  mkdir -p ${wikiname}/${tsubdir}
  echo "${0}: writing to ${outfile}"

  if (-e ${outfile}) rm -f ${outfile}

cat >! ${outfile} << EOF

|Bld|Run|Test| Regr | Compare | Timing | Case |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
EOF

@ ttotl = 0
@ tpass = 0
@ tfail = 0
@ tunkn = 0
@ rpass = 0
@ rfail = 0
@ rothr = 0

foreach case ( ${cases} )
if ( ${case} =~ *_${compiler}_* ) then

# check thata case results are meaningful
  set fbuild = `grep " ${case} " results.log | grep " build"   | cut -c 1-4`
  set frun   = `grep " ${case} " results.log | grep " run"     | cut -c 1-4`
  set ftest  = `grep " ${case} " results.log | grep " test"    | cut -c 1-4`

if ( $fbuild != "" || $frun != "" || $ftest != "" ) then

  set fbuild = `grep " ${case} " results.log | grep " build"   | cut -c 1-4`
  set frun   = `grep " ${case} " results.log | grep " run"     | cut -c 1-4`
  set ftest  = `grep " ${case} " results.log | grep " test"    | cut -c 1-4`
  set fregr  = `grep " ${case} " results.log | grep " compare" | cut -c 1-4`
  set fcomp  = `grep " ${case} " results.log | grep " bfbcomp" | cut -c 1-4`
#  if (${ftest}  == "PASS") set frun   = "PASS"
#  if (${frun}   == "PASS") set fbuild = "PASS"

  set vregr  = `grep " ${case} " results.log | grep " compare" | cut -d " " -f 4 | sed 's/\./ /g' `
  set vcomp  = `grep " ${case} " results.log | grep " bfbcomp" | cut -d " " -f 4`

  set vtime1 = `grep " ${case} " results.log | grep " run" | cut -d " " -f 4`
  set vtime2 = `grep " ${case} " results.log | grep " run" | cut -d " " -f 5`
  set vtime3 = `grep " ${case} " results.log | grep " run" | cut -d " " -f 6`

  set btime1 = `grep " ${case} " results.log | grep " compare" | cut -d " " -f 5`
  set btime2 = `grep " ${case} " results.log | grep " compare" | cut -d " " -f 6`
  set btime3 = `grep " ${case} " results.log | grep " compare" | cut -d " " -f 7`

  if (${btime1} != "") then
    if (`echo "${btime1} < 0.0" | bc`) set btime1 = ""
  endif
  if (${btime2} != "") then
    if (`echo "${btime2} < 0.0" | bc`) set btime2 = ""
  endif
  if (${btime3} != "") then
    if (`echo "${btime3} < 0.0" | bc`) set btime3 = ""
  endif

  set vtime  = ""
  if (${vtime1} != "") set vtime = "$vtime TL=${vtime1}(${btime1})"
  if (${vtime2} != "") set vtime = "$vtime Dyn=${vtime2}(${btime2})"
  if (${vtime3} != "") set vtime = "$vtime Col=${vtime3}(${btime3})"

  set scale1 = 1.2
  set scale2 = 1.5
  set ftime  = ""
  if (${vtime1} != "" && ${btime1} != "") then
    if (`echo "${vtime1} > 0.0" | bc` && `echo "${btime1} > 0.0" | bc`) then
      if (`echo "$vtime1 > $btime1*$scale2" | bc`) then
        set ftime = "FAIL"
      else if (`echo "$vtime1 > $btime1*$scale1" | bc`) then
        set ftime = "NOTSOGOOD"
      else
        set ftime = "PASS"
      endif
    endif
  endif

  @ ttotl = $ttotl + 1
  set tchkpass = 1

  set noglob
  set rbuild = ${yellow}
  set rrun   = ${yellow}
  set rtest  = ${yellow}
  set rregr  = ${yellow}
  set rcomp  = ${yellow}
  set rtime  = ${yellow}

  if (${fbuild} == "PASS") set rbuild = ${green}
  if (${frun}   == "PASS") set rrun   = ${green}
  if (${ftest}  == "PASS") set rtest  = ${green}
  if (${fregr}  == "PASS") set rregr  = ${green}
  if (${fcomp}  == "PASS") set rcomp  = ${green}
  if (${ftime}  == "PASS") set rtime  = ${green}

  if (${fbuild} == "FAIL") set rbuild = ${red}
  if (${frun}   == "FAIL") set rrun   = ${red}
  if (${ftest}  == "FAIL") set rtest  = ${red}
  if (${fregr}  == "FAIL") set rregr  = ${red}
  if (${fcomp}  == "FAIL") set rcomp  = ${red}
  if (${ftime}  == "FAIL") set rtime  = ${red}

  if (${fbuild} == "") set rbuild = ${gray}
  if (${frun}   == "") set rrun   = ${red}
  if (${ftest}  == "") set rtest  = ${red}
  if (${fregr}  == "") set rregr  = ${gray}
  if (${fcomp}  == "") set rcomp  = ${gray}
  if (${ftime}  == "") set rtime  = ${gray}

  if (${fbuild} == "COPY") set rbuild = ${gray}
  if (${fbuild} == "MISS") set rbuild = ${gray}
  if (${frun}   == "MISS") set rrun   = ${gray}
  if (${ftest}  == "MISS") set rtest  = ${gray}
  if (${fregr}  == "MISS") set rregr  = ${gray}
  if (${fcomp}  == "MISS") set rcomp  = ${gray}
  if (${ftime}  == "MISS") set rtime  = ${gray}

  if (${rbuild} == ${yellow}) set tchkpass = 2
  if (${rrun}   == ${yellow}) set tchkpass = 2
  if (${rtest}  == ${yellow}) set tchkpass = 2

  if (${rbuild} == ${red}) set tchkpass = 0
  if (${rrun}   == ${red}) set tchkpass = 0
  if (${rtest}  == ${red}) set tchkpass = 0

  if (${tchkpass} == 1) then
     @ tpass = $tpass + 1
  else
    if (${tchkpass} == 2) then
       @ tunkn = $tunkn + 1
    else
       @ tfail = $tfail + 1
    endif
  endif

  if (${rregr} == ${green}) then
     @ rpass = $rpass + 1
  else if (${rregr} == ${red}) then
     @ rfail = $rfail + 1
  else
     @ rothr = $rothr + 1
  endif

  unset noglob

  # remove final .string which is the testid isn't needed here
  set wvcomp = `echo ${vcomp} | sed  's|^\(.*\)\.[^.]*$|\1|g'`
  set xvcomp = `echo ${wvcomp} | sed 's|_| |g'`
  set xcase  = `echo ${case}   | sed 's|_| |g'`
  #echo "debug | ${rbuild} | ${rrun} | ${rtest} | ${rregr} ${vregr} | ${rcomp} ${vcomp} | ${case} |" 
  echo "| ${rbuild} | ${rrun} | ${rtest} | ${rregr} ${vregr} | ${rcomp} ${xvcomp} | ${rtime} ${vtime} | ${xcase} |" >> ${outfile}

endif
endif
end

set noglob
set tcolor = ${green}
if (${tfail} > 0) set tcolor = ${yellow}
@ chk = ((${ttotl} + 9)/ 10)
if (${tfail} >= ${chk}) set tcolor = ${orange}
@ chk = ((${ttotl} + 4) / 5)
if (${tfail} >= ${chk}) set tcolor = ${red}

set rcolor = ${gray}
if (${rfail} > 0 || ${rpass} > 0) then
  if (${rfail} == 0) then
     set rcolor = ${green}
     if (${rothr} > 0) set rcolor = ${yellow}
  endif
  if (${rfail} > 0) set rcolor = ${yellow}
  @ chk = ((${ttotl} + 9)/ 10)
  if (${rfail} >= ${chk}) set rcolor = ${orange}
  @ chk = ((${ttotl} + 4) / 5)
  if (${rfail} >= ${chk}) set rcolor = ${red}
endif
unset noglob

mv ${outfile} ${outfile}.hold
#- raw results: ${totl} total tests: ${pass} pass, ${fail} fail
cat >! ${outfile} << EOF
- repo = **${repo}** : **${bran}**
- hash = ${hash}
- hash created by ${hashuser} ${hashdate}
- vers = ${vers}
- tested on ${mach}, ${compiler}, ${user}, ${cdat} ${ctim} UTC
- ${ttotl} total tests: ${tpass} pass, ${tfail} fail
- ${ttotl} total regressions: ${rpass} pass, ${rfail} fail, ${rothr} other
EOF
cat ${outfile}.hold >> ${outfile}

cat >> ${outfile} << EOF

--------

EOF

#==============================================================

set hashfile = "${wikiname}/${tsubdir}/${hfile}.md"
set machfile = "${wikiname}/${tsubdir}/${mfile}.md"
set versfile = "${wikiname}/${tsubdir}/${vfile}.md"
set branfile = "${wikiname}/${tsubdir}/${bfile}.md"

foreach xfile ($hashfile $machfile $versfile $branfile)
  if (-e ${xfile}) then
    cp -f ${xfile} ${xfile}.prev
  endif
end

#=====================
# update hashfile
#=====================

set chk = 0
if (-e ${hashfile}) set chk = `grep "\*\*${hash}" ${hashfile} | wc -l`
if ($chk == 0) then
cat >! ${hashfile} << EOF
**${hash}** :

| machine | compiler | version | date | test fail | comp fail | total |
| ------ | ------ | ------ | ------  | ------ | ------ | ------ |
| ${mach} | ${compiler} | ${vers} | ${cdat} | ${tcolor} ${tfail}, ${tunkn} | ${rcolor} ${rfail}, ${rothr} | [${ttotl}](${ofile}) |

EOF
if (-e ${hashfile}.prev) cat ${hashfile}.prev >> ${hashfile}

else
  set oline = `grep -n "\*\*${hash}" ${hashfile} | head -1 | cut -d : -f 1`
  @ nline = ${oline} + 3
  sed -i "$nline a | ${mach} | ${compiler} | ${vers} | ${cdat} | ${tcolor} ${tfail}, ${tunkn} | ${rcolor} ${rfail}, ${rothr} | [${ttotl}](${ofile}) | " ${hashfile}
endif

#=====================
# update versfile
#=====================

set chk = 0
if (-e ${versfile}) set chk = `grep "\*\*${vers}" ${versfile} | wc -l`
if ($chk == 0) then
cat >! ${versfile} << EOF
**${vers}** :

| machine | compiler | hash | date | test fail | comp fail | total |
| ------ | ------ | ------ | ------  | ------ | ------ | ------ |
| ${mach} | ${compiler} | ${shhash} | ${cdat} | ${tcolor} ${tfail}, ${tunkn} | ${rcolor} ${rfail}, ${rothr} | [${ttotl}](${ofile}) |

EOF
if (-e ${versfile}.prev) cat ${versfile}.prev >> ${versfile}

else
  set oline = `grep -n "\*\*${vers}" ${versfile} | head -1 | cut -d : -f 1`
  @ nline = ${oline} + 3
  sed -i "$nline a | ${mach} | ${compiler} | ${shhash} | ${cdat} | ${tcolor} ${tfail}, ${tunkn} | ${rcolor} ${rfail}, ${rothr} | [${ttotl}](${ofile}) | " ${versfile}
endif

#=====================
# update machfile
#=====================

set chk = 0
if (-e ${machfile}) set chk = `grep "\*\*${mach}" ${machfile} | wc -l`
if ($chk == 0) then
cat >! ${machfile} << EOF
**${mach}** :

| version | hash | compiler | date | test fail | comp fail | total |
| ------ | ------ | ------ | ------ | ------  | ------ | ------ |
| ${vers} | ${shhash} | ${compiler} | ${cdat} | ${tcolor} ${tfail}, ${tunkn} | ${rcolor} ${rfail}, ${rothr} | [${ttotl}](${ofile}) |

EOF
if (-e ${machfile}.prev) cat ${machfile}.prev >> ${machfile}

else
  set oline = `grep -n "\*\*${mach}" ${machfile} | head -1 | cut -d : -f 1`
  @ nline = ${oline} + 3
  sed -i "$nline a | ${vers} | ${shhash} | ${compiler} | ${cdat} | ${tcolor} ${tfail}, ${tunkn} | ${rcolor} ${rfail}, ${rothr} | [${ttotl}](${ofile}) | " ${machfile}
endif

#=====================
# update branfile
#=====================

set chk = 0
if (-e ${branfile}) set chk = `grep "\*\*${bran}" ${branfile} | wc -l`
if ($chk == 0) then
cat >! ${branfile} << EOF
**${bran}** **${repo}**:

| machine | compiler | hash | date | test fail | comp fail | total |
| ------ | ------ | ------ | ------  | ------ | ------ | ------ |
| ${mach} | ${compiler} | ${shhash} | ${cdat} | ${tcolor} ${tfail}, ${tunkn} | ${rcolor} ${rfail}, ${rothr} | [${ttotl}](${ofile}) |

EOF
if (-e ${branfile}.prev) cat ${branfile}.prev >> ${branfile}

else
  set oline = `grep -n "\*\*${bran}" ${branfile} | head -1 | cut -d : -f 1`
  @ nline = ${oline} + 3
  sed -i "$nline a | ${mach} | ${compiler} | ${shhash} | ${cdat} | ${tcolor} ${tfail}, ${tunkn} | ${rcolor} ${rfail}, ${rothr} | [${ttotl}](${ofile}) | " ${branfile}
endif

#foreach compiler
end

#=====================
# update wiki
#=====================

cd ${wikiname}
git add ${tsubdir}/${shhash}.${mach}*.md
git add ${tsubdir}/${ofile}.md
git add ${tsubdir}/${hfile}.md
git add ${tsubdir}/${mfile}.md
git add ${tsubdir}/${vfile}.md
git add ${tsubdir}/${bfile}.md
git commit -a -m "update $hash $mach"
git push origin master
cd ../

