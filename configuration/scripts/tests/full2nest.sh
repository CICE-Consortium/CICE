#!/bin/bash
# you need cdo and nco to run this script
# note the Input section below
###### steps to be taken to test full and nest #####
#1) create full set-up
#       ./cice.setup -c {NAME} -g gbox80 -m {ICE_MACHINE} -e {ICE_MACHINE_ENVNAME} -s bdyrestorfull
#2) build full set-up
#       cd {NAME}; 
#       ./cice.build
#3) run full set-up
#       ./cice.submit
#4) create nest set-up
#       ./cice.setup -c {NAME} -g gbox12 -m {ICE_MACHINE} -e {ICE_MACHINE_ENVNAME} -s bdyrestornest
#5) build nest set-up
#       cd {NAME};
#       ./cice.build
#6) use this file to create all needed fields for running the nest
#       cp full2nest.sh /run/dir/nest/
#7) modify {PATH2hist} and {PATH2rst} to point to the history and restart directory
#8) run full2nest.sh
#       ./full2nest.sh
#9) run nest set-up
#       cd /nest/home/dir
#       ./cice.submit
###############END of description #################

# to be executed in the rundirectory of the nest
#------------------input--------------------------

# case/path to full case
CASE="bdyf21"
PATH2hist="/glade/derecho/scratch/tcraig/CICE_RUNS/${CASE}/history/"
PATH2rst="/glade/derecho/scratch/tcraig/CICE_RUNS/${CASE}/restart/"
# ninest and njnest is the grid size of the nest (do not include halo)
ninest=12
njnest=12
# nis, njs are the full grid indices that align with i=1,j=1 in the nest active grid (not the halo)
nis=42
njs=42

#----------------end input------------------------
#
# make the kmt_grid extended grid from full domain 
# set ncks indexing to fortran (-F)
# need to reduce nis and njs by 1 to extract halo

nis=$(($nis - 1))
njs=$(($njs - 1))
nie=$(($nis + $ninest + 1))
nje=$(($njs + $njnest + 1))
echo "cutting extended grid i=$nis,$nie  j=$njs,$nje"

# cut domain from history file of full domain: 
ncks -F -d ni,$nis,$nie -d nj,$njs,$nje ${PATH2hist}/iceh.2005-01-02.nc temp.nc

#cut relevant variables then convert units and rename variables for grid file
ncrename -v tmask,kmt temp.nc
ncks -C -v kmt,HTN,HTE,NLON,NLAT,ANGLET,ANGLE,TLAT,TLON,ULAT,ULON temp.nc temp1.nc
ncap2 -s "HTE=HTE*100.;HTN=HTN*100.;ULON=ULON/(180/3.14159);ULAT=ULAT/(180/3.14159);TLON=TLON/(180/3.14159);TLAT=TLAT/(180/3.14159)" temp1.nc grid_kmt.nc
ncrename -v HTE,hte -v HTN,htn -v ULON,ulon -v ULAT,ulat -v ANGLE,angle -v ANGLET,anglet -v TLAT,tlat -v TLON,tlon grid_kmt.nc
rm temp*.nc

# cut extended restart files from full domain restart files

for dd in {01..06}
do
        for ss in 03600 07200 10800 14400 18000 21600 25200 28800 32400 36000 39600 43200 46800 50400 54000 57600 61200 64800 68400 72000 75600 79200 82800 00000
        do
                ncks --no-abc -F -d ni,$nis,$nie -d nj,$njs,$nje ${PATH2rst}/iced.2005-01-${dd}-${ss}.nc cice_bdy_restart200501${dd}${ss}.nc
        done
done

