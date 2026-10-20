#!/usr/bin/env -S csh -f

# This script can be modified as needed.
# By default, the script is set up to create March and September plots for the global case:
#   ./cice.setup -m derecho -e intel -g gx1 -p 128x1 -s gx1prod,long,run8year --case casename
# The script can be edited to plot output from box cases, e.g.
#   ./cice.setup -m conda -e macos --test smoke -s boxslotcyl --grid gbox80 -p 1x1 --testid 52cb686_remap
# Plots of daily box output are currently set up in the script below.

#### User defined ########################################

source ${MODULESHOME}/init/csh

# Set plot type (box, global, timeseries)
# Set case and case directory
# Set files, notes, fstr, and fields

set plotgrid = "global"
#set plotgrid = "box"
#set plotgrid = "none"
echo "plotgrid = ${plotgrid}"

set plottimeseries = "true"
#set plottimeseries = "false"
echo "plottimeseries = ${plottimeseries}"

if ( ${plotgrid} == 'global' ) then

set case = "CICE6.5.1"
set casedir = "/glade/derecho/scratch/tcraig/CICE_RUNS/cgx1proda"
set histdir = "${casedir}/history"

set files = ("${histdir}/iceh.2012-03.nc" \
             "${histdir}/iceh.2012-09.nc" )
set notes = ("2012 March Mean" \
             "2012 Sept Mean" )
set fstrs = ("Mar12" \
             "Sep12" )

set fields = ("aice" "hi" "hs")

else if ( ${plotgrid} == 'box' ) then

set case = "52cb686_remap"
set casedir = "/Users/eclare/cice-dirs/runs/conda_macos_smoke_gbox80_1x1_boxslotcyl.52cb686_remap"
set histdir = "${casedir}/history"

set files = ("${histdir}/iceh_ic.2005-01-01-00000.nc" \
             "${histdir}/iceh_06h.2005-01-02-00000.nc"\
             "${histdir}/iceh_06h.2005-01-03-00000.nc"\
             "${histdir}/iceh_06h.2005-01-04-00000.nc"\
             "${histdir}/iceh_06h.2005-01-05-00000.nc"\
             "${histdir}/iceh_06h.2005-01-06-00000.nc"\
             "${histdir}/iceh_06h.2005-01-07-00000.nc"\
             "${histdir}/iceh_06h.2005-01-08-00000.nc"\
             "${histdir}/iceh_06h.2005-01-09-00000.nc"\
             "${histdir}/iceh_06h.2005-01-10-00000.nc"\
             "${histdir}/iceh_06h.2005-01-11-00000.nc"\
             "${histdir}/iceh_06h.2005-01-12-00000.nc"\
             "${histdir}/iceh_06h.2005-01-13-00000.nc")
set notes = ("initial condition" \
             "after 1  day " \
             "after 2  day " \
             "after 3  day " \
             "after 4  day " \
             "after 5  day " \
             "after 6  day " \
             "after 7  day " \
             "after 8  day " \
             "after 9  day " \
             "after 10 days" \
             "after 11 days" \
             "after 12 days" )
set fstrs = ("ic" \
             "1day" \
             "2day" \
             "3day" \
             "4day" \
             "5day" \
             "6day" \
             "7day" \
             "8day" \
             "9day" \
             "10days" \
             "11days" \
             "12days" )

set fields = ("aice" "hi")

endif

#### End user defined ########################################

echo " "

# Plot timeseries
if ( ${plottimeseries} == 'true' ) then
    echo ./timeseries.py \"${casedir}\" --case \"${case}\" --grid
    ./timeseries.py "${casedir}" --case "${case}" --grid
endif

# Plot global or box domain
if ( ${plotgrid} == 'global' ) then

#conda config --add channels conda-forge
#conda config --set channel_priority strict
#conda search basemap --channel conda-forge
#conda create -p /glade/u/home/tcraig/conda/envs/basemap -c conda-forge basemap=1.4.1 basemap-data basemap-data-hires netCDF4

module load conda
source ${NCAR_ROOT_CONDA}/etc/profile.d/conda.csh

conda activate /glade/u/home/tcraig/conda/envs/basemap

set cnt = 0
while ($cnt < ${#files})
  @ cnt = $cnt + 1
  set file = "${files[$cnt]}"
  set note = "${notes[$cnt]}"
  set fstr = "${fstrs[$cnt]}"
  foreach field ($fields)
    echo ./ciceplots2d.py \"$field\" \"$file\" \"$case\" \"$note\" \"$fstr\"
    ./ciceplots2d.py "$field" "$file" "$case" "$note" "$fstr"
  end
end

else if ( ${plotgrid} == 'box' ) then

set cnt = 0
while ($cnt < ${#files})
  @ cnt = $cnt + 1
  set file = "${files[$cnt]}"
  set note = "${notes[$cnt]}"
  set fstr = "${fstrs[$cnt]}"
  foreach field ($fields)
    echo ./ciceplots2dbox.py \"$field\" \"$file\" \"$case\" \"$note\" \"$fstr\"
    ./ciceplots2dbox.py "$field" "$file" "$case" "$note" "$fstr"
  end
end

endif

echo " "

echo "DONE"

