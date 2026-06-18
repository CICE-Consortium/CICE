#!/bin/csh -f

#source ${MODULESHOME}/init/csh

# User defined stuff
# Set case and case directory
# Set files, notes, fstr, and fields
# Set plot type (box, global, timeseries)

set case = "CICEv6"
set casedir = "/Users/eclare/cice-dirs/runs/conda_macos_smoke_gbox80"
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

set fields = ("aice" "hi" "hs")

set gridtype = "box"
#set gridtype = "global"
echo "gridtype = ${gridtype}"

set plottimeseries = "false"
#set plottimeseries = "true"
echo "plottimeseries = ${plottimeseries}"

echo " "

# Plot global or box domain
if ( ${gridtype} == 'global' ) then

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

else if ( ${gridtype} == 'box' ) then

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

# Plot timeseries
if ( ${plottimeseries} == 'true' ) then
    echo ./timeseries.py \"${casedir}\" --case \"${case}\" --grid
    ./timeseries.py "${casedir}" --case "${case}" --grid
endif

echo "DONE"

