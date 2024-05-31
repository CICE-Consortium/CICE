#!/bin/csh -f

source ${MODULESHOME}/init/csh

# User defined stuff
# Set case and case directory
# Set files, notes, fstr, and fields

set case = "CICE6.5.1"
set casedir = "/glade/derecho/scratch/tcraig/CICE_RUNS/cgx1proda"

# setup plots

set histdir = "${casedir}/history"

set files = ("${histdir}/iceh.2012-03.nc" \
             "${histdir}/iceh.2012-09.nc" )
set notes = ("2012 March Mean" \
             "2012 Sept Mean" )
set fstrs = ("Mar12" \
             "Sep12" )

set fields = ("aice" "hi" "hs")

#conda config --add channels conda-forge
#conda config --set channel_priority strict
#conda search basemap --channel conda-forge
#conda create -p /glade/u/home/tcraig/conda/envs/basemap -c conda-forge basemap=1.4.1 basemap-data basemap-data-hires netCDF4

module load conda
source ${NCAR_ROOT_CONDA}/etc/profile.d/conda.csh

conda activate /glade/u/home/tcraig/conda/envs/basemap

echo " "
echo " "

echo ./timeseries.py \"${casedir}\" --case \"${case}\" --grid
./timeseries.py "${casedir}" --case "${case}" --grid

echo " "

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

echo "DONE"

