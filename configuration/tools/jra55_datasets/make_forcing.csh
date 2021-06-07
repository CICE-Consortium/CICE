#!/bin/csh
# -----
# This is a script that worked on NCAR's cheyenne in March, 2021.
# It converts raw JRA55 datasets to a format that CICE can use.
# This tools is documented in the CICE user guide.  The
# tool interpolates to a CICE grid and does things like convert units.
# -----
# The interp_jra55_ncdf_bilinar.py script was placed in "scripts_dir"
# The raw JRA55 datasets were placed in "jra55_data_dir"
# The CICE grid files were places in "jra55_data_dir"
# The model output was created in "output_data_dir"
# -----
#PBS -N make_forcing
#PBS -q regular
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -l walltime=06:00:00
#PBS -A P93300665

set scripts_dir = "/glade/work/tcraig/cice-consortium/cice.jra55_tool/configuration/tools/jra55_datasets"
set jra55_data_dir = "/glade/scratch/dbailey/JRA_DATA/"
set output_data_dir = "/glade/scratch/tcraig/JRA_DATA_output"
set grid = "gx3"
set cice_grid_file = "grid_gx3.nc"

module load python/3.7.9
source /glade/u/apps/opt/ncar_pylib/ncar_pylib.csh default
module load nco

mkdir -p ${output_data_dir}
cd ${output_data_dir}

ln -s ${jra55_data_dir}/fcst_*.nc .
ln -s ${jra55_data_dir}/grid_*.nc .

ln -s ${scripts_dir}/interp_jra55_ncdf_bilinear.py .

#foreach year ( 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 )
foreach year ( 1997 )

./interp_jra55_ncdf_bilinear.py ${year}010100_${year}033121 ${cice_grid_file} JRA55_${grid}_03hr_forcing_${year}-q1.nc
./interp_jra55_ncdf_bilinear.py ${year}040100_${year}063021 ${cice_grid_file} JRA55_${grid}_03hr_forcing_${year}-q2.nc
./interp_jra55_ncdf_bilinear.py ${year}070100_${year}093021 ${cice_grid_file} JRA55_${grid}_03hr_forcing_${year}-q3.nc
./interp_jra55_ncdf_bilinear.py ${year}100100_${year}123121 ${cice_grid_file} JRA55_${grid}_03hr_forcing_${year}-q4.nc

ncrcat JRA55_${grid}_03hr_forcing_${year}-??.nc JRA55_${grid}_03hr_forcing_${year}.nc

/bin/rm -f ${jra55_data_dir}/JRA55_${grid}_03hr_forcing_${year}-??.nc

end
