:tocdepth: 3

.. _tools:

Tools
=============


.. _cice4restart:

CICE4 restart conversion
-------------------------

There is a Fortran program in **configuration/tools/cice4_restart_conversion**
that will help convert a CICE4 restart file into a CICE5 restart file.
There is a bit of documentation contained in that source code about how
to build, use, and run the tool.  A few prognostic variables were changed
from CICE4 to CICE5 which fundamentally altered the fields saved to
the restart file.  See 
**configuration/tools/cice4_restart_conversion/convert_restarts.f90** 
for additional information.


.. _jra55datasettool:

JRA55 forcing datasets
------------------------

This section describes how to generate JRA55 forcing data for the CICE model.
Raw JRA55 files have to be interpolated and processed into input files specifically
for the CICE model.  A tool exists in **configuration/tools/jra55_datasets**
to support that process.
The raw JRA55 data is obtained from the NCAR/UCAR Research Data Archive and
the conversion tools are written in python.

Requirements
*********************

Python3 is required, and the following
python packages are required with the tested version number in parenthesis.  These
versions are not necessarily the only versions that work, they just indicate what 
versions were used when the script was recently run.

- python3 (python3.7.9)
- numpy (1.18.5)
- netCDF4 (1.5.5)
- ESMPy (8.0.0)
- xesmf (0.3.0)

NCO is required for aggregating the output files into yearly files.

- netcdf (4.7.4)
- nco (4.9.5)

Raw JRA55 forcing data
*************************

The raw JRA55 forcing data is obtained from the UCAR/NCAR Research Data Archive,
https://rda.ucar.edu/.  You must first register (free) and then sign in.  The
"JRA-55 Reanalysis Daily 3-Hourly and 6-Hourly Data" is ds628.0 and can be found here,
https://rda.ucar.edu/datasets/ds628.0.  

The "Data access" tabs will provide a list of product categories.
The JRA55 data of interest are located in 2 separate products. Winds, air 
temperature, and specific humidity fields are included in "JRA-55 
3-Hourly Model Resolution 2-Dimensional Instantaneous Diagnostic Fields". 
Precipitation and downward radiation fluxes are found in "JRA-55 3-Hourly 
Model Resolution 2-Dimensional Average Diagnostic Fields".  (Note the 
difference between instantaneous and averaged data products. There are several 
JRA55 datasets available, you will likely have to scroll down the page to find 
these datasets.) Data are also available on a coarser 1.25° grid, but the tools
are best used with the native TL319 JRA55 grid.

The fields needed for CICE are

- specific humidity (3-hourly instantaneous), Qa
- temperature (3-hourly instantaneous), Tair
- u-component of wind (3-hourly instantaneous), uatm
- v-component of wind(3-hourly instantaneous), vatm
- downward longwave radiation flux (3 hourly average), flw
- downward solar radiation flux (3 hourly average), fsw
- total precipitation (3 hourly average), fsnow

To customize the dataset for download, choose the “Get a Subset” option. Select 
the desired times in the “Temporal Selection” section, then click on desired parameters
(see list above).  After clicking continue, select Output Format "Converted to NetCDF".

Once the data request is made, an email notification will be sent with a dedicated
URL that will provide a variety of options for downloading the data remotely.
The data will be available to download for 5 days.  
The raw data consists of multiple files, each containing three months of data for
one field.


Data conversion
*************************

The script, **configuration/tools/jra55_datasets/interp_jra55_ncdf_bilinear.py**, 
converts the raw data to CICE input files.

The script uses a bilinear regridding algorithm to regrid from the JRA55 grid to 
the CICE grid. The scripts use the Python package ‘xesmf’ to generate bilinear 
regridding weights, and these regridding weights are written to the file defined by
the variable "blin_grid_name" in **interp_jra55_ncdf_bilinear.py**. This filename
can be modified by editing **interp_jra55_ncdf_bilinear.py**.
The weights file can be re-used if interpolating different data on the same grid. 
Although not tested in this version of the scripts, additional regridding options 
are available by xesmf, including ‘conservative’ and ‘nearest neighbor’. These 
methods have not been tested in the current version of the scripts. The reader 
is referred to the xESMF web page for further documentation 
(https://xesmf.readthedocs.io/en/latest/ last accessed 5 NOV 2020).

To use the **interp_jra55_ncdf_bilinear** script, do ::

  python3 interp_jra55_ncdf_bilinear.py –h

to see the latest interface information ::

  usage: interp_jra55_ncdf_bilinear.py [-h] JRADTG gridout ncout

  Interpolate JRA55 data to CICE grid

  positional arguments:
    JRADTG      JRA55 input file date time group
    gridout     CICE grid file (NetCDF)
    ncout       Output NetCDF filename

  optional arguments:
    -h, --help  show this help message and exit

Sample usage is ::

  ./interp_jra55_ncdf_bilinear.py 1996010100_1996033121 grid_gx3.nc JRA55_gx3_03hr_forcing_1996-q1.nc
  ./interp_jra55_ncdf_bilinear.py 1996040100_1996063021 grid_gx3.nc JRA55_gx3_03hr_forcing_1996-q2.nc
  ./interp_jra55_ncdf_bilinear.py 1996070100_1996093021 grid_gx3.nc JRA55_gx3_03hr_forcing_1996-q3.nc
  ./interp_jra55_ncdf_bilinear.py 1996100100_1996123121 grid_gx3.nc JRA55_gx3_03hr_forcing_1996-q4.nc

In this case, the 4 quarters of 1996 JRA55 data is going to be interpolated to the gx3 grid.
NCO can be used to aggregate these files into a single file ::

  ncrcat JRA55_gx3_03hr_forcing_1996-??.nc JRA55_${grid}_03hr_forcing_1996.nc

NOTES

- The scripts are designed to read a CICE grid file in netCDF format.  This is the "grid_gx3.nc" file above.  The NetCDF grid names are hardcoded in **interp_jra55_ncdf_bilinear.py**. If you are using a different grid file with different variable names, this subroutine needs to be updated. 
- All files should be placed in a common directory.  This includes the raw JRA55 input files, the CICE grid file, and **interp_jra55_ncdf_bilinear.py**.  The output files will be written to the same directory.
- The script **configuration/tools/jra55_datasets/make_forcing.csh** was used on the NCAR cheyenne machine in March, 2021 to generate CICE forcing data.  It assumes the raw JRA55 is downloaded, but then sets up the python environment, links all the data in a common directory, runs **interp_jra55_ncdf_bilinear.py** and then aggregates the quarterly data using NCO.
- The new forcing files can then be defined in the **ice_in** namelist file using the input variables, ``atm_data_type``, ``atm_data_format``, ``atm_data_dir``, ``fyear_init``, and ``ycycle``.  See :ref:`forcing` for more information.
- The total precipitation field is mm/day in JRA55.  This field is initially read in as snow, but prepare_forcing in **ice_forcing.F90** splits that into rain or snow forcing depending on the air temperature.

