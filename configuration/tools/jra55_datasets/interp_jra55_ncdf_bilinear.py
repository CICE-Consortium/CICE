#! /usr/bin/env python3

import xesmf as xe
from netCDF4 import Dataset
import argparse
import os
import numpy as np
from datetime import datetime

######################################################
######################################################
def make_regridder(lon1, lat1, lon2, lat2, method, periodic, grdname,
                   lon1_b=None, lat1_b=None, lon2_b=None, lat2_b=None):
    '''
    make nearest neighbor xESMF regridder object.
    input:
    lon1: source longitudes (degrees)
    lat1: source latitudes  (degrees)
    lon2: target longitudes (degrees)
    lat2: target latitudes (degrees)
    method: regridding  method (bilinear, patch, conservative, nearest_s2d)
    periodic: True if periodic longitudes, false if not
    grdname: filename for regridder (Ugrid or Tgrid)
    '''
    if method != "conservative":
        # define grids for regridder
        grid1 = {'lon' : lon1, 'lat' : lat1}
        grid2 = {'lon' : lon2, 'lat' : lat2}


    else:
        # conservative needs boundary lon/lat
        grid1 = {'lon'   : lon1,   'lat'   : lat1,
                 'lon_b' : lon1_b, 'lat_b' : lat1_b}

        grid2 = {'lon'   : lon2,   'lat'   : lat2,
                 'lon_b' : lon2_b, 'lat_b' : lat2_b}

    # make regridder
    # here specify reuse_weights=False to re-generate weight file.
    # if wanted to reuse file inteas of making int, 
    # check if file exists and change use_file_weights=True. 
    # see commented out example below 
    use_file_weights=False

    # check if regrid file exists.
    # If so, reuse file instead of regenerating.
    # if (os.path.isfile(blin_grid_name)):
    #     use_file_weights = True

    regridder = xe.Regridder(ds_in=grid1,ds_out=grid2,
                             method=method,
                             periodic=periodic,
                             filename=grdname,
                             reuse_weights=use_file_weights)


    return regridder

#########################################
#########################################
def halo_extrapolate(a,ew_bndy_type,ns_bndy_type):
    '''
    Extrapolate to 'halo' cell as in CICE code
    ice_boundary.F90:ice_HaloExtrapolate.
    inputs:
    a: array nx+1, ny+1 (nghost/nhalo hard-coded as 1 for now)
    ew_bndy_type: east/west boundary type (cyclic, regional, etc)
    ns_bndy_type: norh/south boundary type (cyclic, regional, etc)

    return: a with halo applied
    '''

    # get dimension of a
    # expected to be 0:nx+nghost, 0:ny+nghost
    nj, ni = a.shape # note with Python NetCDF is nj, ni order
    # W/E edges
    if ew_bndy_type == 'cyclic':
        a[: ,0] = a[:,-2] # -2, since -1 is ghost cell
        a[:,-1] = a[:, 1] #  1, since  0 is ghost cell
    else: # if (trim(ew_bndy_type) /= 'cyclic') then
        a[:, 0] = 2.0*a[:, 1] - a[:, 2]
        a[:,-1] = 2.0*a[:,-2] - a[:,-3]

    # south edge
    if ns_bndy_type == 'cyclic':
        a[0,:] = a[-2,:] # -2, since -1 is ghost cell
    else:
        a[0,:] = 2.0*a[1,:] - a[2,:]

    # north edge treated a little different, depending
    # on if bndy type is tripole
    if ns_bndy_type == 'cyclic':
        a[-1,:] = a[1,:] # 1, since 0 is ghost cell

    elif (ns_bndy_type != 'cyclic'  and
          ns_bndy_type != 'tripole' and
          ns_bndy_type != 'tripoleT'):

        a[-1,:] = 2.0*a[-2,:] - a[-3,:]

    else:
        pass # do nothing

    # return array with halo upated
    return a

#########################################
#########################################

def Tlatlon(ulat,ulon,ew_bndy_type,ns_bndy_type):
    '''
    Make TLAT/TLON from ULAT/ULON.
    see ice_grid.F90:Tlatlon for method
    Inputs:
    ulat: U grid latitude in degrees
    ulon: U grid longitude in degrees

    output:
    tlat in degrees
    tlon in degrees
    '''

    # method obtained from ice_grid.F90: subroutine Tlatlon
    ulatcos = np.cos(np.deg2rad(ulat))
    ulatsin = np.sin(np.deg2rad(ulat))

    uloncos = np.cos(np.deg2rad(ulon))
    ulonsin = np.sin(np.deg2rad(ulon))

    # initialize array with nghost=1 on each side
    nj, ni = ulatcos.shape # note: Python NetCDF is nj, ni order
    print("Tlatlon nj, ni", nj, ni)

    nghost     = 1
    workdims   = (nj+2*nghost,ni+2*nghost)
    #print("Tlatlon workdims", workdims)

    ulatcos1 = np.zeros(workdims,dtype='f8')
    ulatsin1 = np.zeros(workdims,dtype='f8')
    uloncos1 = np.zeros(workdims,dtype='f8')
    ulonsin1 = np.zeros(workdims,dtype='f8')

    # fill middle of work arrays
    ulatcos1[1:nj+1,1:ni+1] = ulatcos
    ulatsin1[1:nj+1,1:ni+1] = ulatsin

    # fill middle of work arrays
    ulatcos1[1:nj+1,1:ni+1] = ulatcos
    ulatsin1[1:nj+1,1:ni+1] = ulatsin

    uloncos1[1:nj+1,1:ni+1] = uloncos
    ulonsin1[1:nj+1,1:ni+1] = ulonsin

    # fill halos
    ulatcos1 = halo_extrapolate(ulatcos1,ew_bndy_type,ns_bndy_type)
    ulatsin1 = halo_extrapolate(ulatsin1,ew_bndy_type,ns_bndy_type)
    uloncos1 = halo_extrapolate(uloncos1,ew_bndy_type,ns_bndy_type)
    ulonsin1 = halo_extrapolate(ulonsin1,ew_bndy_type,ns_bndy_type)

    # now do computations as in ice_grid.F90:Tlatlon

    # x, y, z are full 2d
    x = uloncos1 * ulatcos1
    y = ulonsin1 * ulatcos1
    z = ulatsin1

    tx = 0.25*(x[0:nj,  0:ni  ]   + # x1
               x[0:nj,  1:ni+1]   + # x2
               x[1:nj+1,0:ni  ]   + # x3
               x[1:nj+1,1:ni+1])    # x4

    #print("Tlonlat: x.shape", x.shape)
    #print("Tlonlat: tx.shape", tx.shape)


    ty = 0.25*(y[0:nj,  0:ni  ]   + # y1
               y[0:nj,  1:ni+1]   + # y2
               y[1:nj+1,0:ni  ]   + # y3
               y[1:nj+1,1:ni+1])    # y4


    tz = 0.25*(z[0:nj,  0:ni  ]   + # z1
               z[0:nj,  1:ni+1]   + # z2
               z[1:nj+1,0:ni  ]   + # z3
               z[1:nj+1,1:ni+1])    # z4

    da = np.sqrt(tx*tx + ty*ty + tz*tz)

    tz = tz/da

    tlon = np.arctan2(ty,tx)
    tlat = np.arcsin(tz)

    # returnd tlat, tlon in degrees
    return np.rad2deg(tlat), np.rad2deg(tlon)

##########################
##########################

def get_command_line_args():
    '''
    argument parser for command line arguments
    '''

    dstr = "Interplate JRA55 data"
    parser = argparse.ArgumentParser(description=dstr)

    # add arguments
    parser.add_argument("JRADTG",  type=str, help="JRA55 file date time group")
    parser.add_argument("dstgrid", type=str, help="Destination grid file (NetCDF)")
    parser.add_argument("ncout",   type=str, help="Output file name (NetCDF)")


    # get the arguments
    args = parser.parse_args()

    # return values
    return args.JRADTG, args.dstgrid, args.ncout


################################
################################

def get_jra55_nc_dict():
    '''
    Create dictionary that links the NetCDF variable name
    with the file prefix. The file prefix is appended by 
    JRADTG from command line
    '''
    # specify dictionary with dataset prefix names 
    jra55dict = {"TPRAT_GDS4_SFC_ave3h" : "fcst_phy2m.061_tprat.reg_tl319", # precip
                 "DSWRF_GDS4_SFC_ave3h" : "fcst_phy2m.204_dswrf.reg_tl319", # downward shortwave
                 "DLWRF_GDS4_SFC_ave3h" : "fcst_phy2m.205_dlwrf.reg_tl319", # downward longwave
                 "TMP_GDS4_HTGL"        : "fcst_surf.011_tmp.reg_tl319"   , # air temp
                 "UGRD_GDS4_HTGL"       : "fcst_surf.033_ugrd.reg_tl319"  , # u velocity
                 "VGRD_GDS4_HTGL"       : "fcst_surf.034_vgrd.reg_tl319"  , # v velocity
                 "SPFH_GDS4_HTGL"       : "fcst_surf.051_spfh.reg_tl319"}   # specify humidity

    
    return jra55dict

################################
################################

def get_jra55_cice_var():
    '''
    Make dictionary relating JRA55 NetCDF variables
    to CICE variables.
    '''
    
    # specify output variable names
    # This is for current CICE expected names
    # it might be better to change CICE in long run
    cice_var = {"TPRAT_GDS4_SFC_ave3h" : "ttlpcp",
                "DSWRF_GDS4_SFC_ave3h" : "glbrad",
                "DLWRF_GDS4_SFC_ave3h" : "dlwsfc",
                "TMP_GDS4_HTGL"        : "airtmp",
                "UGRD_GDS4_HTGL"       : "wndewd",
                "VGRD_GDS4_HTGL"       : "wndnwd",
                "SPFH_GDS4_HTGL"       : "spchmd"}
    
    return cice_var

################################
################################

def init_ncout(ncout,nc1,llat,llon):

    '''
    Initialize output NetCDF file
    with proper units and dimensions.
    '''

    dsout = Dataset(ncout,'w',format='NETCDF3_64BIT_OFFSET')

    # get dimensions from size of lat
    (nlat,nlon) = llat.shape

    # create dimensions
    time  = dsout.createDimension('time',None) # unlimited
    dim_j = dsout.createDimension('dim_j',nlat)
    dim_i = dsout.createDimension('dim_i',nlon)
    
    # create time variable.
    # note is defined as 'times' (with and s) to not conflict
    # with dimension 'time'
    times          = dsout.createVariable('time','f8',('time',))
    times.units    = nc1['initial_time0_hours'].units
    times.calendar = 'gregorian'

    # loop over nc1 times
    dates = []
    dates.append(nc1['initial_time0_hours'][0] + nc1['forecast_time1'][1])
    # loop over remaining
    for h in nc1['initial_time0_hours'][1:-1]:
        for ft in nc1['forecast_time1'][:]:
            dates.append(h + ft)

    # include only first forecast_time of last initial time
    dates.append(nc1['initial_time0_hours'][-1] + nc1['forecast_time1'][0])

    # write dates to file
    times[:] = dates

    # create LON/LAT variables
    LON = dsout.createVariable('LON','f8',('dim_j','dim_i',))
    LON.units = 'degrees_east'

    LAT = dsout.createVariable('LAT','f8',('dim_j','dim_i',))
    LAT.units = 'degrees_north'

    # write LON, LAT to file
    LON[:] = llon[:,:]
    LAT[:] = llat[:,:]
    

    return dsout

################################
################################


# main subroutine
if __name__ == "__main__":

    # get jra dtg and ncout from command line
    JRADTG, dstgrid, ncout = get_command_line_args()

    # get jra55 variable/filename prefix dictionary
    jra55dict = get_jra55_nc_dict()

    # get dictionary linking jra55 variables names
    # and CICE forcing varible names
    cice_var = get_jra55_cice_var()

    # read input grid. 
    # use one of the jra55 files.
    # it is assumed all JRA data are the same grid for later
    fname = f"{jra55dict['TMP_GDS4_HTGL']:s}.{JRADTG:s}.nc"
    print("opening dataset ", fname)
    grid1_ds = Dataset(fname,'r',format='NETCDF3_64BIT_OFFSET')
    lon1     = grid1_ds['g4_lon_3'][:]  # 1D
    lat1     = grid1_ds['g4_lat_2'][:]  # 1D

    # open destination grid 
    # here it is assumed a CICE NetCDF file. 
    # the user can update as appropriate
    print("Opening ", dstgrid)
    grid2_ds = Dataset(dstgrid,'r',format='NETCDF3_64BIT_OFFSET')
    ulon2    = grid2_ds["lon"][:,:] # 2D. Assumed ULON in degrees
    ulat2    = grid2_ds["lat"][:,:] # 2D. Assumed ULAT in degrees
    if np.max(np.abs(ulat2)) < 10. :
       ulon2 = np.rad2deg(ulon2)
       ulat2 = np.rad2deg(ulat2)

    # make tgrid from ugrid
    ew_bndy_type = 'cyclic'
    ns_bndy_type = 'open'
    tlat2, tlon2 = Tlatlon(ulat2,ulon2,ew_bndy_type,ns_bndy_type)

    # make regridders 
    print("making bilinear regridder")
    method = 'bilinear'
    periodic = True
    blin_grid_name = 'bilinear_jra55_gx3.nc'
    rgrd_bilinear = make_regridder(lon1,lat1,tlon2,tlat2,
                                   method,periodic,blin_grid_name)

    # setup output dataset by adding lat/lon
    dsout = init_ncout(ncout,grid1_ds,tlat2,tlon2)
    
    # no longer need grid1, grid2
    grid1_ds.close()
    grid2_ds.close()

    # do the regridding
    # Loop over all the files using regridder from above
    # and add to dataout
    for var, fprefix in jra55dict.items():
        fname = f"{fprefix:s}.{JRADTG:s}.nc"
        print("reading ", fname)
        jra_ds = Dataset(fname,'r',format='NETCDF3_CLASSIC')

        # make output variable
        data = dsout.createVariable(cice_var[var],'f4',('time','dim_j','dim_i'))
        
        # do interpolation
        print("Interpolating ", var)

        if var.find('ave3h') > 0: # ave3r in var
            # use bilinear here
            d = rgrd_bilinear(jra_ds[var][:,:,:,:])

            # write to file in correct time order
            for t in range(d.shape[0]):
                for n in range(d.shape[1]):
                    #print('indx (2*t)+n = ', (2*t)+n)
                    data[(2*t)+n,:,:] = d[t,n,:,:]

        else:
            # instantaneous use bilinear
            d = rgrd_bilinear(jra_ds[var][:,:,:,:])
            
            # write to file in correct time order. 
            # note need to write 2nd forecast_time first.
            # in this case first forecast_time is NaN
            data[0,:,:] = d[0,1,:,:]
            for t in range(1,d.shape[0]-1):
                for n in range(d.shape[1]):
                    #print('indx (2*t)+n-1 = ', (2*t)+n-1)
                    data[(2*t)+n-1,:,:] = d[t,n,:,:]

            # write first forecast time of last initial time
            # second forecast time is NAN
            data[-1,:,:] = d[-1,0,:,:]

        # add coordinates attribute
        data.coordinates = "LON LAT time"
        data.long_name   = jra_ds[var].long_name
        data.units       = jra_ds[var].units
 
        precip_factor = 1. / 86400.

        # Convert mm / day to kg/m^2/s.
        if var.find('PRAT') > 0:
            data[:] = data[:] * precip_factor
            data.units = 'kg/m2/s'
        else:
            data.units       = jra_ds[var].units
        
        # close jra55 file
        jra_ds.close()

    
    # write tou output file
    # close output file
    dsout.close()
        
    print("Done")
    
