#!/usr/bin/env python

# This script performs the t-test validation for non-bit-for-bit results for the 
# CICE model.

# Written by: Matthew Turner
# Date: October, 2017

import netCDF4 as nc
import os
import sys
import numpy as np
import numpy.ma as ma
import logging

def maenumerate(marr):
    '''
    This function provides the enumerate functionality for masked arrays
    '''

    mask = ~marr.mask.ravel()
    try:   # Python 2
        import itertools
        for i,m in itertools.izip(np.ndindex(marr.shape[-2:]),mask):
            if m: yield i
    except:  # Python 3
        for i,m in zip(np.ndindex(marr.shape[-2:]),mask):
            if m: yield i

def read_data(base_dir,test_dir,var):
    '''
    Read the baseline and test data for sea ice thickness.  The calculate
    the difference for all locations where sea ice thickness is greater
    than 0.01 meters.
    '''
    # The path to output files for simulation 'a' (the '-bc' simulation)
    if base_dir.endswith(('history','history/')):
        path_a = base_dir
    else:
        path_a = base_dir + '/history/'
    
    # The path to output files for simulation 'b' (the test simulation)
    if test_dir.endswith(('history','history/')):
        path_b = test_dir
    else:
        path_b = test_dir + '/history/'
    
    # Find the number of output files to be read in
    files_a = [i for i in os.listdir(path_a+'/') if i.startswith('iceh_inst.')]
    files_b = [i for i in os.listdir(path_b+'/') if i.startswith('iceh_inst.')]
    
    if not len(files_a) == len(files_b):
        logger.error("Number of output files for baseline simulation does not match the number" + \
              " of files for the test simulation.  Exiting...\n" + \
              "Baseline directory: {}\n".format(path_a) + \
              "   # of files: {}\n".format(len(files_a)) + \
              "Test directory: {}\n".format(path_b) + \
              "   # of files: {}".format(len(files_b)))
        sys.exit(-1)
      
    num_files = len(files_a)
    logger.info("Number of files: {}".format(num_files))
    
    # Get the array dimensions from the file
    nfid = nc.Dataset("{}/{}".format(path_a,files_a[0]),'r')
    ni = nfid.dimensions['ni'].size
    nj = nfid.dimensions['nj'].size
    nfid.close()
    
    # Pre-allocate a numpy array with a "time" dimension
    data_a = np.zeros((num_files,nj,ni))
    data_b = np.zeros((num_files,nj,ni))
    
    # Read in the data 
    cnt = 0
    for fname in sorted(files_a):
        nfid = nc.Dataset("{}/{}".format(path_a,fname),'r')
        fill_value_a = nfid.variables[var]._FillValue
        data_a[cnt,:,:] = nfid.variables[var][:]
        cnt += 1
        nfid.close()
    data_a[data_a==fill_value_a] = 0.0
    
    cnt = 0
    for fname in sorted(files_b):
        nfid = nc.Dataset("{}/{}".format(path_b,fname),'r')
        fill_value_b = nfid.variables[var]._FillValue
        data_b[cnt,:,:] = nfid.variables[var][:]
        cnt += 1
        nfid.close()
    data_b[data_b==fill_value_b] = 0.0

    # Calculate the difference and mask the points where the the difference at
    #   every timestep is 0, or the sea ice thickness for every timestep is < 0.01 meters
    #   for data_a or data_b
    data_d = data_a - data_b
    mask_d = np.logical_or(\
                  np.logical_or(\
                       np.all(np.equal(data_d,0.),axis=0),np.all(data_a < 0.01,axis=0))\
                  ,np.all(data_b < 0.01,axis=0))
    mask_array_a = np.zeros_like(data_d)
    mask_array_b = np.zeros_like(data_d)
    mask_array_d = np.zeros_like(data_d)
    for x,value in np.ndenumerate(mask_d):
        i,j = x
        mask_array_d[:,i,j] = value
        mask_array_a[:,i,j] = value
        mask_array_b[:,i,j] = value
    data_d = ma.masked_array(data_d,mask=mask_array_d)
    data_a = ma.masked_array(data_a,mask=mask_array_a)
    data_b = ma.masked_array(data_b,mask=mask_array_b)

    return data_a, data_b, data_d, num_files, path_a, fname

def two_stage_test(data_a,data_b,num_files,data_d):
    # Calculate the mean of the difference
    mean_d = np.mean(data_d,axis=0)
    variance_d = np.sum(np.power(data_d - mean_d,2)) / (num_files - 1)
    
    # Calculate the mean from 1:end-1 and 2:end
    mean_nm1_d = np.mean(data_d[:-1,:,:],axis=0)
    mean_2n_d = np.mean(data_d[1:,:,:],axis=0)
    
    # Calculate equation (5) for both simulations
    r1_num = np.zeros_like(mean_d)
    r1_den1 = np.zeros_like(mean_d)
    r1_den2 = np.zeros_like(mean_d)
    for i in np.arange(np.size(data_a,axis=0)-1):
        r1_num = r1_num + (data_d[i,:,:]-mean_nm1_d[:,:])*(data_d[i+1,:,:]-mean_2n_d[:,:])
        r1_den1 = r1_den1 + np.power(data_d[i,:,:]-mean_nm1_d[:,:],2)
    
    for i in np.arange(1,np.size(data_a,axis=0)):
        r1_den2 = r1_den2 + np.power(data_d[i,:,:] - mean_2n_d[:,:],2)
    
    r1 = r1_num / np.sqrt(r1_den1*r1_den2)
    
    # Calculate the effective sample size
    n_eff = num_files*( (1.-r1) / (1.+r1) )
    n_eff[n_eff < 2] = 2
    n_eff[n_eff > num_files] = num_files
    
    # Calculate the t-statistic with n_eff
    t_val = mean_d / np.sqrt(variance_d / n_eff)
    
    # Effective degrees of freedom
    df = n_eff - 1
    
    # Read in t_crit table
    nfid = nc.Dataset("configuration/scripts/tests/QC/CICE_t_critical_p0.8.nc",'r')
    df_table = nfid.variables['df'][:]
    t_crit_table = nfid.variables['tcrit'][:]
    nfid.close()
    t_crit = np.zeros_like(t_val)

    # Calculate critical t-value for each grid cell, based on the t_crit table
    for x in maenumerate(data_d):
        min_val = np.min(np.abs(df[x]-df_table))
        idx = np.where(np.abs(df[x]-df_table)==min_val)
        t_crit[x] = t_crit_table[idx]
    
    # Create an array of Pass / Fail values for each grid cell
    passed_array = abs(t_val) > t_crit

    # Initialize variable that defines whether or not to calculate the area-weighted
    #   fraction of failures to false.
    calc_fval = False

    if np.any(abs(t_val) > t_crit):
        # A non-zero number of grid cells failed
        passed = False
    
    elif np.all(abs(t_val)<=t_crit) and np.all(n_eff >= 30):
        # Every grid cell passed, and has an effective sample size >= 30
        passed = True
    
    elif np.all(abs(t_val) <= t_crit) and np.any(n_eff < 30):
        # Calculate the T-statistic using actual sample size
        t_val = mean_d / np.sqrt(variance_d / num_files)
    
        # Find t_crit from the nearest value on the Lookup Table Test
        nfid = nc.Dataset("configuration/scripts/tests/QC/CICE_Lookup_Table_p0.8_n1825.nc",'r')
        r1_table = nfid.variables['r1'][:]
        t_crit_table = nfid.variables['tcrit'][:]
        nfid.close()
        
        # Fill t_crit based on lookup table
        for x in maenumerate(data_d):
            min_val = np.min(np.abs(r1[x]-r1_table))
            idx = np.where(np.abs(r1[x]-r1_table)==min_val)
            t_crit[x] = t_crit_table[idx]
    
        # Create an array showing locations of Pass / Fail grid cells
        passed_array = abs(t_val) > t_crit

        if np.any(abs(t_val) > t_crit):
            # A non-zero number of grid cells has failed
            passed = False

        elif np.all(abs(t_val) <= t_crit):
            # All grid cells passed the t-test
            passed = True
            # Set variable to make the script calculate the area-weighted fraction
            #   of failure grid cells
            calc_fval = True

        else:
            logger.error('TWO-STAGE TEST NOT CONCLUSIVE')
            passed = False
    
    else:
        logger.error('TEST NOT CONCLUSIVE')
        passed = False

    # Calculate the area-weighted fraction of the test region that failed (f_val).
    #   If f_val is greater than or equal to the critical fraction, the test fails"
    if calc_fval:
        f_val = critical_fraction(data_a,passed_array)
        f_crit = 0.5
        if f_val >= f_crit:
            passed = False
            logger.info('Area-weighted fraction of failures is greater than ' + \
                        'critical fraction.  Test failed.')
            logger.info('Area-weighted fraction of failures = {}'.format(f_val))
            return passed, passed_array

    if passed:
        logger.info('Two-Stage Test Passed')
    else:
        logger.info('Two-Stage Test Failed')

    return passed, passed_array

def critical_fraction(data_a,passed_array):
    # First calculate the weight attributed to each grid point (based on Area)
    nfid = nc.Dataset("{}/{}".format(path_a,fname),'r')
    tarea = nfid.variables['tarea'][:]
    nfid.close()
    tarea = ma.masked_array(tarea,mask=data_a[0,:,:].mask)
    area_weight = tarea / np.sum(tarea)

    # Calculate the area weight of the failing grid cells
    weight_tot = 0
    weight_fail = 0
    for x in maenumerate(data_a):
        weight_tot += area_weight[x]
        if passed_array[x]:
            weight_fail += area_weight[x]

    return weight_fail/weight_tot

def skill_test(path_a,fname,data_a,data_b,num_files,tlat,hemisphere):
    '''Calculate Taylor Skill Score'''
    # First calculate the weight attributed to each grid point (based on Area)
    nfid = nc.Dataset("{}/{}".format(path_a,fname),'r')
    tarea = nfid.variables['tarea'][:]
    nfid.close()
    tarea = ma.masked_array(tarea,mask=data_a[0,:,:].mask)
    area_weight = tarea / np.sum(tarea)

    weighted_mean_a = 0
    weighted_mean_b = 0
    for i in np.arange(num_files):
        weighted_mean_a = weighted_mean_a + np.sum(area_weight*data_a[i,:,:])
        weighted_mean_b = weighted_mean_b + np.sum(area_weight*data_b[i,:,:])
    
    weighted_mean_a = weighted_mean_a / num_files
    weighted_mean_b = weighted_mean_b / num_files
    
    nonzero_weights = np.count_nonzero(area_weight)
    area_var_a = 0
    area_var_b = 0
    for t in np.arange(num_files):
        area_var_a = area_var_a + np.sum(area_weight*np.power(data_a[t,:,:]-weighted_mean_a,2))
        area_var_b = area_var_b + np.sum(area_weight*np.power(data_b[t,:,:]-weighted_mean_b,2))
    
    area_var_a = nonzero_weights / (num_files * nonzero_weights - 1.) * area_var_a
    area_var_b = nonzero_weights / (num_files * nonzero_weights - 1.) * area_var_b
    std_a = np.sqrt(area_var_a)
    std_b = np.sqrt(area_var_b)
    
    combined_cov = 0
    for i in np.arange(num_files):
        combined_cov = combined_cov + np.sum(area_weight*(data_a[i,:,:]-weighted_mean_a)*\
                                            (data_b[i,:,:]-weighted_mean_b))
    
    combined_cov = nonzero_weights / (num_files * nonzero_weights - 1.) * combined_cov
    
    weighted_r = combined_cov / (std_a*std_b)
    
    s = np.power((1+weighted_r)*(std_a*std_b)/\
                 (area_var_a + area_var_b),2)
    
    logger.debug('s = {}'.format(s))
    
    s_crit = 0.99
    if s<0 or s>1:
        logger.error('Skill score out of range for {} Hemisphere'.format(hemisphere))
        passed = False
    elif s > s_crit:
        logger.info('Quadratic Skill Test Passed for {} Hemisphere'.format(hemisphere))
        passed = True
    else:
        logger.info('Quadratic Skill Test Failed for {} Hemisphere'.format(hemisphere))
        passed = False
    return passed

def plot_data(data,lat,lon,units):
    '''This function plots CICE data and creates a .png file (ice_thickness_map.png).'''

    try:
        logger.info('Creating map of the data (ice_thickness_map.png)')
        # Load the necessary plotting libraries
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        from matplotlib.colors import LinearSegmentedColormap
    except:
        logger.warning('Error loading necessary Python modules in plot_data function')

    # Suppress Matplotlib deprecation warnings
    import warnings
    warnings.filterwarnings("ignore",category=UserWarning)

    # Create the figure and axis
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_axes([0.05,0.08,0.9,0.9])
 
    # Create the basemap, and draw boundaries
    m = Basemap(projection='kav7',lon_0=180.,resolution='l')
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines()
    m.drawcountries()

    # Plot the data as a scatter plot
    x,y=m(lon,lat)
    sc = m.scatter(x,y,c=data,cmap='jet',lw=0)
 
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1],fontsize=10)
    m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0],fontsize=10)

    plt.title('CICE Ice Thickness')

    # Create the colorbar and add Pass / Fail labels
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom",size="5%",pad=0.5)
    cb = plt.colorbar(sc,cax=cax,orientation="horizontal",format="%.2f")
    cb.set_label(units,x=1.0)

    plt.savefig('ice_thickness_map.png',dpi=300)

def plot_two_stage_failures(data,lat,lon):
    '''This function plots each grid cell and whether or not it Passed or Failed 
       the two-stage test.  It then either creates a .png file 
       (two_stage_test_failure_map.png), or saves the failure locations to a 
       text file.
    '''

    # Convert the boolean array (data) to an integer array
    int_data = data.astype(int)

    try: 
        logger.info('Creating map of the failures (two_stage_test_failure_map.png)')
        # Load the necessary plotting libraries
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        from matplotlib.colors import LinearSegmentedColormap

        # Suppress Matplotlib deprecation warnings
        import warnings
        warnings.filterwarnings("ignore",category=UserWarning)

        # Create the figure and axis
        fig = plt.figure(figsize=(12,8))
        ax = fig.add_axes([0.05,0.08,0.9,0.9])
 
        # Create the basemap, and draw boundaries
        m = Basemap(projection='kav7',lon_0=180.,resolution='l')
        m.drawmapboundary(fill_color='white')
        m.drawcoastlines()
        m.drawcountries()

        # Create the custom colormap
        colors = [(0,0,1),(1,0,0)]  # Blue, Red
        cmap_name = 'RB_2bins'
        cm = LinearSegmentedColormap.from_list(cmap_name,colors,N=2)

        # Plot the data as a scatter plot
        x,y=m(lon,lat)
        sc = m.scatter(x,y,c=int_data,cmap=cm,lw=0,vmin=0,vmax=1)
 
        m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1],fontsize=10)
        m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0],fontsize=10)

        plt.title('CICE Two-Stage Test Failures')

        # Create the colorbar and add Pass / Fail labels
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom",size="5%",pad=0.5)
        cb = plt.colorbar(sc,cax=cax,orientation="horizontal",format="%.0f")
        cb.set_ticks([])
        cb.ax.text(-0.01,-0.5,'PASS')
        cb.ax.text(0.99,-0.5,'FAIL')

        plt.savefig('two_stage_test_failure_map.png',dpi=300)
    except:
        logger.warning('')
        logger.warning('Unable to plot the data.  Saving latitude and longitude')
        logger.warning('for ONLY failures to two_stage_failure_locations.txt')

        # Create a file and write the failures only to the file
        f = open('two_stage_failure_locations.txt','w')
        f.write('# CICE Two-stage test failures\n')
        f.write('# Longitude,Latitude\n')
        for i in range(data.shape[0]):
          for j in range(data.shape[1]):
            if (not data.mask[i,j]) and data[i,j]:
              f.write('{},{}\n'.format(tlon[i,j],tlat[i,j]))

        f.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script performs the T-test for \
                           CICE simulations that should be bit-for-bit, but are not.')
    parser.add_argument('base_dir',help='Path to the baseline history (iceh_inst*) files.  REQUIRED')
    parser.add_argument('test_dir',help='Path to the test history (iceh_inst*) files.  REQUIRED')
    parser.add_argument('-v','--verbose',dest='verbose',help='Print debug output?', \
                        action='store_true')
    
    parser.set_defaults(verbose=False)
    
    # If no arguments are provided, print the help message
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    
    # Set up the logger
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    # Log to log file as well as stdout
    fh = logging.FileHandler(r'qc_log.txt','w')
    logger = logging.getLogger(__name__)
    logger.addHandler(fh)
    
    var = 'hi'
    data_a, data_b, data_d, num_files, path_a, fname = read_data(args.base_dir,args.test_dir,var)

    if np.ma.all(data_d.mask):
        logger.info("Data is bit-for-bit.  No need to run QC test")
        sys.exit(0)

    # Run the two-stage test
    passed,passed_array = two_stage_test(data_a,data_b,num_files,data_d)
    
    # Get the latitude and longitude information for the domain
    nfid = nc.Dataset("{}/{}".format(path_a,fname),'r')
    tlat = nfid.variables['TLAT'][:]
    tlon = nfid.variables['TLON'][:]
    nfid.close()

    # If test failed, attempt to create a plot of the failure locations
    if not passed:
        plot_two_stage_failures(passed_array,tlat,tlon)

    # Create a northern hemisphere and southern hemisphere mask
    mask_tlat = tlat < 0
    mask_nh = np.zeros_like(data_a)
    mask_sh = np.zeros_like(data_a)
    for (i,j),value in np.ndenumerate(mask_tlat):
        mask_nh[:,i,j] = value
        mask_sh[:,i,j] = not value
    
    # Run skill test on northern hemisphere
    data_nh_a = ma.masked_array(data_a,mask=mask_nh)
    data_nh_b = ma.masked_array(data_b,mask=mask_nh)
    if np.ma.all(data_nh_a.mask) and np.ma.all(data_nh_b.mask):
        logger.info("Northern Hemisphere data is bit-for-bit")
        passed_nh = True
    else:
        passed_nh = skill_test(path_a,fname,data_nh_a,data_nh_b,num_files,tlat,'Northern')
    
    # Run skill test on southern hemisphere
    data_sh_a = ma.masked_array(data_a,mask=mask_sh)
    data_sh_b = ma.masked_array(data_b,mask=mask_sh)
    if np.ma.all(data_sh_a.mask) and np.ma.all(data_sh_b.mask):
        logger.info("Southern Hemisphere data is bit-for-bit")
        passed_sh = True
    else:
        passed_sh = skill_test(path_a,fname,data_sh_a,data_sh_b,num_files,tlat,'Southern')
    
    passed_skill = passed_nh and passed_sh
    
    logger.info('')
    if not passed and not passed_skill:
        logger.error('Quality Control Test FAILED')
        sys.exit(1)  # exit with an error return code
    else:
        logger.info('Quality Control Test PASSED')
        sys.exit(0)  # exit with successfull return code
    
