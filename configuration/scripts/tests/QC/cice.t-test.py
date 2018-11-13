#!/usr/bin/env python
'''
This script performs the t-test validation for non-bit-for-bit results for the
CICE model.

Written by: Matthew Turner
Date: October, 2017
'''

import os
import sys
import logging
import numpy as np
import numpy.ma as ma
import netCDF4 as nc

def maenumerate(marr):
    '''
    This function provides the enumerate functionality for masked arrays
    '''
    mask = ~marr.mask.ravel()
    try:   # Python 2
        import itertools
        for i, m in itertools.izip(np.ndindex(marr.shape[-2:]), mask):
            if m: yield i
    except:  # Python 3
        for i, m in zip(np.ndindex(marr.shape[-2:]), mask):
            if m: yield i

def gen_filenames(base_dir, test_dir):
    '''
    This function is passed the directories of the baseline and test history
    files and generates a list of filenames for each.
    '''
    # The path to output files for simulation 'a' (the '-bc' simulation)
    if base_dir.endswith(('history', 'history/')):
        path_a = base_dir
    else:
        path_a = base_dir + '/history/'

    # The path to output files for simulation 'b' (the test simulation)
    if test_dir.endswith(('history', 'history/')):
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

    logger.info("Number of files: %d", len(files_a))

    return path_a, path_b, files_a, files_b

def get_geom(path, file):
    '''
    This function reads the ni, nj, tlat, and tlon variables from a netcdf file
    '''
    fid = nc.Dataset("{}/{}".format(path, file), 'r')
    tlat = fid.variables['TLAT'][:]
    tlon = fid.variables['TLON'][:]
    ni = fid.dimensions['ni'].size
    nj = fid.dimensions['nj'].size
    fid.close()

    return ni, nj, tlat, tlon

def read_data(path_a, path_b, files_a, files_b, ni, nj):
    '''
    Read the baseline and test data for sea ice thickness.  The calculate
    the difference for all locations where sea ice thickness is greater
    than 0.01 meters.
    '''
    def fill_data_array(path, files, nj, ni):
        '''Function to fill the data arrays'''
        # Initialize the data array
        data = np.zeros((len(files), nj, ni),dtype=np.float32)
        # Read in the data
        logger.debug('Reading in data for files in %s', path)
        cnt = 0
        for fname in sorted(files):
            nfid = nc.Dataset("{}/{}".format(path, fname), 'r')
            fill_value = nfid.variables[var]._FillValue
            data[cnt, :, :] = nfid.variables[var][:]
            cnt += 1
            nfid.close()
        data[data == fill_value] = 0.0

        return data

    def calc_diff(data_a, data_b):
        '''
        Calculate the difference and mask the points where the the difference at
        every timestep is 0, or the sea ice thickness for every timestep is < 0.01 meters
        for data_a or data_b
        '''
        data_d = data_a - data_b
        mask_d = np.logical_or(\
                      np.logical_or(\
                           np.all(np.equal(data_d, 0.), axis=0), np.all(data_a < 0.01, axis=0))\
                      , np.all(data_b < 0.01, axis=0))
        mask_array_a = np.zeros_like(data_d)
        for x, value in np.ndenumerate(mask_d):
            i, j = x
            mask_array_a[:, i, j] = value
        data_a = ma.masked_array(data_a, mask=mask_array_a)
        data_b = ma.masked_array(data_b, mask=mask_array_a)
        data_d = ma.masked_array(data_d, mask=mask_array_a)
        del mask_array_a

        return data_a, data_b, data_d

    var = 'hi'

    data_a = fill_data_array(path_a, files_a, nj, ni)
    data_b = fill_data_array(path_b, files_b, nj, ni)

    data_a, data_b, data_d = calc_diff(data_a, data_b)

    return data_a, data_b, data_d

def two_stage_test(data_a, num_files, data_d, fname, path):
    '''
    This function performs the Two-Stage Paired Thickness Test
    '''
    def stage_one(data_d, num_files, mean_d, variance_d):
        logger.debug('Running step 1 of 2-stage test')

        # Calculate the mean from 1:end-1 and 2:end
        mean_nm1_d = np.mean(data_d[:-1, :, :], axis=0)
        mean_2n_d = np.mean(data_d[1:, :, :], axis=0)

        # Calculate equation (5) for both simulations
        r1_num = np.zeros_like(mean_d)
        r1_den1 = np.zeros_like(mean_d)
        r1_den2 = np.zeros_like(mean_d)
        for i in np.arange(np.size(data_a, axis=0)-1):
            r1_num = r1_num + (data_d[i, :, :]-mean_nm1_d[:, :])*(data_d[i+1, :, :]-mean_2n_d[:, :])
            r1_den1 = r1_den1 + np.power(data_d[i, :, :]-mean_nm1_d[:, :], 2)

        for i in np.arange(1, np.size(data_a, axis=0)):
            r1_den2 = r1_den2 + np.power(data_d[i, :, :] - mean_2n_d[:, :], 2)

        r1 = r1_num / np.sqrt(r1_den1*r1_den2)

        # Calculate the effective sample size
        n_eff = num_files * ((1.-r1) / (1.+r1))
        n_eff[n_eff < 2] = 2
        n_eff[n_eff > num_files] = num_files

        # Calculate the t-statistic with n_eff
        t_val = mean_d / np.sqrt(variance_d / n_eff)

        # Effective degrees of freedom
        df = n_eff - 1

        # Read in t_crit table
        nfid = nc.Dataset("configuration/scripts/tests/QC/CICE_t_critical_p0.8.nc", 'r')
        df_table = nfid.variables['df'][:]
        t_crit_table = nfid.variables['tcrit'][:]
        nfid.close()
        t_crit = np.zeros_like(t_val)

        # Calculate critical t-value for each grid cell, based on the t_crit table
        for x in maenumerate(data_d):
            min_val = np.min(np.abs(df[x]-df_table))
            idx = np.where(np.abs(df[x]-df_table) == min_val)
            t_crit[x] = t_crit_table[idx]

        # Create an array of Pass / Fail values for each grid cell
        H1 = np.abs(t_val) > t_crit

        return n_eff, H1, r1, t_crit

    # Calculate the mean of the difference
    mean_d = np.mean(data_d, axis=0)
    variance_d = np.sum(np.power(data_d - mean_d, 2)) / (num_files - 1)

    n_eff, H1, r1, t_crit = stage_one(data_d, num_files, mean_d, variance_d)

    if np.all(H1 == False) and np.all(n_eff >= 30):
        # H0 confirmed in all cells, and all effective sample size >= 30
        logger.debug('H0 confirmed in all cells, and all effective sample size >= 30')
        logger.info('2 stage test passed')
        return True, H1

    elif np.all(H1):
        # H1 in every grid cell
        logger.debug('H1 in all cells')
        logger.info('2 stage test failed')
        return False, H1

########### H0 confirmed for some grid cells with n_eff < 30 ############
    logger.debug('Number of H1 grid cells after stage 1 = %d', np.sum(H1))

    logger.debug('Running step 2 of 2-stage test')

    # Find the indices where n_eff is less than 30, and H0 is confirmed
    tmp_idx = np.where(n_eff < 30) and np.where(H1 == False)

    # Calculate the T-statistic using actual sample size
    t_val = mean_d / np.sqrt(variance_d / num_files)

    # Find t_crit from the nearest value on the Lookup Table Test
    nfid = nc.Dataset("configuration/scripts/tests/QC/CICE_Lookup_Table_p0.8_n1825.nc", 'r')
    r1_table = nfid.variables['r1'][:]
    t_crit_table = nfid.variables['tcrit'][:]
    nfid.close()

    # Fill t_crit based on lookup table
    for x in maenumerate(data_d):
        min_val = np.min(np.abs(r1[x]-r1_table))
        idx = np.where(np.abs(r1[x]-r1_table) == min_val)
        t_crit[x] = t_crit_table[idx]

    # Create an array showing locations of Pass / Fail grid cells
    H1[tmp_idx] = abs(t_val[tmp_idx]) > t_crit[tmp_idx]

    logger.debug('Number of H1 grid cells after stage 2 = %d', np.sum(H1))

    if np.all(H1):
        # H1 in all grid cells
        logger.debug('H1 in all cells, stage 2')
        logger.info('2 Stage Test Failed')
        return False, H1

    elif np.all(H1 == False):
        # H0 confirmed in all grid cells
        logger.debug('H0 confirmed in all cells with n_eff < 30')
        logger.info('2 Stage Test Passed')
        return True, H1

####### Some grid cells have H0 confirmed, and some do not ######

    # Calculate the area-weighted fraction of the test region that failed (f_val).
    #   If f_val is greater than or equal to the critical fraction, the test fails"
    f_val = critical_fraction(data_a, H1, fname, path)
    f_crit = 0.5
    if f_val >= f_crit:
        logger.info('2 Stage Test Failed')
        logger.debug('Area-weighted fraction of failures is greater than ' + \
                    'critical fraction.  Test failed.')
        logger.debug('Area-weighted fraction of failures = %f', f_val)
        return False, H1
    else:
        logger.info('2 Stage Test Passed')
        logger.debug('Area-weighted fraction of failures = %f', f_val)
        return True, H1

def critical_fraction(data_a, failures, fname, path_a):
    '''
    This function calculates the area-weighted average of cells where H1 is true.
    '''
    logger.debug('Calculating area-weighted average of H1 cells')
    # First calculate the weight attributed to each grid point (based on Area)
    nfid = nc.Dataset("{}/{}".format(path_a, fname), 'r')
    tarea = nfid.variables['tarea'][:]
    nfid.close()
    tarea = ma.masked_array(tarea, mask=data_a[0, :, :].mask)
    area_weight = tarea / np.sum(tarea)

    # Calculate the area weight of the failing grid cells
    weight_tot = 0
    weight_fail = 0
    for x in maenumerate(data_a):
        weight_tot += area_weight[x]
        if failures[x]:
            weight_fail += area_weight[x]

    return weight_fail/weight_tot

def skill_test(path_a, fname, data_a, data_b, num_files, hemisphere):
    '''Calculate Taylor Skill Score'''
    # First calculate the weight attributed to each grid point (based on Area)
    nfid = nc.Dataset("{}/{}".format(path_a, fname), 'r')
    tarea = nfid.variables['tarea'][:]
    nfid.close()
    tarea = ma.masked_array(tarea, mask=data_a[0, :, :].mask)
    area_weight = tarea / np.sum(tarea)

    weighted_mean_a = 0
    weighted_mean_b = 0
    for i in np.arange(num_files):
        weighted_mean_a = weighted_mean_a + np.sum(area_weight*data_a[i, :, :])
        weighted_mean_b = weighted_mean_b + np.sum(area_weight*data_b[i, :, :])

    weighted_mean_a = weighted_mean_a / num_files
    weighted_mean_b = weighted_mean_b / num_files

    nonzero_weights = np.count_nonzero(area_weight)
    area_var_a = 0
    area_var_b = 0
    for t in np.arange(num_files):
        area_var_a = area_var_a + np.sum(area_weight*np.power(data_a[t, :, :]-weighted_mean_a, 2))
        area_var_b = area_var_b + np.sum(area_weight*np.power(data_b[t, :, :]-weighted_mean_b, 2))

    area_var_a = nonzero_weights / (num_files * nonzero_weights - 1.) * area_var_a
    area_var_b = nonzero_weights / (num_files * nonzero_weights - 1.) * area_var_b
    std_a = np.sqrt(area_var_a)
    std_b = np.sqrt(area_var_b)

    combined_cov = 0
    for i in np.arange(num_files):
        combined_cov = combined_cov + np.sum(area_weight*(data_a[i, :, :]-weighted_mean_a)*\
                                            (data_b[i, :, :]-weighted_mean_b))

    combined_cov = nonzero_weights / (num_files * nonzero_weights - 1.) * combined_cov

    weighted_r = combined_cov / (std_a*std_b)

    s = np.power((1+weighted_r)*(std_a*std_b)/\
                 (area_var_a + area_var_b), 2)

    logger.debug('%s Hemisphere skill score = %f', hemisphere, s)

    s_crit = 0.99
    if s < 0 or s > 1:
        logger.error('Skill score out of range for %s Hemisphere', hemisphere)
        return False
    elif s > s_crit:
        logger.info('Quadratic Skill Test Passed for %s Hemisphere', hemisphere)
        return True
    else:
        logger.info('Quadratic Skill Test Failed for %s Hemisphere', hemisphere)
        return False

def plot_data(data, lat, lon, units):
    '''This function plots CICE data and creates a .png file (ice_thickness_map.png).'''

    try:
        logger.info('Creating map of the data (ice_thickness_map.png)')
        # Load the necessary plotting libraries
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        from mpl_toolkits.axes_grid1 import make_axes_locatable
    except ImportError:
        logger.warning('Error loading necessary Python modules in plot_data function')

    # Suppress Matplotlib deprecation warnings
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    # Create the figure and axis
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_axes([0.05, 0.08, 0.9, 0.9])

    # Create the basemap, and draw boundaries
    m = Basemap(projection='kav7', lon_0=180., resolution='l')
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines()
    m.drawcountries()

    # Plot the data as a scatter plot
    x, y = m(lon, lat)
    sc = m.scatter(x, y, c=data, cmap='jet', lw=0)

    m.drawmeridians(np.arange(0, 360, 60), labels=[0, 0, 0, 1], fontsize=10)
    m.drawparallels(np.arange(-90, 90, 30), labels=[1, 0, 0, 0], fontsize=10)

    plt.title('CICE Ice Thickness')

    # Create the colorbar and add Pass / Fail labels
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.5)
    cb = plt.colorbar(sc, cax=cax, orientation="horizontal", format="%.2f")
    cb.set_label(units, x=1.0)

    plt.savefig('ice_thickness_map.png', dpi=300)

def plot_two_stage_failures(data, lat, lon):
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
        warnings.filterwarnings("ignore", category=UserWarning)

        # Create the figure and axis
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_axes([0.05, 0.08, 0.9, 0.9])

        # Create the basemap, and draw boundaries
        m = Basemap(projection='kav7', lon_0=180., resolution='l')
        m.drawmapboundary(fill_color='white')
        m.drawcoastlines()
        m.drawcountries()

        # Create the custom colormap
        colors = [(0, 0, 1), (1, 0, 0)]  # Blue, Red
        cmap_name = 'RB_2bins'
        cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=2)

        # Plot the data as a scatter plot
        x, y = m(lon, lat)
        sc = m.scatter(x, y, c=int_data, cmap=cm, lw=0, vmin=0, vmax=1)

        m.drawmeridians(np.arange(0, 360, 60), labels=[0, 0, 0, 1], fontsize=10)
        m.drawparallels(np.arange(-90, 90, 30), labels=[1, 0, 0, 0], fontsize=10)

        plt.title('CICE Two-Stage Test Failures')

        # Create the colorbar and add Pass / Fail labels
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        cb = plt.colorbar(sc, cax=cax, orientation="horizontal", format="%.0f")
        cb.set_ticks([])
        cb.ax.text(-0.01, -0.5, 'PASS')
        cb.ax.text(0.99, -0.5, 'FAIL')

        plt.savefig('two_stage_test_failure_map.png', dpi=300)
    except:
        logger.warning('')
        logger.warning('Unable to plot the data.  Saving latitude and longitude')
        logger.warning('for ONLY failures to two_stage_failure_locations.txt')

        # Create a file and write the failures only to the file
        f = open('two_stage_failure_locations.txt', 'w')
        f.write('# CICE Two-stage test failures\n')
        f.write('# Longitude,Latitude\n')
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if (not data.mask[i, j]) and data[i, j]:
                    f.write('{},{}\n'.format(lon[i, j], lat[i, j]))

        f.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='This script performs the T-test for \
                           CICE simulations that should be bit-for-bit, but are not.')
    parser.add_argument('base_dir', \
                help='Path to the baseline history (iceh_inst*) files.  REQUIRED')
    parser.add_argument('test_dir', \
                help='Path to the test history (iceh_inst*) files.  REQUIRED')
    parser.add_argument('-v', '--verbose', dest='verbose', help='Print debug output?', \
                        action='store_true')

    parser.set_defaults(verbose=False)

    # If no arguments are provided, print the help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Set up the logger
    global logger
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    # Log to log file as well as stdout
    fh = logging.FileHandler(r'qc_log.txt', 'w')
    logger = logging.getLogger(__name__)
    logger.addHandler(fh)

    logger.info('Running QC test on the following directories:')
    logger.info('  {}'.format(args.base_dir))
    logger.info('  {}'.format(args.test_dir))

    dir_a, dir_b, files_base, files_test = gen_filenames(args.base_dir, args.test_dir)

    nfiles = len(files_base)

    nlon, nlat, t_lat, t_lon = get_geom(dir_a, files_base[0])

    data_base, data_test, data_diff = read_data(dir_a, dir_b, files_base, files_test, nlon, nlat)

    if np.ma.all(data_diff.mask):
        logger.info("Data is bit-for-bit.  No need to run QC test")
        sys.exit(0)

    # Run the two-stage test
    PASSED, H1_array = two_stage_test(data_base, nfiles, data_diff, files_base[0], dir_a)

    # Delete arrays that are no longer necessary
    del data_diff

    # If test failed, attempt to create a plot of the failure locations
    if not PASSED:
        plot_two_stage_failures(H1_array, t_lat, t_lon)
        logger.error('Quality Control Test FAILED')
        sys.exit(-1)

    # Create a northern hemisphere and southern hemisphere mask
    mask_tlat = t_lat < 0
    mask_nh = np.zeros_like(data_base)
    mask_sh = np.zeros_like(data_base)
    for (a, b), val in np.ndenumerate(mask_tlat):
        mask_nh[:, a, b] = val
        mask_sh[:, a, b] = not val

    # Run skill test on northern hemisphere
    data_nh_a = ma.masked_array(data_base, mask=mask_nh)
    data_nh_b = ma.masked_array(data_test, mask=mask_nh)
    if np.ma.all(data_nh_a.mask) and np.ma.all(data_nh_b.mask):
        logger.info("Northern Hemisphere data is bit-for-bit")
        PASSED_NH = True
    else:
        PASSED_NH = skill_test(dir_a, files_base[0], data_nh_a, data_nh_b, nfiles, 'Northern')

    # Run skill test on southern hemisphere
    data_sh_a = ma.masked_array(data_base, mask=mask_sh)
    data_sh_b = ma.masked_array(data_test, mask=mask_sh)
    if np.ma.all(data_sh_a.mask) and np.ma.all(data_sh_b.mask):
        logger.info("Southern Hemisphere data is bit-for-bit")
        PASSED_SH = True
    else:
        PASSED_SH = skill_test(dir_a, files_base[0], data_sh_a, data_sh_b, nfiles, 'Southern')

    PASSED_SKILL = PASSED_NH and PASSED_SH

    logger.info('')
    if not PASSED_SKILL:
        logger.error('Quality Control Test FAILED')
        sys.exit(1)  # exit with an error return code
    else:
        logger.info('Quality Control Test PASSED')
        sys.exit(0)  # exit with successfull return code

if __name__ == "__main__":
    main()
