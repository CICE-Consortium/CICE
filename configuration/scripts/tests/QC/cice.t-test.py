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
    mask = ~marr.mask.ravel()
    try:   # Python 2
        import itertools
        for i,m in itertools.izip(np.ndindex(marr.shape[-2:]),mask):
            if m: yield i
    except:  # Python 3
        for i,m in zip(np.ndindex(marr.shape[-2:]),mask):
            if m: yield i

def read_data(base_dir,test_dir):
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
    var = 'hi'
    cnt = 0
    for fname in sorted(files_a):
        nfid = nc.Dataset("{}/{}".format(path_a,fname),'r')
        fill_value_a = nfid.variables[var]._FillValue
        data_a[cnt,:,:] = nfid.variables[var][:]
        cnt += 1
        nfid.close()
    
    cnt = 0
    for fname in sorted(files_b):
        nfid = nc.Dataset("{}/{}".format(path_b,fname),'r')
        fill_value_b = nfid.variables[var]._FillValue
        data_b[cnt,:,:] = nfid.variables[var][:]
        cnt += 1
        nfid.close()
    
    # Calculate the difference and mask the points where the the difference at
    #   every timestep is 0.
    data_d = data_a - data_b
    mask_d = np.all(np.equal(data_d,0.),axis=0)
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
    for x in maenumerate(data_d):
        min_val = np.min(np.abs(df[x]-df_table))
        idx = np.where(np.abs(df[x]-df_table)==min_val)
        t_crit[x] = t_crit_table[idx]
    
    if np.any(abs(t_val) > t_crit):
        logger.info('Two-Stage Test Failed')
        passed = False
    
    elif np.all(abs(t_val)<=t_crit) and np.all(n_eff >= 30):
        logger.info('Two-Stage Test Passed')
        passed = True
    
    elif np.all(abs(t_val) <= t_crit) and np.any(n_eff < 30):
        # Calculate the T-statistic using actual sample size
        t_val = mean_d / np.sqrt(variance_d / num_files)
    
        # Find t_crit from the nearest value on the Lookup Table Test
        nfid = nc.Dataset("configuration/scripts/tests/QC/CICE_t_lookup_p0.8_n1825.nc",'r')
        r1_table = nfid.variables['r1'][:]
        t_crit_table = nfid.variables['tcrit'][:]
        nfid.close()
        
        for x in maenumerate(data_d):
            min_val = np.min(np.abs(r1[x]-r1_table))
            idx = np.where(np.abs(r1[x]-r1_table)==min_val)
            t_crit[x] = t_crit_table[idx]
    
        if np.any(abs(t_val) > t_crit):
            logger.info('Two-Stage Test Failed')
            passed = False
        elif np.all(abs(t_val) <= t_crit):
            logger.info('Two-Stage Test Passed')
            passed = True
        else:
            logger.error('TWO-STAGE TEST NOT CONCLUSIVE')
            passed = False
    
    else:
        logger.error('TEST NOT CONCLUSIVE')
        passed = False
    return passed

# Calculate Taylor Skill Score
def skill_test(path_a,fname,data_a,data_b,num_files,tlat,hemisphere):
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

def create_cdash_submit(tag):
    with open('submit.cmake','w') as f:
        f.write('cmake_minimum_required(VERSION 2.8)\n')
        f.write('set(CTEST_DASHBOARD_ROOT "$ENV{PWD}")\n')
        f.write('set(CTEST_SOURCE_DIRECTORY "$ENV{PWD}")\n')
        f.write('set(CTEST_BINARY_DIRECTORY "$ENV{PWD}")\n')
        f.write('\n')
        f.write('include(CTestConfig.cmake)\n')
        f.write('\n')
        f.write('ctest_start("Experimental")\n')
        f.write('ctest_submit(FILES "$ENV{PWD}/Testing/'+str(tag)+'/Test.xml")\n')
        f.write('\n')

def post_to_cdash(passed):
    logger.info('Posting QC Test results to CDash')
  
    # Create the CTest config file
    f = open('CTestConfig.cmake','w')
    f.write('set(CTEST_PROJECT_NAME       "myCICE")\n')
    f.write('set(CTEST_DROP_METHOD        "http")\n')
    f.write('set(CTEST_DROP_SITE          "my.cdash.org")\n')
    f.write('set(CTEST_DROP_LOCATION      "/submit.php?project=myCICE")\n')
    f.write('set(CTEST_DROP_SITE_CDASH    TRUE)\n')
    f.write('set(CTEST_NIGHTLY_START_TIME "01:00:00 CET")\n')
    f.close()

    # Create the CTest test file
    f = open('CTestTestfile.cmake','w')
    f.write('add_test(Two-Stage-Test grep \"Two-Stage Test Passed\" qc_log.txt)\n')
    f.write('add_test(Skill-Northern-Hemi grep \"Passed for Northern\" qc_log.txt)\n')
    f.write('add_test(Skill-Southern-Hemi grep \"Passed for Southern\" qc_log.txt)\n')
    f.close()

    # Create the steer.cmake file
    f = open('steer.cmake','w')
    f.write('## -- Set hostname\n')
    f.write('find_program(HOSTNAME_CMD NAMES hostname)\n')
    f.write('exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)\n')
    f.write('set(CTEST_SITE    "${HOSTNAME}")\n')
    f.write('\n')
    f.write('## -- Set site / build name\n')
    f.write('find_program(UNAME NAMES uname)\n')
    f.write('macro(getuname name flag)\n')
    f.write('  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")\n')
    f.write('endmacro(getuname)\n')
    f.write('getuname(osname -s)\n')
    f.write('getuname(osrel -r)\n')
    f.write('getuname(cpu -m)\n')
    f.write('\n')
    f.write('find_program(GIT_CMD NAMES git)\n')
    f.write('exec_program(${GIT_CMD} ARGS rev-parse --short HEAD OUTPUT_VARIABLE GIT_COMMIT_HASH\n)')
    f.write('\n')
    f.write('set(CTEST_BUILD_NAME "${osname}-${cpu}-QC_TEST-${GIT_COMMIT_HASH}")\n')
    f.write('\n')
    f.write('set(CTEST_DASHBOARD_ROOT "$ENV{PWD}")\n')
    f.write('set(CTEST_SOURCE_DIRECTORY "$ENV{PWD}")\n')
    f.write('set(CTEST_BINARY_DIRECTORY "$ENV{PWD}")\n')
    f.write('\n')
    f.write('ctest_start(${MODEL} TRACK ${MODEL})\n')
    f.write('ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)\n')
    f.write('ctest_submit(RETURN_VALUE res)\n')
    f.close()

    import subprocess
    subprocess.call(['ctest','-S','steer.cmake'])

    # Check CTest log to see if submit succeeded
    with open('Testing/TAG','r') as f:
        tag = f.readline()
        tag = tag.rstrip('\n')

    if 'Submission successful' in open('Testing/Temporary/LastSubmit_{}.log'.format(tag)).read():
        success = True
    else:
        success = False

    # If submit failed, create backup
    if not success:
        logger.error('')
        logger.error('')
        logger.error('')
        logger.error('Failed to submit results to CDash.  ')
        logger.error('')
        logger.error('To try the submission again, copy cice_qc_ctest.tgz to another server, ')
        logger.error('and run:')
        logger.error('    ctest -S submit.cmake')
        create_cdash_submit(tag)
        import shutil
        shutil.copyfile('Testing/TAG','Testing/TAG.submit')
        import tarfile
        tar = tarfile.open('cice_qc_ctest.tgz','w:gz')
        for name in ['Testing/TAG.submit','CTestConfig.cmake','Testing/{}/Test.xml'.format(tag),\
                     'submit.cmake']:
            tar.add(name)
        tar.close()

    # Delete the ctest / cdash files and folders
    import os
    os.remove('CTestConfig.cmake')
    os.remove('CTestTestfile.cmake')
    os.remove('steer.cmake')
    import shutil
    shutil.rmtree('Testing/')

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
    
    data_a, data_b, data_d, num_files, path_a, fname = read_data(args.base_dir,args.test_dir)

    if np.ma.all(data_d.mask):
        logger.info("Data is bit-for-bit.  No need to run QC test")
        sys.exit(0)

    # Run the two-stage test
    passed = two_stage_test(data_a,data_b,num_files,data_d)
    
    nfid = nc.Dataset("{}/{}".format(path_a,fname),'r')
    tlat = nfid.variables['TLAT'][:]
    nfid.close()
    
    mask_tlat = tlat < 0
    mask_nh = np.zeros_like(data_a)
    mask_sh = np.zeros_like(data_a)
    for (i,j),value in np.ndenumerate(mask_tlat):
        mask_nh[:,i,j] = value
        mask_sh[:,i,j] = not value
    
    # Run skill test on northern hemisphere
    data_nh_a = ma.masked_array(data_a,mask=mask_nh)
    data_nh_b = ma.masked_array(data_b,mask=mask_nh)
    passed_nh = skill_test(path_a,fname,data_nh_a,data_nh_b,num_files,tlat,'Northern')
    
    # Run skill test on southern hemisphere
    data_sh_a = ma.masked_array(data_a,mask=mask_sh)
    data_sh_b = ma.masked_array(data_b,mask=mask_sh)
    passed_sh = skill_test(path_a,fname,data_sh_a,data_sh_b,num_files,tlat,'Southern')
    
    passed_skill = passed_nh and passed_sh
    
    logger.info('')
    if not passed and not passed_skill:
        logger.error('Quality Control Test FAILED')
        post_to_cdash(False)
        sys.exit(1)  # exit with an error return code
    else:
        logger.info('Quality Control Test PASSED')
        post_to_cdash(True)
        sys.exit(0)  # exit with successfull return code
    
