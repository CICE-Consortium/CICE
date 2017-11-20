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
import itertools
import logging

def maenumerate(marr):
    mask = ~marr.mask.ravel()
    for i,m in itertools.izip(np.ndindex(marr.shape[-2:]),mask):
        if m: yield i

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
logger = logging.getLogger(__name__)

# The path to output files for simulation 'a' (the '-bc' simulation)
if args.base_dir.endswith(('history','history/')):
    path_a = args.base_dir
else:
    path_a = args.base_dir + '/history/'

# The path to output files for simulation 'b' (the test simulation)
if args.test_dir.endswith(('history','history/')):
    path_b = args.test_dir
else:
    path_b = args.test_dir + '/history/'

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

# Calculate the mean of the difference
mean_d = np.mean(data_d,axis=0)
variance_d = np.sum(np.power(data_d - mean_d,2)) / (num_files - 1)

# Calculate the mean from 1:end-1 and 2:end
mean_nm1_d = np.mean(data_d[:-1,:,:],axis=0)
mean_2n_d = np.mean(data_d[1:,:,:],axis=0)

# Calculate equation (5) for both simulations
r1_num = np.zeros((nj,ni))
r1_den1 = np.zeros((nj,ni))
r1_den2 = np.zeros((nj,ni))
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
nfid = nc.Dataset("CICE_t_critical_p0.8.nc",'r')
df_table = nfid.variables['df'][:]
t_crit_table = nfid.variables['tcrit'][:]
nfid.close()
t_crit = np.zeros_like(data_d)
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
    nfid = nc.Dataset("CICE_t_lookup_p0.8_n1825.nc",'r')
    r1_table = nfid.variables['r1'][:]
    t_crit_table = nfid.variables['tcrit'][:]
    nfid.close()
    
    for x in maenumerate(data_d):
        min_val = np.min(np.abs(r1[x]-r1_table))
        idx = np.where(np.abs(r1[x]-r1_table)==min_val)
        t_crit[x] = t_crit_table[idx]

    logger.debug('Max t_crit 1 = {}'.format(np.max(t_crit)))
    logger.debug('Min t_crit 1 = {}'.format(np.min(t_crit)))
    
    if np.any(abs(t_val) > t_crit):
        logger.info('Two-Stage Test Failed')
        passed = False
    elif np.all(abs(t_val) <= t_crit):
        logger.info('Two-Stage Test Passed')
        passed = True
    else:
        logger.error('TEST NOT CONCLUSIVE')

else:
    logger.error('TEST NOT CONCLUSIVE')

# Calculate Taylor Skill Score
#dcorr = np.zeros_like(data_d)
dcorr = np.zeros((nj,ni))
dcorr[:,:] = 1
for i,j in maenumerate(data_d):
    dcorr[i,j] = np.corrcoef(data_a[:,i,j],data_b[:,i,j])[0,1]
s = np.power( \
       (1+dcorr)*np.std(data_a,axis=0)*np.std(data_b,axis=0) / \
       (np.power(np.std(data_a,axis=0),2) + np.power(np.std(data_b,axis=0),2)) \
       , 2)
s_crit = 0.80
if np.any(s<0) or np.any(s>1):
    logger.error('Skill score out of range')
elif np.all(s > s_crit):
    logger.info('Quadratic Skill Test Passed')
else:
    logger.info('Quadratic Skill Test Failed')
    passed = False

if passed == False:
    sys.exit(1)  # exit with an error return code
else:
    sys.exit(0)  # exit with successfull return code

