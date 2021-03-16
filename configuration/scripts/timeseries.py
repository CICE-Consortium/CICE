#!/usr/bin/env python

'''
This script generates timeseries plots of CICE diagnostic output.
It is generated to replicate the previous timeseries.csh script.

Written by: Matthew Turner
Date: August, 2019
'''

import os
import sys
import logging
import numpy as np

def find_logfile(log_dir):
   '''
   This function searches for the most recently created log file in the provided directory.
   '''

   logger.debug('Getting a list of files in {}'.format(log_dir))
   try:
       path = '{}/logs'.format(log_dir.rstrip('/'))
       files = [os.path.join(path,f) for f in os.listdir('{}/logs'.format(log_dir)) \
                if f.startswith('cice.runlog')]
   except:
       path = log_dir
       files = [os.path.join(path,f) for f in os.listdir(log_dir) if f.startswith('cice.runlog')]

   # Check if any files were found.  If not, exit
   if len(files) == 0:
      logger.error('No cice.runlog* files found.  Please make sure you are passing the \
                    correct directory.')
      sys.exit(1)

   # Get the most recently created file
   outfile = max(files, key = os.path.getctime)

   logger.debug('List of files = {}'.format([f for f in files]))
   logger.debug('Most recent file is {}'.format(outfile))

   return outfile

def get_data(logfile,field):
    '''
    This function extracts data from a CICE log file for the specific field.
    '''
    import datetime
    import re

    logger.debug('Extracting data for {}'.format(field))

    # Build the regular expression to extract the data
    field_regex = field.replace('(','\(').replace('^','\^').replace(')','\)')
    number_regex = '[-+]?\d+\.?\d+([eE][-+]?\d+)?'
    my_regex = '^{}\s+=\s+({})\s+({})'.format(field_regex,number_regex,number_regex)

    dtg = []
    arctic = []
    antarctic = []
    with open(logfile) as f:
        for line in f.readlines():
            m1 = re.search('istep1:\s+(\d+)\s+idate:\s+(\d+)\s+sec:\s+(\d+)', line)
            if m1:
                # Extract the current date-time group from the file
                date = m1.group(2)
                seconds = int(m1.group(3))
                hours = seconds // 3600
                minutes = (seconds - hours*3600) // 60
                leftover = seconds - hours*3600 - minutes*60
                curr_date = '{}-{:02d}:{:02d}:{:02d}'.format(date,hours,minutes,leftover)
                dtg.append(datetime.datetime.strptime(curr_date, '%Y%m%d-%H:%M:%S'))
                logger.debug('Currently on timestep {}'.format(dtg[-1]))

            m = re.search(my_regex, line)
            if m:
                # Extract the data from the file
                if 'E' in m.group(1) or 'e' in m.group(1):
                    expon = True
                else:
                    expon = False
                arctic.append(float(m.group(1)))
                antarctic.append(float(m.group(3)))
                logger.debug('    Arctic = {}, Antarctic = {}'.format(arctic[-1], antarctic[-1]))

    return dtg, arctic, antarctic, expon

def latexit(string):
    s = string[::-1].replace('(','($',1)
    return (s.replace(')','$)',1))[::-1]

def plot_timeseries(log, field, dtg, arctic, antarctic, expon, dtg_base=None, arctic_base=None, \
                    antarctic_base=None, base_dir=None, grid=False):
    '''
    Plot the timeseries data from the CICE log file
    '''

    import re
    casename = re.sub(r"/logs", "", os.path.abspath(log).rstrip('/')).split('/')[-1]
    if base_dir:
        base_casename = re.sub(r"/logs", "", os.path.abspath(base_dir).rstrip('/')).split('/')[-1]

    # Load the plotting libraries, but set the logging level for matplotlib
    # to WARNING so that matplotlib debugging info is not printed when running
    # with '-v'
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import matplotlib.ticker as ticker

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_axes([0.05,0.08,0.9,0.9])
    
    # Add the arctic data to the plot
    ax.plot(dtg,arctic,label='Arctic')
    # Add the baseline arctic data to the plot, if available
    if arctic_base:
        ax.plot(dtg_base,arctic_base,label='Baseline Arctic')

    # Add the antarctic data to the plot
    ax.plot(dtg,antarctic,label='Antarctic')
    # Add the baseline antarctic data to the plot, if available
    if antarctic_base:
        ax.plot(dtg_base,antarctic_base,label='Baseline Antarctic')

    ax.set_xlabel('')
    ax.set_title('{} Diagnostic Output'.format(latexit(field)))
    ax.set_ylabel(latexit(field))

    # Format the x-axis labels
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

    # Add a text box that prints the test case name and the baseline case name (if given)
    try:
        text_field = "Test/Case: {}\nBaseline: {}".format(casename,base_casename)
        from matplotlib.offsetbox import AnchoredText
        anchored_text = AnchoredText(text_field,loc=2)
        ax.add_artist(anchored_text)
    except:
        text_field = "Test/Case: {}".format(casename)
        from matplotlib.offsetbox import AnchoredText
        anchored_text = AnchoredText(text_field,loc=2)
        ax.add_artist(anchored_text)

    ax.legend(loc='upper right')

    # Add grid lines if the `--grid` argument was passed at the command line.
    if grid:
        ax.grid(ls='--')

    # Reduce the number of ticks on the y axis
    nbins = 10
    try:  
        minval = min( \
                     min(min(arctic), min(antarctic)), \
                     min(min(arctic_base), min(antarctic_base)))
        maxval = max( \
                     max(max(arctic), max(antarctic)), \
                     max(max(arctic_base), max(antarctic_base)))
    except:
        minval = min(min(arctic), min(antarctic))
        maxval = max(max(arctic), max(antarctic))
    step = (maxval-minval)/nbins
    ax.yaxis.set_ticks(np.arange(minval, maxval+step, step))

    # Format the y-axis tick labels, based on whether or not the values in the log file
    # are in scientific notation or float notation.
    if expon:
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3e'))
    else:
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.5f'))

    # Rotate and right align the x labels
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    # Create an output file and save the figure
    field_tmp = field.split('(')[0].rstrip()
    try:
        outfile = '{}_{}_base-{}.png'.format(field_tmp.replace(' ','_'), casename,base_casename)
    except:
        outfile = '{}_{}.png'.format(field_tmp.replace(' ','_'), casename)
    logger.info('Saving file to {}'.format(outfile))
    plt.savefig(outfile,dpi=300,bbox_inches='tight')

def main():
    import argparse
    parser = argparse.ArgumentParser(description="To generate timeseries plots, this script \
                                     can be passed a directory containing a logs/ subdirectory, \
                                     or it can be run in the directory with the log files, \
                                     without being passed a directory.  It will pull the \
                                     diagnostic data from the most recently modified log file.\
                                     \
                                     If no flags are passed selecting the variables to plot, \
                                     then plots will be created for all available variables.")
    parser.add_argument('log_dir', nargs='?', default=os.getcwd(), \
                        help="Path to diagnostic output log file.  A specific log file can \
                              be passed, or a case directory.  If a directory is passed, \
                              the most recent log file will be used.  If no directory or \
                              file is passed, the script will look for a log file in the \
                              current directory.")
    parser.add_argument('--bdir',dest='base_dir', help='Path to the the log file for a baseline \
                              dataset, if desired.  A specific log file or case directory can \
                              be passed.  If a directory is passed, the most recent log file \
                              will be used.')
    parser.add_argument('-v', '--verbose', dest='verbose', help='Print debug output?', \
                        action='store_true')
    parser.add_argument('--area', dest='area', help='Create a plot for total ice area?', \
                        action='store_true')
    parser.add_argument('--extent', dest='extent', help='Create a plot for total ice extent?', \
                        action='store_true')
    parser.add_argument('--volume', dest='ice_volume', help='Create a plot for total ice volume?', \
                        action='store_true')
    parser.add_argument('--snw_vol', dest='snow_volume', help='Create a plot for total snow \
                        volume?', action='store_true')
    parser.add_argument('--speed', dest='speed', help='Create a plot for rms ice speed?', \
                        action='store_true')
    parser.add_argument('--grid',dest='grid', help='Add grid lines to the figures?', \
                        action='store_true')

    # Set the defaults for the command line options
    parser.set_defaults(verbose=False)
    parser.set_defaults(area=False)
    parser.set_defaults(extent=False)
    parser.set_defaults(ice_volume=False)
    parser.set_defaults(snow_volume=False)
    parser.set_defaults(speed=False)
    parser.set_defaults(grid=False)

    args = parser.parse_args()

    # If no fields are passed, plot all fields
    if not ( args.area or args.extent or args.ice_volume or args.snow_volume or args.speed ):
        args.area = True
        args.extent = True
        args.ice_volume = True
        args.snow_volume = True
        args.speed = True

    # Build the fieldlist based on which fields are passed
    fieldlist = []
    if args.area:
        fieldlist.append('total ice area  (km^2)')
    if args.extent:
        fieldlist.append('total ice extent(km^2)')
    if args.ice_volume:
        fieldlist.append('total ice volume (m^3)')
    if args.snow_volume:
        fieldlist.append('total snw volume (m^3)')
    if args.speed:
        fieldlist.append('rms ice speed    (m/s)')

    # Setup the logger
    global logger
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Find the test and baseline log files, based on the input directories.
    if os.path.isdir(args.log_dir):
        logger.debug('{} is a directory'.format(args.log_dir))
        log = find_logfile(args.log_dir)
        log_dir = args.log_dir
    else:
        logger.debug('{} is a file'.format(args.log_dir))
        log = args.log_dir
        log_dir = args.log_dir.rsplit('/',1)[0]
    logger.info('Log file = {}'.format(log))
    if args.base_dir:
        if os.path.isdir(args.base_dir):
            base_log = find_logfile(args.base_dir)
            base_dir = args.base_dir
        else:
            base_log = args.base_dir
            base_dir = args.base_dir.rsplit('/',1)[0]
        logger.info('Base Log file = {}'.format(base_log))

    # Loop through each field and create the plot
    for field in fieldlist:
        logger.debug('Current field = {}'.format(field))

        # Get the data from the log files
        dtg, arctic, antarctic, expon = get_data(log, field)
        if args.base_dir:
            dtg_base, arctic_base, antarctic_base, expon_base = get_data(base_log,field)

        # Plot the data
        if args.base_dir:
            plot_timeseries(log_dir, field, dtg, arctic, antarctic, expon, dtg_base, \
                            arctic_base, antarctic_base, base_dir, grid=args.grid)
        else:
            plot_timeseries(log_dir, field, dtg, arctic, antarctic, expon, grid=args.grid)

if __name__ == "__main__":
    main()
