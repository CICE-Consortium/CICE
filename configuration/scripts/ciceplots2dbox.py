#!/usr/bin/env python3

#Importing the necessary libraries
import sys
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

if len(sys.argv) != 6:
    print("ciceplots2d.py requires 5 arguments")
    print("  1. field name in file, ie. \"aice\"")
    print("  2. cice history file full path, ie. \"/glade/scratch/user/case/history/iceh.2012-03.nc\"")
    print("  3. case name, used to annotate plot, ie. \"CICE6.5.1\"")
    print("  4. notes, used to annotate plot, ie. 2012 \"March Mean\"")
    print("  5. file string, use to create unique png filenames, ie. \"Mar12\"")
    quit()

field = sys.argv[1]
pathf = sys.argv[2]
casen = sys.argv[3]
notes = sys.argv[4]
fstr  = sys.argv[5]
fname = os.path.basename(pathf)
title = casen + " " + field + " " + notes
cfnam = casen + "_" + fname
#print("field = ",field)
#print("pathf = ",pathf)
#print("casen = ",casen)
#print("notes = ",notes)
#print("fname = ",fname)
#print("title = ",title)
#print("cfnam = ",cfnam)

#Reading the netCDF file
data = Dataset(pathf,'r')
#print (data)

#lons = data.variables['TLON'][:,:]
#lats = data.variables['TLAT'][:,:]
var1 = data.variables[field][:,:,:]
var1 = var1[0,:,:]
var1[ var1==0.00 ] = np.nan
#mask = data.variables['tmask'][:,:]
#mask[ mask>0.5 ] = np.nan

#print("lons.shape = ",lons.shape)
#print("var1.shape = ",var1.shape)

# initialize the plot
fig, ax = plt.subplots(figsize=(10, 6), dpi=300) # 300 DPI ensures crisp text and lines

# generate the filled color contour map
#contour_filled = ax.contourf(lons, lats, var1, levels=20, cmap='coolwarm')
contour_filled = ax.contourf(var1, levels=20, cmap='coolwarm')

# overlay black contour lines
#contour_lines = ax.contour(lons, lats, var1, levels=20, colors='black', linewidths=0.3)
contour_lines = ax.contour(var1, levels=20, colors='black', linewidths=0.3)
ax.clabel(contour_lines, inline=True, fontsize=8, fmt='%.1f')

# include a colorbar
cmap = plt.get_cmap('jet')
barticks = None
fieldunits = ' '
if field in ['hi']:
    barticks=[0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0]
    fieldunits='m'
cbar = fig.colorbar(contour_filled, ax=ax, ticks=barticks)
cbar.set_label(f"{field} ({fieldunits})")

# define labels and layout
ax.set_title(title)
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.tight_layout()

# save as a .png image file
# Use bbox_inches='tight' so labels do not get clipped out of the frame
oname = casen + "_" + field + "_" + fstr + ".png"
print('Saving file to ',oname)
plt.savefig(oname, format='png', dpi=300, bbox_inches='tight')

# finish
plt.close(fig)
print("Success: 'contour plot has been saved to your working directory.")

quit()

