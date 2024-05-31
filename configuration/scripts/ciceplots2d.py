#!/usr/bin/env python3

#Importing the necessary libraries
import sys
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

if len(sys.argv) != 6:
    print("ciceplots.py requires 5 arguments")
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
title = field + " " + notes
cfnam = casen + " " + fname
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

lons = data.variables['TLON'][:,:]
lats = data.variables['TLAT'][:,:]
var1 = data.variables[field][:,:,:]
var1 = var1[0,:,:]
var1[ var1==0.00 ] = np.nan
#mask = data.variables['tmask'][:,:]
#mask[ mask>0.5 ] = np.nan

#print("lons.shape = ",lons.shape)
#print("var1.shape = ",var1.shape)

# Lon/Lat Projection

#print("Plot global")
#m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,
#            llcrnrlon=0,urcrnrlon=360,resolution='c')
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,
            llcrnrlon=0,urcrnrlon=360,resolution='l')
fig, ax = plt.subplots()
#plt.figure(figsize=(6,4))
m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='black',lake_color='white')
#draw parallels and meridians.
m.drawparallels(np.arange(-60.,61.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,361.,45.),labels=[1,0,0,1])
#draw map boundary
m.drawmapboundary(fill_color='white')
#setting colorbar
cmap = plt.get_cmap('jet')
barticks = None
norm = "linear"
if field in ['hi']:
    bounds = np.arange(0,2.05,0.1)
    bounds = np.append(bounds,[2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0])
    norm = mpl.colors.BoundaryNorm(bounds,cmap.N,extend='max')
    barticks=[0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0]
if field in ['hs']:
    bounds = np.arange(0,1.02,0.05)
    bounds = np.append(bounds,[1.5,2.0,2.5,3.0,3.5,4.0])
    norm = mpl.colors.BoundaryNorm(bounds,cmap.N,extend='max')
    barticks=[0,0.25,0.5,0.75,1.0,2.0,3.0,4.0]
#matplotlib scatter-plot 
m.scatter(lons,lats,c=var1,cmap=cmap,marker='o',s=0.2,norm=norm)
m.colorbar(label=field, ticks=barticks)
plt.rcParams["figure.dpi"] = 300
plt.title(title)
plt.text(x=0.0,y=-0.1,s=cfnam,transform=ax.transAxes,horizontalalignment='left',verticalalignment='top',fontsize='x-small')
oname = field + "_gl_" + fstr + ".png"
print('Saving file to ',oname)
plt.savefig(oname)
#plt.show()
plt.close()

# North Polar Stereographic Projection

#print("Plot NH")
#m = Basemap(projection='npstere',boundinglat=45,lon_0=-45,resolution='c')
m = Basemap(projection='npstere',boundinglat=45,lon_0=-45,resolution='l')
fig, ax = plt.subplots()
#plt.figure(figsize=(6,4))
m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='black',lake_color='white')
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,61.,30.),labels=[0,0,0,0])
m.drawmeridians(np.arange(0.,361.,45.),labels=[0,0,0,0])
m.drawmapboundary(fill_color='white')
#setting colorbar (set above)
m.scatter(lons,lats,c=var1,cmap=cmap,marker='o',s=0.2,latlon=True,norm=norm)
#m.colorbar(label=field)
m.colorbar(label=field, ticks=barticks)
plt.rcParams["figure.dpi"] = 300
plt.title (title)
plt.text(x=0.0,y=-0.02,s=cfnam,transform=ax.transAxes,horizontalalignment='left',verticalalignment='top',fontsize='x-small')
oname = field + "_nh_" + fstr + ".png"
print('Saving file to ',oname)
plt.savefig(oname)
#plt.show()
plt.close()

# South Polar Stereographic Projection

#print("Plot SH")
#m = Basemap(projection='npstere',boundinglat=45,lon_0=-45,resolution='c')
m = Basemap(projection='spstere',boundinglat=-45,lon_0=180,resolution='l')
fig, ax = plt.subplots()
#plt.figure(figsize=(6,4))
m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='black',lake_color='white')
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,61.,30.),labels=[0,0,0,0])
m.drawmeridians(np.arange(0.,361.,45.),labels=[0,0,0,0])
m.drawmapboundary(fill_color='white')
#setting colorbar (set above)
m.scatter(lons,lats,c=var1,cmap=cmap,marker='o',s=0.2,latlon=True,norm=norm)
#m.colorbar(label=field)
m.colorbar(label=field, ticks=barticks)
plt.rcParams["figure.dpi"] = 300
plt.title (title)
plt.text(x=0.0,y=-0.02,s=cfnam,transform=ax.transAxes,horizontalalignment='left',verticalalignment='top',fontsize='x-small')
oname = field + "_sh_" + fstr + ".png"
print('Saving file to ',oname)
plt.savefig(oname)
#plt.show()
plt.close()

#print("Done")
quit()

