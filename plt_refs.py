import os, glob
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from wrf import (getvar, to_np, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, CoordPair, vertcross, ll_to_xy)
from netCDF4 import Dataset

import cartopy.crs as crs
import cartopy.feature as cfeat 
from cartopy.mpl.ticker import LongitudeFormatter,LatitudeFormatter
from cartopy.io.shapereader import Reader
#from cartopy.feature import NaturalEarthFeature

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors, BoundaryNorm
import matplotlib.ticker as mticker
import matplotlib.path as mpath

SHP = r'/share/home/zhaokun/.local/share/cartopy/shapefiles/natural_earth/china_shp'

dbz_levels = np.arange(5., 75., 5.)
# Create the color table found on NWS pages.
dbz_rgb = np.array([[4,233,231],
                    [1,159,244], [3,0,244],
                    [2,253,2], [1,197,1],
                    [0,142,0], [253,248,2],
                    [229,188,0], [253,149,0],
                    [253,0,0], [212,0,0],
                    [188,0,0],[248,0,253],
                    [152,84,198]], np.float32) / 255.0
dbz_map, dbz_norm = from_levels_and_colors(dbz_levels, dbz_rgb,
                                           extend="max")

def plt_rf(ax, x, y, variable, cur_time = None, ini_time = None, \
           var_name = 'composite reflectivity (dBZ)', extend = [106.,122.,15.,26.], \
           var_levels = None, var_map = None, isLeft = True, isBottom = True, sub_label = None):

    '''
    plot the reflectivity from wrf file. 
    
    extend: list
        the extend of the lat/lon, in the format with [min lon, max lon, min lat, max lon]
    '''

    ax.add_geometries(Reader(os.path.join(SHP, 'cnhimap.shp')).geometries(),
                      crs.PlateCarree(),facecolor='none',edgecolor='k', linewidth=0.5) 

    dbz_contours = ax.contourf(x, y, variable, 10,
                transform=crs.PlateCarree(),
                levels=var_levels,
                cmap=var_map)
    
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())

    # Add the gridlines
    gl = ax.gridlines(color="black", linestyle="dotted", draw_labels=True)
    gl.top_labels  = False
    gl.right_labels = False
    gl.bottom_labels = True if isBottom else False
    gl.left_labels = True if isLeft else False

    gl.xlocator = mticker.FixedLocator(np.arange(extend[0]+1,extend[1]+3,2))
    gl.ylocator = mticker.FixedLocator(np.arange(extend[2]+1,extend[3]+0.1,2))

    gl.y_inline = False # lat/lon doesn't show on the grid lines.
    gl.x_inline = False

    gl.xlabel_style  = {'rotation':'vertical', 'size':6}
    gl.ylabel_style  = {'rotation':'horizontal', 'size':6}
    
    if sub_label:
        ax.text(0.05, 0.9, sub_label,  transform=ax.transAxes, fontsize=7, color='tab:red')

    ax.set_extent(extend)

def plt_prcp(ax, x, y, variable, cur_time = None, ini_time = None, \
           var_name = 'composite reflectivity (dBZ)', extent = [106.,122.,15.,26.], \
           var_levels = None, norm = None, var_map = None, isLeft = True, isBottom = True):
    
    ax.add_geometries(Reader(os.path.join(SHP, 'cnhimap.shp')).geometries(),
                      crs.PlateCarree(),facecolor='none',edgecolor='k', linewidth=0.5) 

    sc = ax.scatter(x, y, c=variable, cmap=var_map, marker='o', transform=crs.PlateCarree(), s=2, norm=norm)
    cb = fig.colorbar(sc, ax=ax, label='precipitation (mm/h)')
    cb.ax.tick_params(labelsize=8)
    ax.set_extent(extent,crs=crs.PlateCarree())

    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())

    # Add the gridlines
    gl = ax.gridlines(color="black", linestyle="dotted", draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True if isBottom else False
    gl.left_labels = True if isLeft else False


    gl.xlocator = mticker.FixedLocator(np.arange(extent[0]+1,extent[1]+3,2))
    gl.ylocator = mticker.FixedLocator(np.arange(extent[2]+1,extent[3]+0.1,2))

    gl.y_inline = False # lat/lon doesn't show on the grid lines.
    gl.x_inline = False
    gl.xlabel_style  = {'rotation':'horizontal'}
    plt.text(0.72, 0.03, 'max precipitation: ' + str(np.max(variable)), horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes, fontsize=12)

    return ax
    


if  __name__ == '__main__':
    exp_names = ['OBS','3dvar_no_urban','3dvar_urban'] #'radar_1h_cv7', 'hybrid_radar_1h','hybrid_noVinEnkf']
    fcst_hrs = [1, 2, 3, 4, 5, 6]

    cur_date = datetime(2023,9,6,12,0)
    end_date = datetime(2023,9,7,12,0)
    extend = [107.,120.1,17.,26.]

    OBS_dir = '/share/home/zhaokun/haiqin/data/radar1hr/202309'
    pathin = '/share/home/zhaokun/scratch/haiqin/exp/Haikui/'
    pathout = '/share/home/zhaokun/haiqin/scripts/plot/saola_mdbz'

    file_pattern = 'wrfout_d02_'

    n_exp = len(exp_names)
    n_times = len(fcst_hrs)
    while cur_date <= end_date:

        fig, axes = plt.subplots(n_times,n_exp,figsize=(n_exp*1.4,n_times),dpi=600,
        sharey=True,sharex=True,
        constrained_layout=True,
        subplot_kw={'projection':crs.PlateCarree()})

        for i, fcst_hr in enumerate(fcst_hrs):
            for j, i_exp in enumerate(exp_names):
                cur_time = cur_date + timedelta(hours=int(fcst_hr))

                if i_exp == 'OBS':
                    ncfile = Dataset(OBS_dir+'/CINRAD.'+cur_time.strftime('%Y%m%d%H%M')+'.nc')
                    obs_lat, obs_lon = ncfile.variables['lat'], ncfile.variables['lon']
                    lats, lons = obs_lat[:], obs_lon[:]
                    hrv = ncfile.variables['HRV']
                    mdbz = hrv[1,:,:]
                else:
                    wrfout_path = os.path.join(pathin,i_exp, 'fc', cur_date.strftime('%Y%m%d%H'))
                    if not os.path.exists(wrfout_path): 
                        raise Exception("wrfout directory doesn't exist")
                    ncfile = Dataset(wrfout_path+'/'+file_pattern+cur_time.strftime('%Y-%m-%d_%H:%M:%S'))
                    mdbz = np.max(ncfile.variables['REFL_10CM'][0,:,:,:],axis=0)
                    if i == 0:
                        tc = getvar(ncfile, "tc")
                        wrf_lats, wrf_lons = latlon_coords(tc)
                    lats, lons = to_np(wrf_lats), to_np(wrf_lons)

                isLeft = True if j == 0 else False
                isBottom = True if i == (n_times - 1) else False
                plt_rf(ax = axes[i,j],x =lons, y = lats, variable = mdbz,\
                    cur_time = cur_time.strftime('%Y%m%d%H%M'), ini_time = cur_date.strftime('%Y%m%d%H%M'), var_name = 'composite reflectivity (dBZ)', \
                    extend = extend, var_levels = dbz_levels, var_map = dbz_map, isLeft = isLeft, isBottom = isBottom)

        fc = fig.colorbar(
        mpl.cm.ScalarMappable(norm=dbz_norm, cmap=dbz_map),
        ax=axes,
        shrink=0.5,
        extendrect=True,
        orientation='vertical')
        fc.ax.tick_params(labelsize=8)
        
        fig.suptitle("mdBZ initiated from "+ cur_date.strftime('%Y%m%d%H%M'), color='black', fontsize=10)
        plt.savefig(pathout+'/mdbzs_' + cur_date.strftime('%Y%m%d%H%M')+ '.png', bbox_inches='tight', dpi=200)
        plt.show()
        plt.close()

        cur_date += timedelta(hours=12)
