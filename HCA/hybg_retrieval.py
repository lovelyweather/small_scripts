import xarray as xr
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib.colors as colors
import matplotlib
import xarray as xr
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


a1=3.63*1e9
a2=4.26*1e11 #wet snow
a3=9.6*1e8
a4=6.12*1e10
eps_phy = 0.622
rd  = 287.04



'''plot 0.0'''
hid_colors = ['MediumBlue', 'LightPink',
              'Cyan',  'Lime']
cmaphid = colors.ListedColormap(hid_colors)
def plot_lonlat_map(ax, lon, lat, data, transform, extend=None, cmap="CN_ref", bounds=np.arange(0.5,4.6,1),
                    cbar=True, orientation="vertical",cbar_ticks=None, cbar_ticklabels=None, **kwargs):
    """
    :param ax:cartopy.mpl.geoaxes.GeoAxesSubplot, it should get from cartopy, eg:plt.axes(projection=ccrs.PlateCarree())
    :param lon: lon mesh grid for data units: degree
    :param lat: lat mesh grid for data units: degree
    :param data: radar data ,dims like lat, lon
    :param transform: The transform argument to plotting functions tells Cartopy what coordinate system your data are defined in.
    :param extend: (min_lon, max_lon, min_lat, max_lat), Latitude and longitude range, units:degrees
    :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
    :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
    :param cmap_bins:bins of cmaps, int
    :param cbar: bool, if True, plot with colorbar,
    :param orientation: vertical or horizontal, it is vaild when cbar is True
    :param kwargs: kwargs: other arguments for pcolormesh!
    :return:  pcolor result
    """
    assert isinstance(ax, cartopy.mpl.geoaxes.GeoAxesSubplot), "axes is not cartopy axes!"
    if extend is None:
        min_lon = np.min(lon)
        max_lon = np.max(lon)
        min_lat = np.min(lat)
        max_lat = np.max(lat)
    else:
        min_lon, max_lon, min_lat, max_lat = extend

  #  ax.set_aspect("equal")
    cmaps = plt.get_cmap(cmap)
    norm = BoundaryNorm(bounds, ncolors=cmaps.N, clip=True)
    pm = ax.pcolormesh(lon, lat, data, transform=transform, cmap=cmap, norm=norm, zorder=1, **kwargs)
#    ax.add_feature(cfeature.OCEAN.with_scale('50m'), zorder=1)
#    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', \
#                                                edgecolor='none', facecolor="white"), zorder=2)
    ax.add_feature(cfeature.LAKES.with_scale('50m'), zorder=1)
    ax.add_feature(cfeature.RIVERS.with_scale('50m'), zorder=1)
    ax.add_feature(cfeature.STATES,zorder=1)

    #ax.add_feature(cfeature.ShapelyFeature(CN_shp_info.geometries(), transform, \
    #                                       edgecolor='k', facecolor='none'), linewidth=0.5, \
    #               linestyle='-', zorder=5, alpha=0.8)
    parallels = np.arange(int(min_lat), np.ceil(max_lat) + 1, 1)
    meridians = np.arange(int(min_lon), np.ceil(max_lon) + 1, 1)
    ax.set_xticks(meridians, crs=transform)
    ax.set_yticks(parallels, crs=transform)
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_extent(extend, crs=transform)
    ax.spines['left'].set_visible(True)
    ax.spines['top'].set_visible(True)

    if cbar:
        cb = plt.colorbar(mappable=pm, ax=ax, orientation=orientation)
        if cbar_ticks is None:
            ticks = bounds
        else:
            ticks = cbar_ticks  #1,2,3,4,5...
        cb.set_ticks(ticks)

    if cbar_ticklabels is not None:
        if orientation == "vertical":
            cb.ax.set_yticklabels(cbar_ticklabels)
        else:
            cb.ax.set_xticklabels(cbar_ticklabels)
    return pm


def plot_xy(ax, x, y, data, cmap="CN_ref", bounds=np.arange(0.5,4.6,1), cbar=True, orientation="vertical",
                cbar_ticks=None, cbar_ticklabels=None, **kwargs):
    """
    :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
    :param x: mesh grid x for data units: m
    :param y: y: mesh grid y for data units: m
    :param data: radar data ,dims like x,y
    :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
    :param orientation: vertical or horizontal, if cbar is True , this is vaild!, colorbar oriention!
    :param cbar: if True, plot with colorbar, else not!
    :param bounds: Monotonically increasing sequence of boundaries
    :param cbar_ticks: Set the locations of the tick marks from sequence ticks
    :param cbar_ticklabels: Set the text values of the tick labels.
    :param kwargs:
    :return:
    """
    assert isinstance(ax,  matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"

    #ax.set_aspect("equal")
    cmaps = plt.get_cmap(cmap)
    norm = BoundaryNorm(bounds, ncolors=cmaps.N, clip=True)
   # gci = ax.pcolormesh(x / 1000., y / 1000., data, cmap=cmaps, \
    gci = ax.pcolormesh(x , y , data, cmap=cmaps, \
                        zorder=0, norm=norm, **kwargs)

    ax.set_xlabel('Latitude (degrees)')
    ax.set_ylabel('Model Level')

    if cbar:
        cb = plt.colorbar(mappable=gci, ax=ax, orientation=orientation)
        if cbar_ticks is not None and cbar_ticklabels is not None:
            cb.set_ticks(cbar_ticks)
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)
        else:
            cb.set_ticks(bounds)
    return gci


refdata= nc.Dataset('./NEXRAD_REF3D_')
REF3d=refdata['REFMOSAIC3D'][0,:,:,:]
print(REF3d[22,390:420,230])
XLAT=refdata['XLAT_M'][0,:,:]
XLON=refdata['XLONG_M'][0,:,:]

frain = xr.open_dataset('./qr_rtv_100.nc')
fsnow = xr.open_dataset('./qs_rtv_100.nc')
fgrau = xr.open_dataset('./qg_rtv_100.nc')
wrfdata=nc.Dataset('./wrfvar_output')
rain=frain.qr_rtv[1:51,1:601,1:601]
snow=fsnow.qs_rtv[1:51,1:601,1:601]
grau=fgrau.qg_rtv[1:51,1:601,1:601]

'''calculate temperature C from perturbed potential temperature (K)'''
rp0 = 1.E-5   # p0=100000 的倒数
Rd = 287.04
Cp = 7.*Rd/2.
RCP = Rd/Cp
TEMP=wrfdata['T'][0,:,:,:]
Pres=wrfdata['P'][0,:,:,:]+wrfdata['PB'][0,:,:,:]
T= (TEMP + 300.) * ( Pres * rp0 )**RCP -273.15

qv=wrfdata['QVAPOR'][0,:,:,:]
Pres=wrfdata['P'][0,:,:,:]+wrfdata['PB'][0,:,:,:]

rho = Pres / ((T+273.15) * (eps_phy + qv) / (eps_phy *(1.+qv))) / rd   #air density
print(rho[2,200,200])
zerr= a1 * (rho * rain)**1.75
zegr= a4 * (rho * grau)**1.75
print(T[5,390:420,230])
zews =np.where( T <= 0, 0.0, a2 * (rho * snow)**1.75) # T<=-1, zews=0, T>-1. zews=...
zeds=np.where( T > 0, 0.0, a3 * (rho * snow)**1.75)

ret=np.zeros([4,50,600,600],np.float)

ret[0,:,:,:]=zerr
ret[1,:,:,:]=zeds
ret[2,:,:,:]=zews
ret[3,:,:,:]=zegr
'''
ret[0,:,:,:]=rain
ret[1,:,:,:]=snow
ret[2,:,:,:]=grau
'''
#ret_index=np.argmax(ret, axis=0)
ret_index=np.where(np.all(ret==0.0,axis=0),np.nan,np.argmax(ret, axis=0))


ret_index=np.where(REF3d<-1, np.nan, ret_index)
ret_index=np.where(ret_index == 3, 4, ret_index) #LD graupel
ret_index=np.where(ret_index == 2, 3, ret_index) #wet snow
ret_index=np.where(ret_index == 1, 2, ret_index) #dry snow
#ret_index=np.where(ret_index == 1 and T < 0, 4, ret_index)
ret_index=np.where(ret_index == 0, 1, ret_index)  #rain
#print(ret_index[22,:,194])
print(ret_index.shape)
#np.where(np.any(np.isnan(HCL), axis=0), np.nan, np.argmax(HCL, axis=0) + 1)
'''
print('start')
for k in range(0,5):
    for j in range(0,500):
        for i in range(0,500):
            if max([rain[k,j,i],snow[k,j,i],grau[k,j,i]]) > 0.0:
                if max([rain[k,j,i],snow[k,j,i],grau[k,j,i]]) == rain[k,j,i]:
                    ret[k,j,i]=2
                elif max([rain[k,j,i],snow[k,j,i],grau[k,j,i]]) == grau[k,j,i]:
                    ret[k,j,i]=7
                elif max([rain[k,j,i],snow[k,j,i],grau[k,j,i]]) == snow[k,j,i]:
                    if TT[k,j,i] >= 0 :
                        ret[k,j,i]=5
                    else:
                        ret[k,j,i]=4
print('end')
'''
x = range(0,600)
y = range(0,600)
#hcl = hca3d[20,:,:]
hcl = ret_index[19,:,:]
ticks = np.arange(1, 5, 1)
ticklabels = ['Rain', 'Dry Snow',
                           'Wet Snow', 'Graupel']
ax = plt.axes(projection=ccrs.PlateCarree())
#ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=-99.76, central_latitude=35.84))
#fig, ax = plt.subplots()

extend=[-99.2,-90,33.8,42.4] #min_lon, max_lon, min_lat, max_lat = extend
plot_lonlat_map(ax, XLON,XLAT, hcl,transform=ccrs.PlateCarree(),extend=extend,
        cmap=cmaphid, bounds=np.arange(0.5,4.6,1),
        cbar_ticks=ticks, cbar_ticklabels=ticklabels)
#plot_xy(ax, x, y, hcl, cmap=cmaphid, bounds=np.arange(0.5,10.6,1), cbar_ticks=ticks, cbar_ticklabels=ticklabels)
plt.savefig('./ret20level_723.jpg')

plt.clf()
#x = range(0,500)
y = range(0,50)
x=XLAT[:,255]
hcl = ret_index[:,:,255]
#x = XLON[200,:]
#hcl = hca3d[:,200,:]
ticks = np.arange(1, 5, 1)
ticklabels = ['Rain', 'Dry Snow',
                           'Wet Snow',  'Graupel']

fig, ax = plt.subplots()
plot_xy(ax, x, y, hcl, cmap=cmaphid, bounds=np.arange(0.5,4.6,1), cbar_ticks=ticks, cbar_ticklabels=ticklabels)
plt.savefig('./ret_vert_723.jpg')
