import HID as hid
import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib.colors as colors
import matplotlib
import xarray as xr
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


frain = xr.open_dataset('qr_rtv_100.nc')
fsnow = xr.open_dataset('qs_rtv_100.nc')
fgrau = xr.open_dataset('qg_rtv_100.nc')
rain=frain.qr_rtv[1:51,1:601,1:601]
refdata= nc.Dataset('./NEXRAD_REF3D_')
zdrfile= nc.Dataset('./zdr_100.nc')
kdpfile= nc.Dataset('./kdp_100.nc')
rhvfile= nc.Dataset('./rhv_100.nc')
#aaa=radardata.variables.keys()
wrfdata=nc.Dataset('./wrfvar_output')

REF3d=refdata['REFMOSAIC3D']
XLAT=refdata['XLAT_M'][0,:,:]
XLON=refdata['XLONG_M'][0,:,:]
height=refdata['height']
ZDR3d=zdrfile['zdr'][1:51,1:601,1:601]  #last no not included
KDP3d=kdpfile['kdp'][1:51,1:601,1:601]
RHV3d=rhvfile['rhv'][1:51,1:601,1:601]
#print(XLON[200,194:205])

'''calculate temperature C from perturbed potential temperature (K)'''
rp0 = 1.E-5   # p0=100000 的倒数
Rd = 287.04
Cp = 7.*Rd/2.
RCP = Rd/Cp
TEMP=wrfdata['T'][0,:,:,:]
Pres=wrfdata['P'][0,:,:,:]+wrfdata['PB'][0,:,:,:]
T= (TEMP + 300.) * ( Pres * rp0 )**RCP -273.15
#print(T[22,200,:])

hca3d=[]
for i in range(0,50):
    REF = REF3d[0, i, :, :]
    ZDR = ZDR3d[i, :, :]
    KDP = KDP3d[i, :, :]
    RHV = RHV3d[i, :, :]
    TT  = T[i,:,:]
 #   print(ZDR.shape, KDP.shape, RHV.shape, REF.shape)
    hca_product=hid.fhc_HCL(dBZ=REF, ZDR=ZDR, KDP=KDP, CC=RHV, method="hybrid",T=TT,
            band="S")
    hca3d=np.append(hca3d,hca_product)

hca3d=hca3d.reshape(50,600,600)
hca3d=np.where(hca3d==1, np.nan , hca3d)

'''write HCA product to nc file'''
level=range(0,50)
foo=xr.DataArray(hca3d,dims=['level','XLAT','XLON'])
foo.to_netcdf('hcatest.nc','w')

refdata.close()
zdrfile.close()
kdpfile.close()
rhvfile.close()
wrfdata.close()

'''plot 0.0'''
hid_colors = ['LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
              'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
cmaphid = colors.ListedColormap(hid_colors)

def plot_xy(ax, x, y, data, cmap="CN_ref", bounds=np.arange(0.5,10.6,1), cbar=True, orientation="vertical",
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


def plot_lonlat_map(ax, lon, lat, data, transform, extend=None, cmap="CN_ref", bounds=np.arange(0.5,10.6,1),
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


x = range(0,600)
y = range(0,600)
#hcl = hca3d[20,:,:]
hcl = hca3d[19,:,:]
ticks = np.arange(1, 11, 1)
ticklabels = ['Drizzle', 'Rain', 'Ice Crystals', 'Dry Snow',
                           'Wet Snow', 'Vertical Ice', 'LD Graupel',
                           'HD Graupel', 'Hail', 'Big Drops']
ax = plt.axes(projection=ccrs.PlateCarree())
#ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=-99.76, central_latitude=35.84))
#fig, ax = plt.subplots()

extend=[-99.2,-90,33.8,42.4] #min_lon, max_lon, min_lat, max_lat = extend
plot_lonlat_map(ax, XLON, XLAT, hcl,transform=ccrs.PlateCarree(),extend=extend,
        cmap=cmaphid, bounds=np.arange(0.5,10.6,1),
        cbar_ticks=ticks, cbar_ticklabels=ticklabels)
#plot_xy(ax, x, y, hcl, cmap=cmaphid, bounds=np.arange(0.5,10.6,1), cbar_ticks=ticks, cbar_ticklabels=ticklabels)
plt.savefig('./hca30level.jpg')

plt.clf()
#x = range(0,500)
y = range(0,50)
x=XLAT[:,255]
hcl = hca3d[:,:,255]
#x = XLON[200,:]
#hcl = hca3d[:,200,:]
ticks = np.arange(1, 11, 1)
ticklabels = ['Drizzle', 'Rain', 'Ice Crystals', 'Dry Snow',
                           'Wet Snow', 'Vertical Ice', 'LD Graupel',
                           'HD Graupel', 'Hail', 'Big Drops']

fig, ax = plt.subplots()
plot_xy(ax, x, y, hcl, cmap=cmaphid, bounds=np.arange(0.5,10.6,1), cbar_ticks=ticks, cbar_ticklabels=ticklabels)
plt.savefig('./hca.jpg')
