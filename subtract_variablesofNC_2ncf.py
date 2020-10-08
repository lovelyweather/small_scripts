import xarray as xr
import netCDF4 as nc
import sys

wrfdata=xr.open_dataset(sys.argv[1])
T=wrfdata['T'][0,:,:,:]
P=wrfdata['P'][0,:,:,:]
PB=wrfdata['PB'][0,:,:,:]
QVAPOR=wrfdata['QVAPOR'][0,:,:,:]
ds=xr.Dataset({'T':T,'P':P,'PB':PB,'QVAPOR':QVAPOR})
ds.to_netcdf(sys.argv[2])
