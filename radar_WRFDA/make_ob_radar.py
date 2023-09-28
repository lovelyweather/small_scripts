import os, yaml
import warnings
warnings.filterwarnings('ignore') # setting ignore as a parameter

import numpy as np
import pyart
from pyart.correct import dealias_region_based
from netCDF4 import num2date
from datetime import date, datetime,timedelta

import func_inc
from func_inc import Radar_prep, height_exist

def make_obs_radar(yaml_data:dict):

    '''setting for domain and time'''
    station, extend_ll, resolution  = yaml_data['station'], yaml_data['extend_ll'], \
        yaml_data['resolution']

    time_ana   = yaml_data['time_ana']
    time_begin = time_ana - timedelta(minutes=yaml_data['time_window'])
    time_end   = time_ana + timedelta(minutes=yaml_data['time_window'])
    #time_begin = datetime.strptime(str(yaml_data['time_begin']),"%Y%m%d%H%M" )
    #time_end   = datetime.strptime(str(yaml_data['time_end']),"%Y%m%d%H%M" )


    Min_Lat, Max_Lat = extend_ll[0], extend_ll[1]
    Min_Lon, Max_Lon = extend_ll[2], extend_ll[3]     
    Lat_des_1D = np.arange( Max_Lat, Min_Lat, -1*resolution )  # 生成插值后的纬度
    Lon_des_1D = np.arange( Min_Lon, Max_Lon,  resolution )  # 生成插值后的经度

    '''some customized setting for output'''
    radar_type='WSR98D'
    fileOut = 'ob.radar'+time_ana.strftime("%Y%m%d%H%M")
    dbz_qc, vel_qc = 1, 1
    dbz_err, vel_err = 2, 1

    '''select the existing stations'''

    filenames = []
    cur_time = time_begin
    while cur_time <= time_end:
        for i_radar in station:
            i_file = 'qc_Z_RADR_I_Z' + i_radar +  '_' + cur_time.strftime("%Y%m%d%H%M%S") \
                + '_O_DOR_SAD_CAP_FMT.bin.nc'
            if os.path.exists(i_file):
                filenames.append(i_file)
            else:
                print(f"The file '{i_file}' does not exist.")
        cur_time += timedelta(minutes=6)

    n_radar = len(filenames)

    with open(fileOut, 'w') as f:
    # formatted print. Use {} and : to replace % in Fortran. 
        f.write("total_number= {:3d}\n".format(n_radar))
        f.write('#------------------------#\n')

    for i_radar in filenames:

        # 1. read radar file
        try: 
            radar = pyart.io.read(i_radar)
            print('Performing operations for file:'+i_radar)
        except FileNotFoundError:
            print(f"Error: The file '{i_radar}' does not exist.")

        times = int(radar.time["data"][0])
        units = radar.time["units"]
        calendar = radar.time["calendar"]
        obs_date = num2date(
                times,
                units,
                calendar,
                only_use_cftime_datetimes=False,
                only_use_python_datetimes=True,
            )
        date_str = datetime.fromtimestamp(obs_date.timestamp()).strftime("%Y_%m_%d_%H:%M:%S") #2023_09_02_00:00:00

        # reorder the azimuths
        radar = Radar_prep.ordered_az(radar, 'dBZ')
        radar = Radar_prep.ordered_az(radar, 'VEL')
        radar = Radar_prep.ordered_az(radar, 'lat')
        radar = Radar_prep.ordered_az(radar, 'lon')
        radar = Radar_prep.ordered_az(radar, 'alt')

        lat = radar.fields['lat_reorder']['data']
        lon = radar.fields['lon_reorder']['data']
        alt = radar.fields['alt_reorder']['data'] + radar.altitude['data'][0]

        nlevel_ValidValue = np.zeros((len(Lat_des_1D), len(Lon_des_1D))).astype(int)
        maxlev = 20
        values = np.zeros((3, maxlev, len(Lat_des_1D), len(Lon_des_1D))) #第一维高度，第二维反射率因子，第三维速度
        values[:, :, :,:] = -999.0

        # lat & lon before, size is (nlat*nlon,2)
        LatLonAlt_Before = np.hstack(
            (lat.reshape(-1, 1), lon.reshape(-1, 1), alt.reshape(-1, 1)) ) # 按水平方向进行叠加，形成两列
        refl_stack = radar.fields['dBZ_reorder']['data'].reshape(-1, 1)
        vel_stack  = radar.fields['VEL_reorder']['data'].reshape(-1, 1)


        # 2. interpolate to equal LatLon
        n_count = 0
        for i in range(0,len(LatLonAlt_Before)):
            lat_point, lon_point, alt_point = LatLonAlt_Before[i]
            if lat_point == 0.0 or lon_point == 0.0:
                    continue
            if refl_stack[i] >= 0 and abs(vel_stack[i] < 100): 
                i_loc = int((Max_Lat - lat_point) / resolution) #row i
                j_loc = int((lon_point - Min_Lon) / resolution) #column j

                if not height_exist(alt_point, values[0,:, i_loc, j_loc], 500):
                    level_toPut = nlevel_ValidValue[i_loc, j_loc]
                    #print(level_toPut)
                    if level_toPut >= 10:
                        print(level_toPut, alt_point, refl_stack[i], vel_stack[i])
                    values[:, level_toPut, i_loc, j_loc] = [alt_point, refl_stack[i][0], vel_stack[i][0]]
                    nlevel_ValidValue[i_loc, j_loc] += 1 
                    n_count += 1


        # 3. write to file with WRFDA format
        with open(fileOut, 'a') as f:
            f.write('\n')
            # 5d:5个字符 12d:12个字符 

            f.write("{:>5}  {:>12}{:8.3f}  {:8.3f}  {:8.1f}  {:>19}{:6}{:6}\n".format(radar_type, i_radar[13:17], radar.latitude['data'][0], \
                            radar.longitude['data'][0], radar.altitude['data'][0],\
                                date_str, n_count, maxlev ))
            
            for i in range(len(Lat_des_1D)):
                for j in range(len(Lon_des_1D)):
                    if nlevel_ValidValue[i,j] > 0:
                        f.write("FM-128 RADAR   {:>19}  {:12.3f}  {:12.3f}  {:8.1f}  {:>6}\n"\
                                    .format(date_str, Lat_des_1D[i], Lon_des_1D[j], \
                                    radar.altitude['data'][0], nlevel_ValidValue[i,j])) 
                        for k in range(nlevel_ValidValue[i,j]):
                            f.write("   {:12.1f}{:12.3f}{:4}{:12.3f}  {:12.3f}{:4}{:12.3f}\n"\
                                    .format(values[0, k, i_loc, j_loc], \
                                    values[1, k, i_loc, j_loc],dbz_qc, dbz_err,\
                                    values[2, k, i_loc, j_loc],vel_qc, vel_err ))
                            
    print(f"OBS for {time_ana}  is done ^_^ ")

if __name__=='__main__':

    '''parse the yaml config file'''
    yaml_file = 'ob_radar.yaml'
    with open(yaml_file, 'r') as file:
        yaml_data = yaml.safe_load(file)
    
    time_begin = datetime.strptime(str(yaml_data['time_begin']),"%Y%m%d%H%M" )
    time_end   = datetime.strptime(str(yaml_data['time_end']),"%Y%m%d%H%M" )

    time_ana   = time_begin
    while time_ana <= time_end:
        yaml_data['time_ana'] = time_ana
        make_obs_radar(yaml_data)
        time_ana += timedelta(minutes=yaml_data['time_intval'])