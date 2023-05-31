import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import xarray as xr
from osgeo import gdal
from pykdtree.kdtree import KDTree
import metpy.calc
from metpy.units import units
import climlab

def hour_rounder(t):
    # Rounds to nearest hour by adding a timedelta hour if minute >= 30
    return (t.replace(second=0, microsecond=0, minute=0, hour=t.hour)+timedelta(hours=t.minute//30))

def get_data(file, ERA5single, ERA5press, ERA5press2):
    date = datetime.strptime(file.split('/')[-1].split('_')[2][3:], '%m%d%Y%H%M%S')
    dt64 = np.datetime64(hour_rounder(date))

    aster = gdal.Open(file)
    meta = aster.GetMetadata()
    lats = []
    lats.append(np.float32(meta['UPPERLEFT'].split(',')[0]))
    lats.append(np.float32(meta['LOWERRIGHT'].split(',')[0]))
    lons = []
    lons.append(np.float32(meta['UPPERLEFT'].split(',')[1]))
    lons.append(np.float32(meta['LOWERRIGHT'].split(',')[1]))
    
    tree = KDTree(ERA5single['latitude'].data)
    dist, ind_lat = tree.query(np.array(lats))
    tree = KDTree(ERA5single['longitude'].data)
    dist, ind_lon = tree.query(np.array(lons))

    t = np.where(ERA5single['time'].data == dt64)[0]
    sst = np.nanmean(ERA5single['sst'].data[t,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    tcwv = np.nanmean(ERA5single['tcwv'].data[t,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    u10 = units.Quantity(ERA5single['u10'].data[t,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v10 = units.Quantity(ERA5single['v10'].data[t,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    uv10 = np.nanmean(metpy.calc.wind_speed(u10, v10)).magnitude
    wdir10 = np.nanmean(metpy.calc.wind_direction(u10, v10)).magnitude
    rh1000 = np.nanmean(ERA5press['r'].data[t,-1,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    rh975 = np.nanmean(ERA5press['r'].data[t,-2,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    rh950 = np.nanmean(ERA5press['r'].data[t,-3,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    rh925= np.nanmean(ERA5press['r'].data[t,-4,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    rh900= np.nanmean(ERA5press['r'].data[t,-5,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    rh875= np.nanmean(ERA5press['r'].data[t,-6,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    rh850= np.nanmean(ERA5press['r'].data[t,-7,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    w1000 = np.nanmean(ERA5press2['w'].data[t,-1,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    w975 = np.nanmean(ERA5press2['w'].data[t,-2,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    w950 = np.nanmean(ERA5press2['w'].data[t,-3,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    w925 = np.nanmean(ERA5press2['w'].data[t,-4,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    w900 = np.nanmean(ERA5press2['w'].data[t,-5,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    w875 = np.nanmean(ERA5press2['w'].data[t,-6,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    w850 = np.nanmean(ERA5press2['w'].data[t,-7,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1])
    u1000 = units.Quantity(ERA5press2['u'].data[t,-1,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    u975 = units.Quantity(ERA5press2['u'].data[t,-2,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    u950 = units.Quantity(ERA5press2['u'].data[t,-3,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    u925 = units.Quantity(ERA5press2['u'].data[t,-4,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    u900 = units.Quantity(ERA5press2['u'].data[t,-5,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    u875 = units.Quantity(ERA5press2['u'].data[t,-6,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    u850 = units.Quantity(ERA5press2['u'].data[t,-7,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v1000 = units.Quantity(ERA5press2['v'].data[t,-1,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v975 = units.Quantity(ERA5press2['v'].data[t,-2,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v950 = units.Quantity(ERA5press2['v'].data[t,-3,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v925 = units.Quantity(ERA5press2['v'].data[t,-4,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v900 = units.Quantity(ERA5press2['v'].data[t,-5,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v875 = units.Quantity(ERA5press2['v'].data[t,-6,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v850 = units.Quantity(ERA5press2['v'].data[t,-7,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    uv1000 = np.nanmean(metpy.calc.wind_speed(u1000, v1000)).magnitude
    uv975 = np.nanmean(metpy.calc.wind_speed(u975, v975)).magnitude
    uv950 = np.nanmean(metpy.calc.wind_speed(u950, v950)).magnitude
    uv925 = np.nanmean(metpy.calc.wind_speed(u925, v925)).magnitude
    uv900 = np.nanmean(metpy.calc.wind_speed(u900, v900)).magnitude
    uv875 = np.nanmean(metpy.calc.wind_speed(u875, v875)).magnitude
    uv850 = np.nanmean(metpy.calc.wind_speed(u850, v850)).magnitude
    wdir1000 = np.nanmean(metpy.calc.wind_direction(u1000, v1000)).magnitude
    wdir975 = np.nanmean(metpy.calc.wind_direction(u975, v975)).magnitude
    wdir950 = np.nanmean(metpy.calc.wind_direction(u950, v950)).magnitude
    wdir925 = np.nanmean(metpy.calc.wind_direction(u925, v925)).magnitude
    wdir900 = np.nanmean(metpy.calc.wind_direction(u900, v900)).magnitude
    wdir875 = np.nanmean(metpy.calc.wind_direction(u875, v875)).magnitude
    wdir850 = np.nanmean(metpy.calc.wind_direction(u850, v850)).magnitude
    u600 = units.Quantity(ERA5press2['u'].data[t,-14,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    v600 = units.Quantity(ERA5press2['v'].data[t,-14,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1], 'm/s')
    uv600 = np.nanmean(metpy.calc.wind_speed(u600, v600)).magnitude
    wdir600 = np.nanmean(metpy.calc.wind_direction(u600, v600)).magnitude
    T1000 = ERA5press['t'].data[t,-1,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1] * units.kelvin
    T975 = ERA5press['t'].data[t,-2,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1] * units.kelvin
    T950 = ERA5press['t'].data[t,-3,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1] * units.kelvin
    T925 = ERA5press['t'].data[t,-4,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1] * units.kelvin
    T900 = ERA5press['t'].data[t,-5,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1] * units.kelvin
    T875 = ERA5press['t'].data[t,-6,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1] * units.kelvin
    T850 = ERA5press['t'].data[t,-7,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1] * units.kelvin
    T700 = ERA5press['t'].data[t,-12,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1] * units.kelvin
    theta1000 = metpy.calc.potential_temperature(1000 * units.mbar, T1000)
    theta975 = metpy.calc.potential_temperature(975 * units.mbar, T975)
    theta950 = metpy.calc.potential_temperature(950 * units.mbar, T950)
    theta925 = metpy.calc.potential_temperature(925 * units.mbar, T925)
    theta900 = metpy.calc.potential_temperature(900 * units.mbar, T900)
    theta875 = metpy.calc.potential_temperature(875 * units.mbar, T875)
    theta850 = metpy.calc.potential_temperature(1000 * units.mbar, T850)
    theta700 = metpy.calc.potential_temperature(700 * units.mbar, T700)
    LTS700 = np.nanmean((theta700 - theta1000)).magnitude
    LTS850 = np.nanmean((theta850 - theta1000)).magnitude
    LTS875 = np.nanmean((theta875 - theta1000)).magnitude
    LTS900 = np.nanmean((theta900 - theta1000)).magnitude
    LTS925 = np.nanmean((theta925 - theta1000)).magnitude
    LTS950 = np.nanmean((theta950 - theta1000)).magnitude
    LTS975 = np.nanmean((theta975 - theta1000)).magnitude
    EIS700 = np.nanmean(climlab.utils.thermo.EIS(T1000/units.kelvin, T700/units.kelvin)).magnitude
    thetae1000 = np.nanmean(metpy.calc.equivalent_potential_temperature(1000*units.hPa, T1000, metpy.calc.dewpoint_from_relative_humidity(T1000, rh1000/100))).magnitude
    thetae975 = np.nanmean(metpy.calc.equivalent_potential_temperature(975*units.hPa, T975, metpy.calc.dewpoint_from_relative_humidity(T975, rh975/100))).magnitude
    thetae950 = np.nanmean(metpy.calc.equivalent_potential_temperature(950*units.hPa, T950, metpy.calc.dewpoint_from_relative_humidity(T950, rh950/100))).magnitude
    thetae925 = np.nanmean(metpy.calc.equivalent_potential_temperature(925*units.hPa, T925, metpy.calc.dewpoint_from_relative_humidity(T925, rh925/100))).magnitude
    thetae900 = np.nanmean(metpy.calc.equivalent_potential_temperature(900*units.hPa, T900, metpy.calc.dewpoint_from_relative_humidity(T900, rh900/100))).magnitude
    thetae875 = np.nanmean(metpy.calc.equivalent_potential_temperature(875*units.hPa, T875, metpy.calc.dewpoint_from_relative_humidity(T875, rh875/100))).magnitude
    thetae850 = np.nanmean(metpy.calc.equivalent_potential_temperature(850*units.hPa, T850, metpy.calc.dewpoint_from_relative_humidity(T850, rh850/100))).magnitude

    p = np.flip(ERA5press['level'].data[-7:]*units.hPa)
    T = np.flip(np.nanmean(ERA5press['t'].data[t,-7:,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1],axis=(2,3))*units.K)[0]
    rh = np.flip(np.nanmean(ERA5press['r'].data[t,-7:,ind_lat[0]:ind_lat[1]+1,ind_lon[0]:ind_lon[1]+1],axis=(2,3))/100)[0]
    Td = metpy.calc.dewpoint_from_relative_humidity(T, rh)
    prof = metpy.calc.parcel_profile(p, T[0], Td[0]).to('degC')
    CAPE, CIN = metpy.calc.cape_cin(p, T, Td, prof)

    return([file.split('/')[-1].split('.')[0], file.split('/')[-2], sst, tcwv, np.nanmean(u10).magnitude, np.nanmean(v10).magnitude, 
    uv10, wdir10, rh1000, rh975, rh950, rh925, rh900, rh875, rh850, w1000, w975, w950, w925, w900, w875, w850, 
    np.nanmean(u1000).magnitude, np.nanmean(u975).magnitude, np.nanmean(u950).magnitude, np.nanmean(u925).magnitude, np.nanmean(u900).magnitude, np.nanmean(u875).magnitude, np.nanmean(u850).magnitude, 
    np.nanmean(v1000).magnitude, np.nanmean(v975).magnitude, np.nanmean(v950).magnitude, np.nanmean(v925).magnitude, np.nanmean(v900).magnitude, np.nanmean(v875).magnitude, np.nanmean(v850).magnitude, 
    uv1000, uv975, uv950, uv925, uv900, uv875, uv850, wdir1000, wdir975, wdir950, wdir925, wdir900, wdir875, wdir850, np.nanmean(u600).magnitude, np.nanmean(v600).magnitude, uv600, wdir600, 
    LTS700, LTS850, LTS875, LTS900, LTS925, LTS950, LTS975, EIS700,
    thetae1000, thetae975, thetae950, thetae925, thetae900, thetae875, thetae850, CAPE.magnitude, CIN.magnitude])

all_file_list = pd.read_csv('/data/keeling/a/mdevera2/ASTER_cloud_mask/ocean_no_cirrus_glint2.csv')
all_file_list['desc'][all_file_list['desc'].isnull()] = 'good'
print(all_file_list)

description = ['good', 'sun_glint']
sun_add = ['okay', 'maybe']

out_dir = '/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/ERA5_data/'

CAMP2Expressure = xr.open_dataset('/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/ERA5_data/CAMP2Ex_RH_T_SH.nc')
CAMP2Expressure2 = xr.open_dataset('/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/ERA5_data/CAMP2Ex_w_U_V.nc')
CAMP2Exsingle = xr.open_dataset('/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/ERA5_data/CAMP2Ex_singlelevels.nc')

data_list = []

for desc in description:
    if desc == 'good':
        files = all_file_list[all_file_list['desc']==desc].copy().reset_index(drop=True)
        print(files)
        for file in files['filename']:
            print(file)
            data_list.append(get_data(file, CAMP2Exsingle, CAMP2Expressure, CAMP2Expressure2))

    else:
        files = all_file_list[(all_file_list['desc']==desc) & (all_file_list['add_desc']==sun_add[0])].copy().reset_index(drop=True)
        for file in files['filename']:
            print(file)
            data_list.append(get_data(file, CAMP2Exsingle, CAMP2Expressure, CAMP2Expressure2))

df = pd.DataFrame(data_list,columns=['file','day','sst','tcwv','u10','v10','uv10','wdir10','rh1000','rh975','rh950','rh925','rh900','rh875','rh850',
'w1000','w975','w950','w925','w900','w875','w850','u1000','u975','u950','u925','u900','u875','u850','v1000','v975','v950','v925','v900','v875','v850',
'uv1000','uv975','uv950','uv925','uv900','uv875','uv850','wdir1000','wdir975','wdir950','wdir925','wdir900','wdir875','wdir850','u600', 'v600', 'uv600', 'wdir600', 
'LTS700', 'LTS850', 'LTS875', 'LTS900', 'LTS925', 'LTS950', 'LTS975', 'EIS700',
'thetae1000','thetae975','thetae950','thetae925','thetae900','thetae875','thetae850', 'CAPE850', 'CIN850'])
print(df)
df.to_csv(out_dir+'/'+'CAMP2Ex_ASTER_ERA5_data.csv', index=False)
