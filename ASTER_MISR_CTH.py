from osgeo import gdal, osr
import numpy as np
import os, glob, sys, getopt, argparse
import MisrToolkit as Mtk
import pandas as pd
from pyresample import create_area_def
import ASTER_cloud_mask_funcs as Acm

def getaster1tgeo(file_name,resolution):
    # print('Processing File: ' + file_name + ' (' + str(k+1) + ' out of ' 
    # + str(len(file_list)) + ')')
    # Read in the file and metadata
    aster = gdal.Open(file_name)
    # aster_sds = aster.GetSubDatasets()
    meta = aster.GetMetadata()
    # Define UL, LR, UTM zone    
    ll = [float(x) for x in meta['LOWERLEFTM'].split(', ')]
    ur = [float(x) for x in meta['UPPERRIGHTM'].split(', ')]
    utm = int(meta['UTMZONENUMBER'])
    # Create UTM zone code numbers    
    #n_s = np.float(meta['NORTHBOUNDINGCOORDINATE'])
    # Define UTM zone based on North or South
    # if n_s < 0:
    #     ur[1] = ur[1] + 10000000
    #     ll[1] = ll[1] + 10000000  
    # Define extent for UTM North zones             
    proj_dict = {'proj': 'utm','zone': utm, 'datum':'WGS84', \
                 'towgs84':'0,0,0','ellps': 'WGS84', 'no_defs':'', 'units': 'm'} 
    area_extent = (ll[1]-7.5,ll[0]-7.5, ur[1]+7.5, ur[0]+7.5)
    area_def = create_area_def('aster', proj_dict, area_extent=area_extent,resolution=resolution,units='m')
    lons, lats = area_def.get_lonlats()
    return lats, lons

def getbls(asterfile):
    ast_resolution = 15
    ast_lats, ast_lons = getaster1tgeo(asterfile, 15)
    aster = gdal.Open(asterfile)
    meta = aster.GetMetadata()
    orbit = meta['ORBITNUMBER']

    misrfile = glob.glob('/data/gdi/e/MISR/TCCloud/'+yy+'.'+mm+'.'+dd+'/*'+orbit+'*.hdf')[0]
    m = Mtk.MtkFile(misrfile)
    path = m.path
    misr_resolution = 1100
    block, line, sample = Mtk.latlon_to_bls(path, misr_resolution, ast_lats.flatten(), ast_lons.flatten())

    title = asterfile.split('/')[-1].split('.hdf')[0]

    line = line + 0.5
    if(block[0] != block[-1]):
        block2 = np.where(block==block[-1])
        line[block2] = line[block2] + 128
    line = line.astype(int)
    line = line.reshape(ast_lats.shape[0], ast_lats.shape[1])

    sample = sample + 0.5
    sample = sample.astype(int) + 1
    sample = sample.reshape(ast_lats.shape[0], ast_lats.shape[1])

    block = block.reshape(ast_lats.shape[0], ast_lats.shape[1])

    return block,line,sample


def getMISR_CTH_hist(asterfile):
    ast_resolution = 15
    ast_lats, ast_lons = getaster1tgeo(asterfile, 15)
    aster = gdal.Open(asterfile)
    meta = aster.GetMetadata()
    orbit = meta['ORBITNUMBER']

    block, line, sample = getbls(asterfile)

    date = asterfile.split('/')[-2]

    misrfile = glob.glob('/data/gdi/e/MISR/TCCloud/'+date+'/*'+orbit+'*.hdf')[0]
    m = Mtk.MtkFile(misrfile)
    path = m.path
    misr_resolution = 1100

    b1 = int(block[0])
    b2 = int(block[-1])

    b3N = Acm.get_dn(asterfile, '3N')
    b3N[b3N==0.]=np.nan

    blk = block.reshape(line.shape[0], line.shape[1]).astype(float)
    row = line.astype(float)
    col = sample.astype(float)

    blk[np.isnan(b3N)]=np.nan
    row[np.isnan(b3N)]=np.nan
    col[np.isnan(b3N)]=np.nan

    r1 = Mtk.MtkRegion(path, b1, b1)
    h1 = m.grid('Stereo_1.1_km').field('CloudTopHeight').read(r1).data()

    row1 = row[blk==b1]
    col1 = col[blk==b1]

    ls1 = np.hstack((np.array([row1]).T, np.array([col1]).T))

    uniquels1 = np.unique(ls1, axis=0)

    x1 = uniquels1[:,0].astype(int)
    y1 = uniquels1[:,1].astype(int)

    MISR_CTH1 = h1[x1,y1]

    if b1 != b2:
        r2 = Mtk.MtkRegion(path, b2, b2)
        h2 = m.grid('Stereo_1.1_km').field('CloudTopHeight').read(r2).data()

        row2 = row[blk==b2]
        col2 = col[blk==b2]

        ls2 = np.hstack((np.array([row2]).T, np.array([col2]).T))

        uniquels2 = np.unique(ls2, axis=0)

        x2 = uniquels2[:,0].astype(int) - 128
        y2 = uniquels2[:,1].astype(int)

        MISR_CTH2 = h2[x2,y2]
        MISR_CTH_all = np.concatenate((MISR_CTH1, MISR_CTH2))
    else:
        MISR_CTH_all = MISR_CTH1

    binrange = np.arange(0,10100,100)
    histall, _ = np.histogram(MISR_CTH_all[MISR_CTH_all>0], bins=binrange)

    return histall


out_dir = '/data/gdi/c/mdevera2/ASTER_MISR/'

all_file_list = pd.read_csv('/data/keeling/a/mdevera2/ASTER_cloud_mask/ocean_no_cirrus_glint2.csv')
all_file_list['desc'][all_file_list['desc'].isnull()]='good'
good_files = all_file_list[all_file_list['desc']=='good']
sun_glint_files = all_file_list[(all_file_list['desc']=='sun_glint') & (all_file_list['add_desc']=='okay')]

file_list = pd.concat([good_files, sun_glint_files])

print(len(file_list))

i = 1

heighthistall = np.zeros(100)

for asterfile in file_list['filename']:
    print(i)
    print(asterfile)
    hall = getMISR_CTH_hist(asterfile)
    heighthistall += hall
    
    i = i + 1

fig, ax = plt.subplots()
all = np.sum(heighthistall)
bincenter = np.arange(50,10050,100)
plt.plot(heighthistall/all, bincenter, label='all')

plt.xlim(0.001,0.2)
plt.xscale('log')
plt.xlabel('Normalized Frequency')
plt.ylabel('Cloud Height (m)')
plt.legend(title='Cloud Diameter (km)')
plt.savefig('/data/keeling/a/mdevera2/macrost/cth_misr_dist_test_log.png', dpi=300)
plt.clf()