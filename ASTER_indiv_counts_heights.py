import aster_subroutines2 as ast
import netCDF4 as nc
import numpy as np
import collections
import pandas as pd
import os, glob, sys, getopt, argparse, re
import matplotlib.pyplot as plt
from osgeo import gdal, osr
import ASTER_cloud_mask_funcs as Acm
import MisrToolkit as Mtk

def get_counts(ds, num):
    totalpixels = np.count_nonzero(~np.isnan(ds['lowresmask'][:]))
    outmask, totalflag, number = ast.flag_each_cloud(ds['lowresmask'][:])
    onns, nnpair = ast.nns(outmask,totalflag)
    edge = ast.perimeter(outmask,totalflag)
    c = collections.Counter(outmask.flatten())
    del c[0]
    cloudypixels = sum(c.values())
    cloudfraction = cloudypixels/totalpixels

    df = pd.DataFrame(sorted(c.items()), columns=['cloud', 'numpixels'])
    df['nndist'] = onns
    df['nnpair'] = nnpair
    df['numedge'] = edge
    df['diam'] = np.sqrt(4 * df['numpixels'] * 90 * 90 / np.pi)/1000

    bins_diam = np.arange(0,7,0.1)
    diamhist, _ = np.histogram(df['diam'], bins=bins_diam)

    binrange = np.arange(0,10100,100)
    #binrange = np.arange(0,10250,250)

    height = list(ds.variables.keys())

    heightdata = ds[height[num]][:]

    heightbin0 = df['cloud'][df['diam']<0.5].values
    bin0 = heightdata[np.where(np.isin(outmask, heightbin0))]
    hist0, _ = np.histogram(bin0, bins=binrange)

    heightbin1 = df['cloud'][df['diam']<1].values
    bin1 = heightdata[np.where(np.isin(outmask, heightbin1))]
    hist1, _ = np.histogram(bin1, bins=binrange)
        
    heightbin2 = df['cloud'][df['diam']<2].values
    bin2 = heightdata[np.where(np.isin(outmask, heightbin2))]
    hist2, _ = np.histogram(bin2, bins=binrange)

    heightbin3 = df['cloud'][df['diam']<3].values
    bin3 = heightdata[np.where(np.isin(outmask, heightbin3))]
    hist3, _ = np.histogram(bin3, bins=binrange)

    heightbin4 = df['cloud'][df['diam']<4].values
    bin4 = heightdata[np.where(np.isin(outmask, heightbin4))]
    hist4, _ = np.histogram(bin4, bins=binrange)

    heightbinall = df['cloud'].values
    binall = heightdata[np.where(np.isin(outmask, heightbinall))]
    histall, _ = np.histogram(binall, bins=binrange)

    return hist0, hist1, hist2, hist3, hist4, histall, diamhist

def getMISR_CTH_hist(file):
    fileID = file.split('/')[-1].split('.')[0][:-8]
    mm = fileID.split('_')[2][3:5]
    dd = fileID.split('_')[2][5:7]
    yy = fileID.split('_')[2][7:11]
    date = yy+'.'+mm+'.'+dd
    asterfile = '/data/gdi/e/ASTER/camp2ex/L1T/'+date+'/'+fileID+'.hdf'
    print(asterfile)

    aster = gdal.Open(asterfile)
    meta = aster.GetMetadata()
    orbit = meta['ORBITNUMBER']

    b3N = Acm.get_dn(asterfile, '3N')
    b3N[b3N==0.]=np.nan

    linefile = '/data/keeling/a/mdevera2/c/ASTER_MISR/line/'+fileID+'_line.bin'
    l = open(linefile, "rb")
    line = np.fromfile(l, dtype=int)
    line = line.reshape(b3N.shape[0], b3N.shape[1])

    samplefile = '/data/keeling/a/mdevera2/c/ASTER_MISR/sample/'+fileID+'_sample.bin'
    s = open(samplefile, "rb")
    sample = np.fromfile(s, dtype=int)
    sample = sample.reshape(b3N.shape[0], b3N.shape[1])

    blockfile = '/data/keeling/a/mdevera2/c/ASTER_MISR/blocks/'+orbit+'_'+fileID+'_block.txt'
    with open(blockfile) as f:
        blocks = f.readlines()
    b1 = int(blocks[0].split(' ')[0])
    b2 = int(blocks[0].split(' ')[1])

    misrfile = glob.glob('/data/gdi/e/MISR/TCCloud/'+date+'/*'+orbit+'*.hdf')[0]
    m = Mtk.MtkFile(misrfile)
    path = m.path
    misr_resolution = 1100

    #blk = block.reshape(line.shape[0], line.shape[1]).astype(float)
    row = line.astype(float) - 1
    col = sample.astype(float) - 1

    #blk[np.isnan(b3N)]=np.nan
    row[np.isnan(b3N)]=np.nan
    col[np.isnan(b3N)]=np.nan

    r1 = Mtk.MtkRegion(path, b1, b1)
    h1 = m.grid('Stereo_1.1_km').field('CloudTopHeight').read(r1).data()

    row1 = row[row < 128]
    col1 = col[row < 128]

    ls1 = np.hstack((np.array([row1]).T, np.array([col1]).T))

    uniquels1 = np.unique(ls1, axis=0)

    x1 = uniquels1[:,0].astype(int)
    y1 = uniquels1[:,1].astype(int)

    MISR_CTH1 = h1[x1,y1]

    if b1 != b2:
        r2 = Mtk.MtkRegion(path, b2, b2)
        h2 = m.grid('Stereo_1.1_km').field('CloudTopHeight').read(r2).data()

        row2 = row[row > 127]
        col2 = col[row > 127]

        ls2 = np.hstack((np.array([row2]).T, np.array([col2]).T))

        uniquels2 = np.unique(ls2, axis=0)

        x2 = uniquels2[:,0].astype(int) - 128
        y2 = uniquels2[:,1].astype(int)

        MISR_CTH2 = h2[x2,y2]
        MISR_CTH_all = np.concatenate((MISR_CTH1, MISR_CTH2))
    else:
        MISR_CTH_all = MISR_CTH1

    print(MISR_CTH_all)
    np.savetxt('/data/keeling/a/mdevera2/c/ASTER_macrostats/ocean_no_cirrus_scenes2/MISR_heights/'+fileID+'.txt', MISR_CTH_all[MISR_CTH_all>0], delimiter =", ")

    binrange = np.arange(0,10100,100)
    #binrange = np.arange(0,10250,250)
    histall, _ = np.histogram(MISR_CTH_all[MISR_CTH_all>0], bins=binrange)

    return histall

height_file_dir = '/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/heights_w_interpolated_min_p1K/'
#out_dir = '/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/heights_w_interpolated_min/'

description = ['good', 'sun_glint']
sun_add = ['okay']

i = 1

heighthist0 = np.zeros(100)
heighthist1 = np.zeros(100)
heighthist2 = np.zeros(100)
heighthist3 = np.zeros(100)
heighthist4 = np.zeros(100)
heighthistall = np.zeros(100)
MISRheighthistall = np.zeros(100)
diameterhist = np.zeros(69)

'''heighthist0 = np.zeros(40)
heighthist1 = np.zeros(40)
heighthist2 = np.zeros(40)
heighthist3 = np.zeros(40)
heighthist4 = np.zeros(40)
heighthistall = np.zeros(40)
MISRheighthistall = np.zeros(40)'''

for desc in description:
    if desc == 'good':
        file_list = glob.glob(height_file_dir+'/'+desc+'/AST_L1T**.nc')
        #file_list = glob.glob(height_file_dir+'/'+desc+'/AST_L1T_00309**.nc')
        #file_list.extend(glob.glob(height_file_dir+'/'+desc+'/AST_L1T_00310**.nc'))
        for file in file_list:
            print(i)
            print(file)
            height_nc_file = nc.Dataset(file)

            h0, h1, h2, h3, h4, hall, d = get_counts(height_nc_file, 3)
            heighthist0 += h0
            heighthist1 += h1
            heighthist2 += h2
            heighthist3 += h3
            heighthist4 += h4
            heighthistall += hall
            diameterhist += d

            #MISRhall = getMISR_CTH_hist(file)
            #MISRheighthistall += MISRhall

            #print(MISRheighthistall)

            i += 1

    else:
        for add_des in sun_add:
            file_list = glob.glob(height_file_dir+'/'+desc+'/'+add_des+'/AST_L1T**.nc')
            #file_list = glob.glob(height_file_dir+'/'+desc+'/'+add_des+'/AST_L1T_00309**.nc')
            #file_list.extend(glob.glob(height_file_dir+'/'+desc+'/'+add_des+'/AST_L1T_00310**.nc'))
            for file in file_list:
                print(i)
                print(file)
                height_nc_file = nc.Dataset(file)

                h0, h1, h2, h3, h4, hall, d = get_counts(height_nc_file, 3)
                heighthist0 += h0
                heighthist1 += h1
                heighthist2 += h2
                heighthist3 += h3
                heighthist4 += h4
                heighthistall += hall
                diameterhist += d

                #MISRhall = getMISR_CTH_hist(file)
                #MISRheighthistall += MISRhall

                i += 1

#plt.rcParams.update({'font.size': 22})
fig, ax = plt.subplots()
all = np.nansum(heighthistall)
bincenter = np.arange(50,10050,100)
#bincenter = np.arange(125,10125,250)
plt.plot(heighthist0/all, bincenter, label='0-0.5', linewidth = 0.7)
plt.plot(heighthist1/all, bincenter, label='0-1.0', linewidth = 0.7)
plt.plot(heighthist2/all, bincenter, label='0-2.0', linewidth = 0.7)
plt.plot(heighthist3/all, bincenter, label='0-3.0', linewidth = 0.7)
plt.plot(heighthist4/all, bincenter, label='0-4.0', linewidth = 0.7)
plt.plot(heighthistall/all, bincenter, label='all', color='black')

plt.ylim(0,6000)
plt.xlim(0.001,0.2)
plt.xscale('log')
plt.xlabel('Normalized Frequency')
plt.ylabel('Cloud Height (m)')
#lgd=plt.legend(title='Cloud Diameter (km)',bbox_to_anchor=(1.02, 1), loc='upper left')
plt.legend(title='Cloud Diameter (km)')
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_cth1_dist_test_interp_min_log_6000_p1K.png', bbox_inches='tight', dpi=300)
#bbox_extra_artists=(lgd,), 
plt.clf()

heighthistall_array = np.array(heighthistall)
np.savetxt('/data/keeling/a/mdevera2/macrost/nolarge_cth1_heighthistall_interp_min_p1K.txt', heighthistall_array, delimiter =", ")

'''fig, ax = plt.subplots()
all = np.nansum(heighthistall)
MISRall = np.nansum(MISRheighthistall)
#bincenter = np.arange(50,10050,100)
bincenter = np.arange(125,10125,250)

plt.plot(heighthistall/all, bincenter, label='all')
plt.plot(MISRheighthistall/MISRall, bincenter, label='MISR all')

plt.xlim(0.001,0.2)
plt.xscale('log')
plt.xlabel('Normalized Frequency')
plt.ylabel('Cloud Height (m)')
plt.legend(title='Cloud Diameter (km)')
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_cth1_misr_dist_test_log.png', dpi=300)
plt.clf()

MISRheighthistall_array = np.array(MISRheighthistall)
np.savetxt('/data/keeling/a/mdevera2/macrost/nolarge_cth1_misr_heighthistall.txt', MISRheighthistall_array, delimiter =", ")'''


'''fig, ax = plt.subplots()
binrange = np.arange(0,7,0.1)
print(diameterhist)
diameterhist = np.insert(diameterhist, 0, diameterhist[0]/np.sum(diameterhist))
plt.step(binrange, diameterhist/np.sum(diameterhist))

plt.grid()

#plt.ylim(0,1)
#plt.xlim(0,10)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Normalized Frequency')
plt.xlabel('Cloud Diameter (km)')
plt.savefig('/data/keeling/a/mdevera2/macrost/cloud_size_90_test.png', dpi=300)'''