#from osgeo import gdal, osr
import numpy as np
import os, glob, sys, getopt, argparse
import MisrToolkit as Mtk
#import pandas as pd
#from pyresample import create_area_def
#import ASTER_cloud_mask_funcs as Acm
import matplotlib.pyplot as plt

orbits=[104351,104365,104366,104380,104394,104395,104409,104424,104438,104439,104453,104467,104468,104482,104496,104497,104511,104526,104540,104541,104555,104569,104570,104584,104598,104599,104613,104627,104628,104642,104657,104671,104672,104686,104700,104701,104715,104729,104730,104744,104759,104773,104774,104788]
MISRall = np.zeros(100)

for orbit in orbits:
    print(orbit)
    r = Mtk.MtkRegion(25,110,0,135)
    yy = Mtk.orbit_to_time_range(orbit)[0][0:4]
    mm = Mtk.orbit_to_time_range(orbit)[0][5:7]
    dd = Mtk.orbit_to_time_range(orbit)[0][8:10]
    date = yy + '.' + mm + '.' + dd
    print(date)

    try:
        tccloudfile = glob.glob('/data/gdi/e/MISR/TCCloud/'+date+'/*'+str(orbit)+'*.hdf')[0]
    except:
        tccloudfile = glob.glob('/data/gdi/c/mdevera2/MISR_TCCloud/*'+str(orbit)+'*.hdf')[0]
    #    continue
    m = Mtk.MtkFile(tccloudfile)
    hbig = m.grid('Stereo_1.1_km').field('CloudTopHeight').read(r).data()

    binrange = np.arange(0,10100,100)
    histall, _ = np.histogram(hbig[hbig>0], bins=binrange)

    MISRall += histall

fig, ax = plt.subplots()
total = np.sum(MISRall)
bincenter = np.arange(50,10050,100)

plt.plot(MISRall/total, bincenter, label='MISR whole region')

plt.xlim(0.001,0.2)
plt.xscale('log')
plt.xlabel('Normalized Frequency')
plt.ylabel('Cloud Height (m)')
#plt.legend(title='Cloud Diameter (km)')
plt.savefig('/data/keeling/a/mdevera2/macrost/cth_misr_big_dist_test_log.png', dpi=300)
plt.clf()

MISRall_array = np.array(MISRall)
np.savetxt('/data/keeling/a/mdevera2/macrost/cth_misr_big_heighthistall.txt', MISRall_array, delimiter =", ")