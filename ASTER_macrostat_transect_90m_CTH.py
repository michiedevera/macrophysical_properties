import netCDF4 as nc
import numpy as np
import collections
import pandas as pd
import os, glob, sys, getopt, argparse, re
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit, fsolve
import seaborn as sns

#folder = '/AGU_2022_poster/'
#plt.rcParams.update({'font.size': 22})
folder='/'

def power_law(x, a, b):
    return a*np.power(x, -b)

def power_law2(x, a):
    return a*np.power(x, -root)

def line_law(x, a, b):
    return (a - b*(x))

def direct(l):
    return (((1.0-l)*(np.power(du,2-l) - np.power(d0,2-l)))/((2.0-l)*(np.power(du,1-l) - np.power(d0,1-l)))-dmean)


in_dir = '/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/transects_'

#df = pd.read_csv(in_dir+'/good/AST_L1T_00308122019021254_20190813144107_30191_cf.txt')
#print(df)

description = ['good', 'sun_glint/okay']
#description = ['good']

cloudy_li=[]
clear_li=[]
trans_length = []

cloudy_nofil = []
clear_nofil = []
trans_length_nofil = []

direction = 'col'
for desc in description:
    print(desc)
    cloudy_file_list = glob.glob(in_dir+direction+'_90m/'+desc+'/cloudy_AST_L1T**.csv')
    for cloudy_file in cloudy_file_list:
        cloudy_df = pd.read_csv(cloudy_file, index_col=None, header=0)
        cloudy_nofil.append(cloudy_df)

        AST_file = cloudy_file.split('/')[-1].split('cloudy_')[-1]
        clear_file = in_dir+direction+'_90m/'+desc+'/clear_'+AST_file

        clear_df = pd.read_csv(clear_file, index_col=None, header=0)
        clear_nofil.append(clear_df)

        merged = pd.concat([cloudy_df, clear_df])
        sum_df = merged.groupby([direction])['count'].sum()
        row_col = sum_df[sum_df> 634].index

        trans_length_nofil.append(sum_df)
        trans_length.append(sum_df[sum_df> 634])
        cloudy_li.append(cloudy_df[cloudy_df[direction].isin(row_col)])
        clear_li.append(clear_df[clear_df[direction].isin(row_col)])

all_cloudy_nofil = pd.concat(cloudy_nofil, axis=0, ignore_index=True)
all_cloudy_nofil['transect_length'] = all_cloudy_nofil['count'] * 90e-3
all_cloudy = pd.concat(cloudy_li, axis=0, ignore_index=True)
all_cloudy['transect_length'] = all_cloudy['count'] * 90e-3
#all_cloudy['equiv_diameter'] = np.sqrt(4 * all_cloudy['count'] * 15 * 15 / np.pi)/1000
print(all_cloudy)
print(all_cloudy['transect_length'].max())

'''all_trans_length_nofil = pd.concat(trans_length_nofil, axis=0, ignore_index=True) * 90e-3
all_trans_length = pd.concat(trans_length, axis=0, ignore_index=True) * 90e-3
print(all_trans_length)

bins_list = list(np.arange(0,65.5,0.5))
n, bins, patches = plt.hist(x=all_trans_length_nofil, bins=bins_list)
plt.yscale('log')
plt.xlabel('Length (km)')
plt.ylabel('Counts')
plt.savefig('/data/keeling/a/mdevera2/macrost/nofil_nolarge_all_transect_'+direction+'_length_test_90m.png', bbox_inches='tight')
plt.close()

bins_list = list(np.arange(57,65.5,0.5))
n, bins, patches = plt.hist(x=all_trans_length, bins=bins_list)
plt.yscale('log')
plt.xlabel('Length (km)')
plt.ylabel('Counts')
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_all_transect_'+direction+'_length_test_90m.png', bbox_inches='tight')
plt.close()'''

#bins_list = list(np.arange(0,10.3,0.1))
#bin_center = list(np.arange(0.05,10.25,0.1))

bins_list = list(np.arange(0,10.1,0.1))
bin_center = list(np.arange(0.05,10.05,0.1))

#no filter
#cloudy_less7_nofil = all_cloudy_nofil[all_cloudy_nofil['transect_length'] < 7].copy()

cloudy_less7_nofil = all_cloudy_nofil.copy()
print(len(cloudy_less7_nofil['height1']))
weights = np.ones_like(cloudy_less7_nofil['height1'])/float(len(cloudy_less7_nofil['height1']))
n, bins, patches = plt.hist(x=cloudy_less7_nofil['height1']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='1')
print(n.sum())
n, bins, patches = plt.hist(x=cloudy_less7_nofil['height2']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='2')
print(n.sum())
n, bins, patches = plt.hist(x=cloudy_less7_nofil['height3']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='3')
print(n.sum())

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nofil_nolarge_cloudy_transect_'+direction+'_cth_dist_test_90m.png', bbox_inches='tight')
plt.close()


#filtered
#cloudy_less7 = all_cloudy[all_cloudy['transect_length'] < 7].copy()
cloudy_less7 = all_cloudy.copy()
print(len(cloudy_less7['height1']))
weights = np.ones_like(cloudy_less7['height1'])/float(len(cloudy_less7['height1']))
n, bins, patches = plt.hist(x=cloudy_less7['height1']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='ASTER')
#n, bins, patches = plt.hist(x=cloudy_less7['height2']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='2')
#n, bins, patches = plt.hist(x=cloudy_less7['height3']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='3')

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_'+direction+'_cth_dist_test_90m.png', bbox_inches='tight', dpi=300)
plt.close()

hsrl = pd.read_csv('/data/keeling/a/mdevera2/macrost/17_RF_HSRL_2Hz_Cloud_element.csv', index_col=None, header=0)
rsp = pd.read_csv('/data/keeling/a/mdevera2/macrost/17RF_RSP_cloud_element.csv', index_col=None, header=0)

hsrl_cloud = hsrl[hsrl['clear0/cloudy1'] == 1].copy()
hsrl_cloud_less7 = hsrl_cloud[hsrl_cloud['horizontal_length(km)'] < 7]['Mean_CTH_HSRL_2Hz']*1e-3
hsrl_cloud_less7_less4km = hsrl_cloud[(hsrl_cloud['horizontal_length(km)'] < 7) & (hsrl_cloud['Mean_CTH_HSRL_2Hz'] < 4000)]['Mean_CTH_HSRL_2Hz']*1e-3
print(hsrl_cloud_less7.max())

hsrl_cloud_nocir = hsrl[(hsrl['clear0/cloudy1'] == 1) & (hsrl['Mean_SPNs_Trans'] > 0.95)].copy()
hsrl_cloud_nocir_less7 = hsrl_cloud_nocir[hsrl_cloud_nocir['horizontal_length(km)'] < 7]['Mean_CTH_HSRL_2Hz']*1e-3
hsrl_cloud_nocir_less7_less4km = hsrl_cloud_nocir[(hsrl_cloud_nocir['horizontal_length(km)'] < 7) & (hsrl_cloud_nocir['Mean_CTH_HSRL_2Hz'] < 4000)]['Mean_CTH_HSRL_2Hz']*1e-3

rsp_cloud = rsp[rsp['clear0/cloudy1'] == 1].copy()
rsp_cloud_less7 = rsp_cloud[rsp_cloud['horizontal_length(km)'] < 7]['Mean_CTH_stereo']
rsp_cloud_less7_less4km = rsp_cloud[(rsp_cloud['horizontal_length(km)'] < 7) & (rsp_cloud['Mean_CTH_stereo'] < 4)]['Mean_CTH_stereo']
print(rsp_cloud_less7)

rsp_cloud_nocir = rsp[(rsp['clear0/cloudy1'] == 1) & (rsp['Mean_SPNs_Trans'] > 0.95)].copy()
rsp_cloud_nocir_less7 = rsp_cloud_nocir[rsp_cloud_nocir['horizontal_length(km)'] < 7]['Mean_CTH_stereo']
rsp_cloud_nocir_less7_less4km = rsp_cloud_nocir[(rsp_cloud_nocir['horizontal_length(km)'] < 7) & (rsp_cloud_nocir['Mean_CTH_stereo'] < 4)]['Mean_CTH_stereo']

weights = np.ones_like(hsrl_cloud_less7)/float(len(hsrl_cloud_less7))
nhsrl, bins, patches = plt.hist(x=hsrl_cloud_less7, bins=bins_list, histtype=u'step', weights=weights, color='blue', label='HSRL-2')

plt.yscale('log')

cumulativehsrl = np.cumsum(nhsrl)
cumulativehsrl = np.insert(cumulativehsrl,0,cumulativehsrl[0])

plt.step(bins, cumulativehsrl, linestyle='--', color='blue')

plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_hsrl_cth_dist_test_90m.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(hsrl_cloud_less7_less4km)/float(len(hsrl_cloud_less7_less4km))
nhsrl, bins, patches = plt.hist(x=hsrl_cloud_less7_less4km, bins=bins_list, histtype=u'step', weights=weights, color='blue', label='HSRL-2')

plt.yscale('log')

cumulativehsrl = np.cumsum(nhsrl)
cumulativehsrl = np.insert(cumulativehsrl,0,cumulativehsrl[0])

plt.step(bins, cumulativehsrl, linestyle='--', color='blue')

plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_hsrl_cth_dist_test_90m_4km.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(hsrl_cloud_nocir_less7)/float(len(hsrl_cloud_nocir_less7))
nhsrl, bins, patches = plt.hist(x=hsrl_cloud_nocir_less7, bins=bins_list, histtype=u'step', weights=weights, color='blue', label='HSRL-2')

plt.yscale('log')

cumulativehsrl = np.cumsum(nhsrl)
cumulativehsrl = np.insert(cumulativehsrl,0,cumulativehsrl[0])

plt.step(bins, cumulativehsrl, linestyle='--', color='blue')

plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_hsrl_nocir_cth_dist_test_90m.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(hsrl_cloud_nocir_less7_less4km)/float(len(hsrl_cloud_nocir_less7_less4km))
nhsrl, bins, patches = plt.hist(x=hsrl_cloud_nocir_less7_less4km, bins=bins_list, histtype=u'step', weights=weights, color='blue', label='HSRL-2')

plt.yscale('log')

cumulativehsrl = np.cumsum(nhsrl)
cumulativehsrl = np.insert(cumulativehsrl,0,cumulativehsrl[0])

plt.step(bins, cumulativehsrl, linestyle='--', color='blue')

plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_hsrl_nocir_cth_dist_test_90m_4km.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(rsp_cloud_less7)/float(len(rsp_cloud_less7))
nrsp, bins, patches = plt.hist(x=rsp_cloud_less7, bins=bins_list, histtype=u'step', weights=weights, color='red', label='RSP')

cumulativersp = np.cumsum(nrsp)
cumulativersp =np.insert(cumulativersp,0,cumulativersp[0])

plt.step(bins, cumulativersp, linestyle='--', color='red')

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_rsp_cth_dist_test_90m.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(rsp_cloud_less7_less4km)/float(len(rsp_cloud_less7_less4km))
nrsp, bins, patches = plt.hist(x=rsp_cloud_less7_less4km, bins=bins_list, histtype=u'step', weights=weights, color='red', label='RSP')

cumulativersp = np.cumsum(nrsp)
cumulativersp =np.insert(cumulativersp,0,cumulativersp[0])

plt.step(bins, cumulativersp, linestyle='--', color='red')

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_rsp_cth_dist_test_90m_4km.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(rsp_cloud_nocir_less7)/float(len(rsp_cloud_nocir_less7))
nrsp, bins, patches = plt.hist(x=rsp_cloud_nocir_less7, bins=bins_list, histtype=u'step', weights=weights, color='red', label='RSP')

cumulativersp = np.cumsum(nrsp)
cumulativersp =np.insert(cumulativersp,0,cumulativersp[0])

plt.step(bins, cumulativersp, linestyle='--', color='red')

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_rsp_nocir_cth_dist_test_90m.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(rsp_cloud_nocir_less7_less4km)/float(len(rsp_cloud_nocir_less7_less4km))
nrsp, bins, patches = plt.hist(x=rsp_cloud_nocir_less7_less4km, bins=bins_list, histtype=u'step', weights=weights, color='red', label='RSP')

cumulativersp = np.cumsum(nrsp)
cumulativersp =np.insert(cumulativersp,0,cumulativersp[0])

plt.step(bins, cumulativersp, linestyle='--', color='red')

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_rsp_nocir_cth_dist_test_90m_4km.png', bbox_inches='tight', dpi=300)
plt.close()

print(n.sum())
weights = np.ones_like(cloudy_less7['height1'])/float(len(cloudy_less7['height1']))
n, bins, patches = plt.hist(x=cloudy_less7['height1']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='ASTER')
#n, bins, patches = plt.hist(x=cloudy_less7['height2']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='2')
#n, bins, patches = plt.hist(x=cloudy_less7['height3']*1e-3, bins=bins_list, histtype=u'step', weights=weights, label='3')
weights = np.ones_like(hsrl_cloud_less7)/float(len(hsrl_cloud_less7))
nhsrl, bins, patches = plt.hist(x=hsrl_cloud_less7, bins=bins_list, histtype=u'step', weights=weights, color='blue', label='HSRL-2')
weights = np.ones_like(rsp_cloud_less7)/float(len(rsp_cloud_less7))
nrsp, bins, patches = plt.hist(x=rsp_cloud_less7, bins=bins_list, histtype=u'step', weights=weights, color='red', label='RSP')

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_'+direction+'_hsrl_rsp_cth_dist_test_90m.png', bbox_inches='tight', dpi=300)
plt.close()

print(n.sum())
print(nhsrl.sum())
print(nrsp.sum())

weights = np.ones_like(cloudy_less7['height1'])/float(len(cloudy_less7['height1']))
n, _ = np.histogram(cloudy_less7['height1']*1e-3, bins=bins_list, weights=weights)
weights = np.ones_like(hsrl_cloud_less7)/float(len(hsrl_cloud_less7))
nhsrl, _ = np.histogram(hsrl_cloud_less7, bins=bins_list, weights=weights)
weights = np.ones_like(rsp_cloud_less7)/float(len(rsp_cloud_less7))
nrsp, _ = np.histogram(rsp_cloud_less7, bins=bins_list, weights=weights)

plt.plot(n, bin_center, label='ASTER', linewidth = 0.7)
plt.plot(nhsrl, bin_center, label='HSRL-2', color='blue', linewidth = 0.7)
plt.plot(nrsp, bin_center, label='RSP', color='red', linewidth = 0.7)

plt.xscale('log')
plt.legend()
plt.grid()
plt.ylabel('Height (km)')
plt.xlabel('Normalized Frequency')
#plt.ylim(0,1)
plt.ylim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_'+direction+'_hsrl_rsp_cth_dist_test_90m_vertical.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(cloudy_less7['height1'])/float(len(cloudy_less7['height1']))
n, _ = np.histogram(cloudy_less7['height1']*1e-3, bins=bins_list, weights=weights)
weights = np.ones_like(hsrl_cloud_less7_less4km)/float(len(hsrl_cloud_less7_less4km))
nhsrl, _ = np.histogram(hsrl_cloud_less7_less4km, bins=bins_list, weights=weights)
weights = np.ones_like(rsp_cloud_less7_less4km)/float(len(rsp_cloud_less7_less4km))
nrsp, _ = np.histogram(rsp_cloud_less7_less4km, bins=bins_list, weights=weights)

plt.plot(n, bin_center, label='ASTER', linewidth = 0.7)
plt.plot(nhsrl, bin_center, label='HSRL-2', color='blue', linewidth = 0.7)
plt.plot(nrsp, bin_center, label='RSP', color='red', linewidth = 0.7)

plt.xscale('log')
plt.legend()
plt.grid()
plt.ylabel('Height (km)')
plt.xlabel('Normalized Frequency')
#plt.ylim(0,1)
plt.ylim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_'+direction+'_hsrl_rsp_cth_dist_test_90m_4km_vertical.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(cloudy_less7['height1'])/float(len(cloudy_less7['height1']))
n, _ = np.histogram(cloudy_less7['height1']*1e-3, bins=bins_list, weights=weights)
weights = np.ones_like(hsrl_cloud_nocir_less7)/float(len(hsrl_cloud_nocir_less7))
nhsrl, _ = np.histogram(hsrl_cloud_nocir_less7, bins=bins_list, weights=weights)
weights = np.ones_like(rsp_cloud_nocir_less7)/float(len(rsp_cloud_nocir_less7))
nrsp, _ = np.histogram(rsp_cloud_nocir_less7, bins=bins_list, weights=weights)

plt.plot(n, bin_center, label='ASTER', linewidth = 0.7)
plt.plot(nhsrl, bin_center, label='HSRL-2', color='blue', linewidth = 0.7)
plt.plot(nrsp, bin_center, label='RSP', color='red', linewidth = 0.7)

plt.xscale('log')
plt.legend()
plt.grid()
plt.ylabel('Height (km)')
plt.xlabel('Normalized Frequency')
#plt.ylim(0,1)
plt.ylim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_'+direction+'_hsrl_rsp_nocir_cth_dist_test_90m_vertical.png', bbox_inches='tight', dpi=300)
plt.close()

weights = np.ones_like(cloudy_less7['height1'])/float(len(cloudy_less7['height1']))
n, _ = np.histogram(cloudy_less7['height1']*1e-3, bins=bins_list, weights=weights)
weights = np.ones_like(hsrl_cloud_nocir_less7_less4km)/float(len(hsrl_cloud_nocir_less7_less4km))
nhsrl, _ = np.histogram(hsrl_cloud_nocir_less7_less4km, bins=bins_list, weights=weights)
weights = np.ones_like(rsp_cloud_nocir_less7_less4km)/float(len(rsp_cloud_nocir_less7_less4km))
nrsp, _ = np.histogram(rsp_cloud_nocir_less7_less4km, bins=bins_list, weights=weights)

plt.plot(n, bin_center, label='ASTER', linewidth = 0.7)
plt.plot(nhsrl, bin_center, label='HSRL-2', color='blue', linewidth = 0.7)
plt.plot(nrsp, bin_center, label='RSP', color='red', linewidth = 0.7)

plt.xscale('log')
plt.legend()
plt.grid()
plt.ylabel('Height (km)')
plt.xlabel('Normalized Frequency')
#plt.ylim(0,1)
plt.ylim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_'+direction+'_hsrl_rsp_nocir_cth_dist_test_90m_4km_vertical.png', bbox_inches='tight', dpi=300)
plt.close()

cth1 = np.loadtxt('/data/keeling/a/mdevera2/macrost/nolarge_cth1_heighthistall_interp_min.txt')

fig, ax = plt.subplots()

plt.step(bins_list[:-1], cth1/np.nansum(cth1), label='ASTER', where='post')

hsrl_cloud_full = hsrl_cloud['Mean_CTH_HSRL_2Hz']*1e-3
rsp_cloud_full = rsp_cloud['Mean_CTH_stereo']

weights = np.ones_like(hsrl_cloud_full)/float(len(hsrl_cloud_full))
nhsrl, bins, patches = plt.hist(x=hsrl_cloud_full, bins=bins_list, histtype=u'step', weights=weights, color='blue', label='HSRL-2')
weights = np.ones_like(rsp_cloud_full)/float(len(rsp_cloud_full))
nrsp, bins, patches = plt.hist(x=rsp_cloud_full, bins=bins_list, histtype=u'step', weights=weights, color='red', label='RSP')

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_cth1_hsrl_rsp_cth_dist_test_90m.png', bbox_inches='tight', dpi=300)
plt.close()


weights = np.ones_like(hsrl_cloud_full)/float(len(hsrl_cloud_full))
nhsrl, _ = np.histogram(hsrl_cloud_full, bins=bins_list, weights=weights)
weights = np.ones_like(rsp_cloud_full)/float(len(rsp_cloud_full))
nrsp, _ = np.histogram(rsp_cloud_full, bins=bins_list, weights=weights)

plt.plot(cth1/np.nansum(cth1), bin_center, label='ASTER', linewidth = 0.7)
plt.plot(nhsrl, bin_center, label='HSRL-2', color='blue', linewidth = 0.7)
plt.plot(nrsp, bin_center, label='RSP', color='red', linewidth = 0.7)

plt.xscale('log')
plt.legend()
plt.grid()
plt.ylabel('Height (km)')
plt.xlabel('Normalized Frequency')
#plt.ylim(0,1)
plt.ylim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_cth1_hsrl_rsp_cth_dist_test_90m_vertical.png', bbox_inches='tight', dpi=300)
plt.close()


bins_list = list(np.arange(0,10,0.15))
hsrl_cloud_full = hsrl_cloud['Mean_CTH_HSRL_2Hz']*1e-3
rsp_cloud_full = rsp_cloud['Mean_CTH_stereo']

weights = np.ones_like(hsrl_cloud_full)/float(len(hsrl_cloud_full))
nhsrl, bins, patches = plt.hist(x=hsrl_cloud_full, bins=bins_list, histtype=u'step', weights=weights, color='blue', label='HSRL-2')
weights = np.ones_like(rsp_cloud_full)/float(len(rsp_cloud_full))
nrsp, bins, patches = plt.hist(x=rsp_cloud_full, bins=bins_list, histtype=u'step', weights=weights, color='red', label='RSP')

plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Height (km)')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloudy_transect_Fu_hsrl_rsp_cth_dist_test_90m.png', bbox_inches='tight', dpi=300)
plt.close()




'''
#bins_list = list(np.arange(0,10.3,0.1))
weights = np.ones_like(clear_less7['transect_length'])/float(len(clear_less7['transect_length']))
n, bins, patches = plt.hist(x=clear_less7['transect_length'], bins=bins_list, histtype=u'step', weights=weights, label='Observations')
plt.xscale('log')
plt.yscale('log')

plt.grid()
plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_clear_transect_row_test_90m.png', bbox_inches='tight')
plt.close()
'''
'''with open('/data/keeling/a/mdevera2/macrost/rsp_chl.out') as f:
    rsp = f.readlines()
with open('/data/keeling/a/mdevera2/macrost/hsrl_chl.out') as f:
    hsrl = f.readlines()

rsp_trans_length = [float(i.strip()) for i in rsp]
hsrl_trans_length = [float(i.strip()) for i in hsrl]

bins_list = list(np.arange(0,10.3,0.1))
weights = np.ones_like(rsp_trans_length)/float(len(rsp_trans_length))
nrsp, bins, patches = plt.hist(x=rsp_trans_length, bins=bins_list, histtype=u'step', weights=weights, color='red')

cumulativersp = np.cumsum(nrsp)
cumulativersp =np.insert(cumulativersp,0,cumulativersp[0])

plt.step(bins, cumulativersp, linestyle='--', color='red')

plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_cloudy_transect_rsp_test.png', bbox_inches='tight')
plt.close()'''



'''all['cloud_area_km2'] = all['numpixels'] * 15 * 15e-6
all['equiv_diameter'] = np.sqrt(4 * all['numpixels'] * 15 * 15 / np.pi)/1000
all['cloud_perimeter_km'] = all['numedge'] * 15e-3
print(all)
#2181059

#all= all[all['numpixels'] > 1].copy().reset_index(drop=True)
#glit = all[all['numpixels'] == 1].copy()
#print(glit)
#878635

#all.to_csv('/data/keeling/a/mdevera2/macrost/all_macrost_nolarge.csv', index=False)

less7 = all[all['equiv_diameter'] < 7]
dmean=less7['equiv_diameter'].mean()
d0=less7['equiv_diameter'].min()
du=less7['equiv_diameter'].max()

bins_list = list(np.arange(0,7,0.1))
weights = np.ones_like(less7['equiv_diameter'])/float(len(less7['equiv_diameter']))
n, bins, patches = plt.hist(x=less7['equiv_diameter'], bins=bins_list, histtype=u'step', weights=weights, label='Observations')
plt.xscale('log')
plt.yscale('log')

bin_center = list(np.arange(0.05,6.95,0.1))

n[n==0]=1

#mu, sigma = stats.norm.fit(all['equiv_diameter'])
#best_fit_line = stats.norm.pdf(bins, mu, sigma)
#plt.plot(bins, best_fit_line)

#pars, cov = curve_fit(f=power_law, xdata=bin_center, ydata=n)
#print(pars)

#plt.plot(bin_center, power_law(bin_center, *pars), linestyle='--', linewidth=1.5, color='black', label=str(pars[1]))

pars2, cov2 = curve_fit(f=line_law, xdata=np.log(bin_center), ydata=np.log(n))
#print(np.log(n))
print(pars2)
print(cov2)
#print(line_law(np.array(np.log(bin_center)), *pars2))
print('line ' + str(pars2[1]))

plt.plot(bin_center, np.exp(line_law(np.array(np.log(bin_center)), *pars2)), linestyle='--', linewidth=1.5, color='black', label='Line Fit')

#m,b = np.polyfit(np.log(bin_center), np.log(n), 1)
#plt.plot(bin_center, np.exp(m*np.log(bin_center)+b), linestyle='--', linewidth=1.5, color='purple', label=str(m))
#print(b)

root = fsolve(direct, 2.9)
print(root)
print(n)
print('direct ' + str(root[0]))

pars3, cov = curve_fit(f=power_law2, xdata=bin_center, ydata=n)
print(pars3)
print(cov)

print(power_law(bin_center, pars3[0], root[0]))
y = power_law(bin_center, pars3[0], root[0])
y = np.insert(y, 0, power_law(0.05, pars3[0], root[0]))

plt.step(bins, y, linestyle='--', linewidth=1.5, color='green', label='Direct Power-Law Fit')

plt.grid()
plt.ylim(0,1)
plt.xlim(0,10)
plt.xlabel('Cloud Equivalent Diameter (km)')
plt.ylabel('Normalized Frequency')
lgd=plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

#plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_cloud_size_test.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

residuals1 = np.log(n) - line_law(np.log(bin_center), *pars2)
ss_res1 = np.sum(residuals1**2)
ss_tot1 = np.sum((np.log(n)-np.mean(np.log(n)))**2)
r_squared1 = 1 - (ss_res1 / ss_tot1)
print(r_squared1)
print(np.sqrt(r_squared1))

residuals2 = n - power_law(bin_center, pars3[0], root[0])
ss_res2 = np.sum(residuals2**2)
ss_tot2 = np.sum((n-np.mean(n))**2)
r_squared2 = 1 - (ss_res2 / ss_tot2)
print(r_squared2)
print(np.sqrt(r_squared2))

print(all[all['equiv_diameter'] > 30])
print(all['equiv_diameter'].max())
total_cf=cf['cloudy'].sum()/cf['total'].sum()
print(cf['cloudy'].sum()/cf['total'].sum())
print(all['numpixels'].sum()/cf['total'].sum())
#print((cf['cloudy'].sum()/cf['total'].sum())-(cf['cloudy'][all['equiv_diameter'].argmax()]/cf['total'][all['equiv_diameter'].argmax()]))
big = all['numpixels'][all['equiv_diameter'].argmax()]
print((cf['cloudy'].sum()-big)/(cf['total'].sum()-big))
print(all['numpixels'][all['equiv_diameter'].argmax()])
#glitter=glit['numpixels'].sum()
#print((cf['cloudy'].sum()-glitter)/(cf['total'].sum()))
bins_list = list(np.arange(0,30,0.1))
bin_center = list(np.arange(0.05,29.95,0.1))

sums, edges = np.histogram(all['equiv_diameter'], bins=bins_list, weights=all['numpixels'])
counts, _ = np.histogram(all['equiv_diameter'], bins=bins_list)
cf_dist = sums / cf['total'].sum()
print(cf_dist)
print(np.cumsum(cf_dist))


plt.plot(bin_center, cf_dist, label='Cloud Fraction')
plt.plot(bin_center, np.cumsum(cf_dist), linestyle='--', linewidth=2, color='black', label='Cumulative Cloud Fraction')
plt.yscale('log')
plt.xlabel('Cloud Equivalent Diameter (km)')
plt.ylabel('Fraction')
#plt.xticks(np.arange(0, 31, 1))
#plt.axhline(y=total_cf/2, color='r', linestyle='--', linewidth=1)
plt.grid()
lgd = plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
#plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_cloud_frac_dist_test.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()


greater12=all[all['numpixels'] > 12]
print(greater12['numedge'].min())
print(greater12['cloud_perimeter_km'].min())
plt.scatter(greater12['cloud_perimeter_km'], greater12['cloud_area_km2'], s=2)

pars4, cov4 = curve_fit(f=line_law, xdata=np.log(greater12['cloud_perimeter_km']), ydata=np.log(greater12['cloud_area_km2']))
print(cov4)
print('d=' + str(-pars4[1]))

residuals = np.log(greater12['cloud_area_km2'])- line_law(np.log(greater12['cloud_perimeter_km']), *pars4)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((np.log(greater12['cloud_area_km2'])-np.mean(np.log(greater12['cloud_area_km2'])))**2)
r_squared = 1 - (ss_res / ss_tot)
print(r_squared)
print(np.sqrt(r_squared))

peri = np.arange(0.2,1e4,1)
plt.plot(peri, np.exp(line_law(np.array(np.log(peri)), *pars4)), linestyle='--', linewidth=1.5, color='black', label='d='+str(-pars4[1])[:4])

plt.xscale('log')
plt.yscale('log')
plt.grid()

plt.ylabel('Cloud Area (km$^2$)')
plt.xlabel('Perimeter (km)')
plt.legend()
plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_cloud_peri_area_test.png')
plt.close()

bins_list = list(np.arange(0,1050,50))
weights = np.ones_like(all['nndist']*15)/float(len(all['nndist']))
n, bins, patches = plt.hist(x=all['nndist']*15, bins=bins_list, histtype=u'step', weights=weights)
plt.grid()
plt.xlabel('NN Distance (m)')
plt.ylabel('Normalized Frequency')
plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_cloud_nn_dist_test.png')
plt.close()

all['nndist/radius'] = all['nndist']*15/((all['equiv_diameter']*1e3)/2.)
bins_list = list(np.arange(0,110,10))
weights = np.ones_like(all['nndist/radius'])/float(len(all['nndist/radius']))
n, bins, patches = plt.hist(x=all['nndist/radius'], bins=bins_list, histtype=u'step', weights=weights)
plt.grid()
plt.xlabel('NN Distance / Radius')
plt.ylabel('Normalized Frequency')
plt.ylim(0,1)
plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/nolarge_cloud_nn_dist_radius_test.png')'''