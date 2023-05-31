import netCDF4 as nc
import numpy as np
import collections
import pandas as pd
import os, glob, sys, getopt, argparse, re
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit, fsolve
import mpl_scatter_density
import seaborn as sns
import matplotlib.colors

#plt.rcParams.update({'font.size': 22})
#folder = '/AGU_2022_poster/'
folder = ''

def power_law(x, a, b):
    return a*np.power(x, -b)

def power_law2(x, a):
    return a*np.power(x, -root[0])

def line_law(x, a, b):
    return (a - b*(x))

def direct(l):
    return (((1.0-l)*(np.power(du,2-l) - np.power(d0,2-l)))/((2.0-l)*(np.power(du,1-l) - np.power(d0,1-l)))-dmean)


in_dir = '/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/cloudfraction'

#df = pd.read_csv(in_dir+'/good/AST_L1T_00308122019021254_20190813144107_30191_cf.txt')
#print(df)

description = ['good', 'sun_glint/okay']
#description = ['good']

cf = pd.DataFrame(columns = ['cloudy', 'total', 'cloud_fraction'])

for desc in description:
    file_list = glob.glob(in_dir+'/'+desc+'/AST_L1T**.txt')
    for file in file_list:
        with open(file) as f:
            contents = f.read().split()

        cf = cf.append({'cloudy' : int(contents[0]), 'total' : int(contents[1]), 'cloud_fraction' : float(contents[2])}, ignore_index = True)

bins_list = list(np.arange(0,0.7,0.05))
n, bins, patches = plt.hist(x=cf['cloud_fraction'], bins=bins_list)
plt.xlabel('Cloud Fraction')
plt.ylabel('Number of ASTER Scenes')
plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cf_test.png')
plt.close()

print(n.sum())
#print(cf['cloud_fraction'].max())
#max1 = cf['cloud_fraction'].max()
#cf.drop(cf.index[cf['cloud_fraction'] == max1], inplace=True)
#print(cf['cloud_fraction'].max())
#max1 = cf['cloud_fraction'].max()
#cf.drop(cf.index[cf['cloud_fraction'] == max1], inplace=True)
#print(cf['cloud_fraction'].max())


li=[]

in_dir = '/data/gdi/c/mdevera2/ASTER_macrostats/ocean_no_cirrus_scenes2/indiv_counts'
for desc in description:
    print(desc)
    file_list = glob.glob(in_dir+'/'+desc+'/AST_L1T**.csv')
    for file in file_list:
        df = pd.read_csv(file, index_col=None, header=0)
        li.append(df)

all = pd.concat(li, axis=0, ignore_index=True)
print(all)

all['cloud_area_km2'] = all['numpixels'] * 15 * 15e-6
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

np.savetxt('/data/keeling/a/mdevera2/macrost/ASTER_cloud_equiv_diameter_hist.txt', n)

bin_center = list(np.arange(0.05,6.95,0.1))
idx = np.where(n > 0)
    
pars, cov = curve_fit(f=power_law, xdata=np.array(bin_center)[idx], ydata=n[idx])
print(pars)

a,b = np.polyfit(np.log(np.array(bin_center)[idx]), np.log(n[idx]), 1)
print(a,b)

#print(n)
#n[n==0]=1
#print(n)
#print(np.log(n))
#mu, sigma = stats.norm.fit(all['equiv_diameter'])
#best_fit_line = stats.norm.pdf(bins, mu, sigma)
#plt.plot(bins, best_fit_line)

#plt.plot(bin_center, power_law(bin_center, *pars), linestyle='--', linewidth=1.5, color='black', label=str(pars[1]))

pars2, cov2 = curve_fit(f=line_law, xdata=np.log(np.array(bin_center)[idx]), ydata=np.log(n[idx]))
#print(np.log(n))
print(pars2)
print(cov2)
#print(line_law(np.array(np.log(bin_center)), *pars2))
print('line ' + str(pars2[1]))

plt.plot(bin_center, np.exp(line_law(np.array(np.log(bin_center)), *pars2)), linestyle='--', linewidth=1.5, color='black', label='Line Fit ($\lambda$=%.2f)' % pars2[1])

#m,b = np.polyfit(np.log(bin_center), np.log(n), 1)
#plt.plot(bin_center, np.exp(m*np.log(bin_center)+b), linestyle='--', linewidth=1.5, color='purple', label=str(m))
#print(b)

root = fsolve(direct, 2.9)
print(root)
print(n)
print('direct ' + str(root[0]))

pars3, cov = curve_fit(f=power_law2, xdata=np.array(bin_center)[idx], ydata=n[idx])
print(pars3)
print(cov)

print(power_law(bin_center, pars3[0], root[0]))
y = power_law(bin_center, pars3[0], root[0])
y = np.insert(y, 0, power_law(0.05, pars3[0], root[0]))

plt.step(bins, y, linestyle='--', linewidth=1.5, color='green', label='Direct Power-Law Fit ($\lambda$=%.2f)' % root[0])

plt.grid()
plt.ylim(0,1)
plt.xlim(0,10)
plt.xlabel('Cloud Equivalent Diameter (km)')
plt.ylabel('Normalized Frequency')
#lgd=plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
plt.legend()

#plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloud_size_test.png', bbox_inches='tight') #bbox_extra_artists=(lgd,),
plt.close()

print(np.sum(less7['equiv_diameter']*less7['equiv_diameter'])/np.sum(less7['equiv_diameter']))

n_counts, _ = np.histogram(less7['equiv_diameter'], bins=bins_list)
weights = np.ones_like(less7['equiv_diameter'])/float(len(less7['equiv_diameter']))
n_counts_d, _ = np.histogram(less7['equiv_diameter'], bins=bins_list, weights=weights)
A_counts, _ = np.histogram(less7['equiv_diameter'], bins=bins_list, weights=less7['cloud_area_km2'])

sigma = n_counts_d*(A_counts/n_counts)
print(np.sum(sigma*bin_center*0.1)/np.sum(sigma*0.1))
#1.6389856267898057

sigma2 = A_counts/(0.1*170*60*60)
print(np.sum(sigma2*bin_center*0.1)/np.sum(sigma2*0.1)) #1.6389856267917147
plt.stairs(sigma, bins_list)
plt.stairs(sigma2, bins_list)
#plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.xlabel('Length (km)')
plt.ylabel('Normalized Frequency')
#plt.ylim(0,1)
plt.xlim(0,10)
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/sigma_nolarge_cloudy_test.png', bbox_inches='tight')
plt.close()

'''dens, _, _ = plt.hist(x=less7['equiv_diameter'], bins=bins_list, histtype=u'step', label='Observations')
n = dens/(0.1*170*60*60)
np.savetxt('/data/keeling/a/mdevera2/macrost/ASTER_cloud_equiv_diameter_PDF.txt', n)'''

residuals1 = np.log(n) - line_law(np.log(bin_center), *pars2)
ss_res1 = np.sum(residuals1**2)
ss_tot1 = np.sum((np.log(n)-np.mean(np.log(n)))**2)
r_squared1 = 1 - (ss_res1 / ss_tot1)
print(r_squared1)
print('line r '+ str(np.sqrt(r_squared1)))

residuals2 = n - power_law(bin_center, pars3[0], root[0])
ss_res2 = np.sum(residuals2**2)
ss_tot2 = np.sum((n-np.mean(n))**2)
r_squared2 = 1 - (ss_res2 / ss_tot2)
print(r_squared2)
print('direct r '+ str(np.sqrt(r_squared2)))

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
#lgd = plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
plt.legend()
#plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloud_frac_dist_test.png', bbox_inches='tight') #bbox_extra_artists=(lgd,),
plt.close()


greater12=all[all['numpixels'] > 12]
print(greater12['numedge'].min())
print(greater12['cloud_perimeter_km'].min())
plt.scatter(greater12['cloud_area_km2'], greater12['cloud_perimeter_km'], s=0.5)

pars4, cov4 = curve_fit(f=line_law, xdata=np.log(greater12['cloud_area_km2']), ydata=np.log(greater12['cloud_perimeter_km']))
print(cov4)
print('d=' + str(-pars4[1]))

residuals = np.log(greater12['cloud_perimeter_km'])- line_law(np.log(greater12['cloud_area_km2']), *pars4)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((np.log(greater12['cloud_perimeter_km'])-np.mean(np.log(greater12['cloud_perimeter_km'])))**2)
r_squared = 1 - (ss_res / ss_tot)
print(r_squared)
print(np.sqrt(r_squared))

peri = np.arange(0.002,1e4,1)
plt.plot(peri, np.exp(line_law(np.array(np.log(peri)), *pars4)), linestyle='--', linewidth=1.5, color='black', label='d='+str(-2*pars4[1])[:4])

plt.xscale('log')
plt.yscale('log')
plt.grid()

plt.xlabel('Cloud Area (km$^2$)')
plt.ylabel('Perimeter (km)')
plt.legend()
plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloud_peri_area_test.png')
plt.close()

y = greater12['cloud_perimeter_km'].copy()
x = greater12['cloud_area_km2'].copy()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='scatter_density')

ybin = np.logspace(np.log10(0.1), np.log10(10000.0), 50)
xbin = np.logspace(np.log10(0.001), np.log10(10000.0), 50)
weights = np.ones_like(greater12['cloud_area_km2'])/float(len(greater12['cloud_area_km2']))

d=ax.hist2d(x, y, bins=[xbin, ybin], weights=weights, norm=matplotlib.colors.LogNorm())
ax.plot(peri, np.exp(line_law(np.array(np.log(peri)), *pars4)), linestyle='--', linewidth=1.5, color='black', label='d='+str(-2*pars4[1])[:4])

print(d[0].sum())

plt.xscale('log')
plt.yscale('log')
plt.grid()

plt.colorbar(d[3],ax=ax, label='Normalized Frequency')
plt.xlabel('Cloud Area (km$^2$)')
plt.ylabel('Perimeter (km)')
plt.legend()
plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloud_peri_area_density_test.png')
plt.close()


bins_list = list(np.arange(0,1050,50))
weights = np.ones_like(all['nndist']*15)/float(len(all['nndist']))
n, bins, patches = plt.hist(x=all['nndist']*15, bins=bins_list, histtype=u'step', weights=weights)
plt.grid()
plt.xlabel('NN Distance (m)')
plt.ylabel('Normalized Frequency')
plt.tight_layout()
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloud_nn_dist_test.png')
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
plt.savefig('/data/keeling/a/mdevera2/macrost/'+folder+'/nolarge_cloud_nn_dist_radius_test.png')

print(n)