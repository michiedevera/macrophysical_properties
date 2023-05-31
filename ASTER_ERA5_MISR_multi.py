import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from scipy.optimize import curve_fit, fsolve
from scipy import stats as st
import netCDF4 as nc
from statistics import mode, mean, stdev
import matplotlib.dates as mdates
from datetime import datetime
from sklearn import linear_model
import statsmodels.api as sm
import seaborn as sns
from statsmodels.stats.outliers_influence import variance_inflation_factor

df = pd.read_csv('/data/keeling/a/mdevera2/macrost/df_nolarge_CAMP2Ex_ASTER_ERA5_MISR_vars.csv')

properties = ['line', 'direct', 'fractal_dim', 'height_mode', 'height_mean', 'cf']

for prop in properties:
    x = df[['sst','tcwv','uv10','wdir10','rh1000','rh975','rh950','rh925','rh900','rh875','rh850',
    'w1000','w975','w950','w925','w900','w875','w850','u1000','u975','u950','u925','u900','u875','u850','v1000','v975','v950','v925','v900','v875','v850',
    'uv1000','uv975','uv950','uv925','uv900','uv875','uv850','wdir1000','wdir975','wdir950','wdir925','wdir900','wdir875','wdir850','u600', 'v600', 'uv600', 'wdir600', 
    'LTS700', 'LTS850', 'LTS875', 'LTS900', 'LTS925', 'LTS950', 'LTS975', 'EIS700',
    'thetae1000','thetae975','thetae950','thetae925','thetae900','thetae875','thetae850', 'CAPE850', 'CIN850',
    'uv1000-975', 'uv1000-950', 'uv1000-925', 'uv1000-900', 'uv1000-875', 'uv1000-850', 'uv950-925', 'uv950-900', 'uv950-875', 'uv950-850', 'uv850-600',
    'thetae1000-975', 'thetae1000-950', 'thetae1000-925', 'thetae1000-900', 'thetae1000-875', 'thetae1000-850','thetae950-925', 'thetae950-900', 'thetae950-875', 'thetae950-850']]

    x_std = st.zscore(x)
    y = df[prop]

    # with statsmodels
    x_std = sm.add_constant(x_std) # adding a constant

    model = sm.OLS(y, x_std).fit()

    if (np.argmax(model.pvalues) == 0):
        max_p = model.pvalues.sort_values()[-2]
    else:
        max_p = np.max(model.pvalues)

    while max_p > 0.05:
        idx = np.where(model.pvalues == max_p)[0][0]
        x = x.drop(x.columns[idx-1], axis=1)
        x_std = st.zscore(x)
        y = df[prop]

        # with statsmodels
        x_std = sm.add_constant(x_std) # adding a constant
        
        model = sm.OLS(y, x_std).fit()

        if (np.argmax(model.pvalues) == 0):
            max_p = model.pvalues.sort_values()[-2]
        else:
            max_p = np.max(model.pvalues)

    with open('/data/keeling/a/mdevera2/macrost/multi_'+prop+'_std_auto.txt', 'w') as fh:
        fh.write(model.summary().as_text())
    fh.close()

    r2 = []
    adj_r2 = []

    for i in range(len(x.columns)):
        x_std = st.zscore(x.drop(x.columns[i], axis=1))
        y = df[prop]

        # with statsmodels
        x_std = sm.add_constant(x_std) # adding a constant
        
        model = sm.OLS(y, x_std).fit()
        r2.append(model.rsquared)
        adj_r2.append(model.rsquared_adj)

    res = sorted(range(len([i - model.rsquared for i in r2])), key = lambda sub: [i - model.rsquared for i in r2][sub])[:6]
    res_adj = sorted(range(len([i - model.rsquared_adj for i in adj_r2])), key = lambda sub: [i - model.rsquared_adj for i in adj_r2][sub])[:6]

    top_r2 = x.columns[res]
    top_r2_change = []
    for idx in res:
        top_r2_change.append([i - model.rsquared for i in r2][idx])

    top_r2_adj = x.columns[res_adj]
    top_r2_adj_change = []
    for idx in res_adj:
        top_r2_adj_change.append([i - model.rsquared_adj for i in adj_r2][idx])

    psort = sorted(range(len(model.pvalues[1:])), key = lambda sub: model.pvalues[1:][sub])[:6]
    top_pval = x.columns[psort]

    df_rank = pd.DataFrame(zip(top_r2, top_r2_change, top_r2_adj, top_r2_adj_change, top_pval), columns=['rank_r2', 'change_r2', 'rank_r2_adj', 'change_r2_adj', 'rank_pval'])
    df_rank.to_csv('/data/keeling/a/mdevera2/macrost/multi_'+prop+'_std_auto_rank.csv', index=False)