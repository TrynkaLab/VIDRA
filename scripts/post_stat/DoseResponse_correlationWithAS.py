#%% 
import os
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import seaborn as sns
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import pearsonr, spearmanr
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

#%% 
# import clinical trials data
clinical_trials = pd.read_csv('/container/Resources/ManuallyCurated_DoseResponseClinicalTrials.csv', encoding='latin1')
#%% Calculate the regression oin the clinical trials data

clinical_trials_slopes = clinical_trials\
    .groupby(['study', 'GeneId', 'PhenotypeId'])\
    .apply(lambda x: sm.RLM(endog=x['Response outcome'], exog=x['Drug concentration']).fit().params,
                    #    sm.OLS(endog=x['Response outcome'], exog=x['Drug concentration']).fit().pvalues 
                    ) \
    .reset_index()\
    .rename(columns={'Drug concentration': "stats_regression"})

#%% import AS estimates
dir_path = "/container/AS_stats_out_20240325/"
AS = [i for i in os.listdir(dir_path)]
AS = [os.path.join(dir_path, i) for i in AS]
as_estimes = [pd.read_csv(i) for i in AS]
as_estimes = pd.concat(as_estimes)
as_estimes['method'] = 'AS'
as_estimes = as_estimes[as_estimes['Unnamed: 1'] == 'slope']
as_estimes = as_estimes.dropna(axis=1, how='all')
#%% 
# Select only confident ones
def treshold_df(df, treshold1 = 0.5, threshold2 = 0.5):
    filter_list = []
    for _,row in df.iterrows():
        if ( (row['qtl'] == '2')|(row['qtl'] == '[2]') & ( (row['PP_slope<0'] > treshold1) | (row['PP_slope>0'] > treshold1) ) ):
                filter_list.append(True)
        else:
            if ( (row['qtl'] != '2')&(row['qtl'] != '[2]') & ( (row['PP_slope<0'] > threshold2) | (row['PP_slope>0'] > threshold2) ) ):
                    filter_list.append(True)
            else:
                filter_list.append(False)
    return filter_list

as_estimes = as_estimes[treshold_df(as_estimes, .8, .5)]
#%% Merge clinical trials and AS estimates
merged = pd.merge(clinical_trials_slopes, as_estimes, left_on=['PhenotypeId', 'GeneId'], right_on=['as_disease', 'gene'], how='inner')
# Taking only the significant pvalue may be too stringent
# group by phenotype, gene and take only the lowest p-value
# merged = merged.groupby(['phenotype','gene']).apply(lambda x: x[x['p-value'] == x['p-value'].min()])
merged.reset_index(drop=True)
merged.drop_duplicates(inplace=True)
# correlation coefficient between Mean and effect
pearsonr(merged['stats_regression'], merged['mean'])
spearmanr(merged['stats_regression'], merged['mean'])

#%%  formally test it with regression and robust regression
# Ordinary least squares regression
model = smf.ols('stats_regression ~ mean', data=merged).fit()
model.summary()
# Robust regression
model2 = smf.rlm('stats_regression ~ mean', data=merged).fit()
model2.summary()
merged['predictY'] = model2.predict(exog=merged['mean'])
# round to the 3rd digit
slope_estimate = np.around( model2.params['mean'], 3)
pval = np.around( model2.pvalues['mean'], 3)
#%% visualise
# plot bootstrap
width_px = 2500
height_px = 1500

# DPI for your system. For most systems, it's 96, 100, or 120. You may need to adjust this value.
dpi = 600

# Calculate size in inches
width_in = width_px / dpi
height_in = height_px / dpi

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
fig = plt.figure(figsize=(width_in, height_in))
sns.set_theme(style="ticks", rc=custom_params)
sns.regplot(data=merged, x='mean', y = 'stats_regression', scatter=True )
plt.text(-1.25,2, f"beta coefficient: {slope_estimate} \np-value: {pval}")
# add vertical and horizontal lines
# plt.axvline(0, color='red')
# plt.axhline(0, color='red')
plt.savefig('/container/Results/AS_Drug_dosing_correlation.pdf', dpi=dpi)
