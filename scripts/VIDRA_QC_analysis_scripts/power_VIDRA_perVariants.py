import os 
import pandas as pd
import sys
sys.path.append('/container/scripts/')
import modules_py.moduleAS as als
import numpy as np
import statsmodels.formula.api as smf
import statsmodels.api as sm
import seaborn as sns   
from matplotlib import pyplot as plt

# Read the full model data
dir = '/container/AS_stats_out/as_bayes_models_mix_robust_beta'
files = os.listdir(dir)
# paste di to file to get the full path
filesFP = [os.path.join(dir, f) for f in files]
# Read the files
data = [pd.read_csv(f) for f in filesFP]
fullModelDF = pd.concat(data)

fullModelDFwide = fullModelDF[fullModelDF['Unnamed: 0'].isin(['slope'])] \
                .pivot(index=['gene','phenotype'], columns='Unnamed: 0', values=['Mean', 'StdDev']) \
                .reset_index()

# Read the subsampled model data
dir = '/container/AS_stats_out/subsampled_as_bayes_models_mix_robust_beta'
files = os.listdir(dir)
# paste di to file to get the full path
filesFP = [os.path.join(dir, f) for f in files]
# Read the files
data = [pd.read_csv(f) for f in filesFP]
subsModelDF = pd.concat(data)

# Reorganise the subsampled model data
subsModelDFwide = subsModelDF[subsModelDF['Unnamed: 0'].isin(['slope'])] \
                .pivot(index=['gene','phenotype', 'iteration', 'fraction'], columns='Unnamed: 0', values=['Mean', 'StdDev']) \
                .reset_index()
# group by gene, phenotype and fraction and get the mean of the slopes
subsModelDFmean = subsModelDFwide\
        .drop(['iteration'],axis=1)\
        .groupby(['phenotype', 'gene', 'fraction'])\
        .mean('slope')\
        .reset_index()\
        .pivot(index=['gene', 'phenotype'], columns='fraction', values=['Mean'])\
        .reset_index()

subsModelDFsd = subsModelDFwide\
        .drop(['iteration'],axis=1)\
        .groupby(['phenotype', 'gene', 'fraction'])\
        .mean('slope')\
        .reset_index()\
        .pivot(index=['gene', 'phenotype'], columns='fraction', values=['StdDev'])\
        .reset_index()

# Merge to the full model data
meanDF = subsModelDFmean.merge(fullModelDFwide.drop('StdDev',axis=1), on=['gene', 'phenotype'])
sdDF = subsModelDFsd.merge(fullModelDFwide.drop('Mean',axis=1), on=['gene', 'phenotype'])

# pivot to longer format
pltSdDF = sdDF.melt(id_vars=[('gene', ''), ('phenotype', '')], 
                value_vars=[('StdDev', 0.25), ('StdDev', 0.5), ('StdDev', 0.75), ('StdDev', 'slope')])

pltSdDF.rename(columns={('gene', ''):'gene', ('phenotype', ''):'phenotype'}, inplace=True)

# joint hte gene and phenotype columns to create the index for colouring
pltSdDF['color'] = pltSdDF['gene'] + pltSdDF['phenotype']

ax, fig = plt.subplots(1, 1, figsize=(5,5))
# plot with dots linked by lines
# each dot is a measure, lines link gene-phenotype pairs
sns.pointplot(data=pltSdDF, 
              x='variable_1', 
              y='value',
              hue='color',
              markersize=.2,
              alpha=.2,
              legend=False, 
              palette=sns.color_palette("light:#5A9"))
sns.pointplot(data=pltSdDF, 
              x='variable_1', 
              y='value',
              markersize=1,
              color='#2986cc')
plt.savefig('/container/subsampled_AS_sd.png')
plt.clf()

# pivot to longer format
meanDFDF = meanDF.melt(id_vars=[('gene', ''), ('phenotype', '')], 
                value_vars=[('Mean', 0.25), ('Mean', 0.5), ('Mean', 0.75), ('Mean', 'slope')])

meanDFDF.rename(columns={('gene', ''):'gene', ('phenotype', ''):'phenotype'}, inplace=True)

# joint hte gene and phenotype columns to create the index for colouring
meanDFDF['color'] = meanDFDF['gene'] + meanDFDF['phenotype']
# group the df by ['gene','phenotype'] and return the entry in the 'value' column where 'variable_1' == 'slope'
slope_value = meanDFDF.groupby(['gene','phenotype']).apply(lambda x: x.loc[x.variable_1 == 'slope', 'value']).reset_index()
meanDFDF = meanDFDF.merge(slope_value, on=['gene', 'phenotype'])
value = []
for _,i in meanDFDF.iterrows():
    value.append(i.value_x - i.value_y)

meanDFDF['value'] = value
meanDFDF['abs_value'] = np.log(np.abs(value))
meanDFDF['abs_value_x'] = np.log(np.abs(meanDFDF.value_x))

# have 3 axis, one plot for each arrangement
ax, fig = plt.subplots(1, 3, figsize=(15,5))
# plot with dots linked by lines
# each dot is a measure, lines link gene-phenotype pairs
sns.pointplot(data=meanDFDF, 
              x='variable_1', 
              y='value',
              hue='color',
              ax=fig[0],
              markersize=.2,
              alpha=.2,
              legend=False, 
              palette=sns.color_palette("light:#5A9"))
sns.pointplot(data=meanDFDF, 
              x='variable_1', 
              y='value',
              ax=fig[0],
              markersize=1,
              color='#2986cc')
# title for axis 0
fig[0].set_title('Slope estimates centered on 0')
sns.pointplot(data=meanDFDF, 
              x='variable_1', 
              y='abs_value',
              hue='color',
              ax=fig[1],
              markersize=.2,
              alpha=.2,
              legend=False,
              palette=sns.color_palette("light:#5A9"))
sns.pointplot(data=meanDFDF, 
              x='variable_1', 
              y='abs_value',
              ax=fig[1],
              markersize=1,
              color='#2986cc')
# y axis lavel for axis 1
fig[1].set_ylabel('log(abs(slope))')
# title for axis 1
fig[1].set_title('Slope estimates centered on 0')
sns.pointplot(data=meanDFDF, 
              x='variable_1', 
              y='abs_value_x',
              hue='color',
              ax=fig[2],
              markersize=.2,
              alpha=.2,
              legend=False, 
              palette=sns.color_palette("light:#5A9"))
sns.pointplot(data=meanDFDF, 
              x='variable_1', 
              y='abs_value_x',
              ax=fig[2],
              markersize=1,
              color='#2986cc')
# y axis lavel for axis 2
fig[2].set_ylabel('log(abs(slope))')
# title for axis 2
fig[2].set_title('Slope estimates NOT centered')
plt.savefig('/container/subsampled_AS_MeanSlope.png')

mergedDF = pd.merge(meanDFDF, pltSdDF, on=['gene', 'phenotype', 'variable_1'], suffixes=('_mean', '_sd'))
mergedDF.variable_1.replace('slope', '1', inplace=True)
mergedDF.variable_1 = mergedDF.variable_1.astype(float)
fitted = sm.Poisson(mergedDF.value_sd, mergedDF.variable_1, missing='drop').fit(method='basinhopping')
fitted.summary()
mergedDF['color'] = mergedDF['gene'] + mergedDF['phenotype']
sns.pointplot(x=mergedDF.variable_1, 
              y=mergedDF.value_sd,
              markersize=1,
              hue=mergedDF.color,
              legend=False,
              color='#2986cc')
sns.pointplot(x=mergedDF.variable_1, 
              y=mergedDF.value_sd,
              markersize=1,
              legend=False,
              color='Blue')
sns.pointplot(x=[0.25, 0.5, 0.75, 1, 1.5, 2, 5, 10], 
              y=fitted.predict([0.25, 0.5, 0.75, 1, 1.5, 2, 5, 10]),
              markersize=1,
              color='red')
plt.savefig('/container/test_std_.png', dpi=300)
plt.clf()
