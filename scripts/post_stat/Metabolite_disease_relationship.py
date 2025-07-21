#%% 
# import modules
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import seaborn as sns
import ast
from sklearn.mixture import GaussianMixture
from sklearn.linear_model import TheilSenRegressor, LinearRegression
from sklearn.metrics import mean_squared_error
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

#%%  Read metabolites and diseases tuples
# Define traits and metabolites
# The list traits contains phenotypes and metabolites
# the metabolites list has only the metabolites
# the traits can be obtained subtracting metabolites from traits
# The first split command devide the lines in 2 and take the one before the comment
# the if skip the lines that start with #
with open('/container/Resources/metabolite_disease_pairs.txt', 'r') as f:
    MetDis = [ast.literal_eval(line.split('#')[0]) for line in f if not line.strip().startswith('#')]

#---------------------------------------
#%% 
# read in the AS data
dir = "/container/AS_stats_out_20240325/"
files = os.listdir(dir)
pathFiles = [os.path.join(dir, f) for f in files]
df = pd.concat([pd.read_csv(f) for f in pathFiles])
# select only the slopes
df = df[ df['Unnamed: 1'] == 'slope']

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

df = df[treshold_df(df, .8, .55)]

#%% 
# Create diseases and metabolites df
# for every metabolite disease value pair in the MetDis list look at what metabolites are in the df and what disease are in the df
# if the metabolite is in the df then then check also the disease is in the df
# then for each metabolite in the df - and each disease in the df - (there could be multiple associated with differnt genes) get the mean values
# concatenate these two df by column
def metabolite_disease_matching(df, MetDis):
    listout = []
    for met, dis in MetDis:
        if met in df['as_disease'].values:
            if dis in df['as_disease'].values:
                for _, metabolite in df[df.as_disease == met].iterrows():
                    for _, disease in df[df.as_disease == dis].iterrows():
                        if metabolite['gene'] == disease['gene']:
                            listout.append(
                                pd.DataFrame(
                                    {'metabolite': metabolite['as_disease'],
                                     'disease': disease['as_disease'],
                                     'gene': metabolite['gene'],
                                     'Mean_metabolite': metabolite['mean'],
                                     'Mean_disease': disease['mean'],
                                     'conf_metabolite': metabolite['PP_slope<0'],
                                     'conf_disease': disease['PP_slope<0'],},
                                    index=[0]))
            else:
                continue
        else:
            continue
    return pd.concat(listout)


#%% Get the metabolites and diseases pairs
MetDiseDF = metabolite_disease_matching(df, MetDis)

#%% 
# define the quadrant of the plot
sign = []
for _,line in MetDiseDF.iterrows():
    if (( float(line['Mean_metabolite']) > 0) & ( float(line['Mean_disease']) > 0)):
        sign.append("++")
    if (( float(line['Mean_metabolite']) < 0) & ( float(line['Mean_disease']) < 0)):
        sign.append("--")
    if (( float(line['Mean_metabolite']) > 0) & ( float(line['Mean_disease']) < 0)):
        sign.append("+-")
    if (( float(line['Mean_metabolite']) < 0) & ( float(line['Mean_disease']) > 0)):
        sign.append("-+")

# add the sign to the df
MetDiseDF['sign'] = sign
sign_general = []
for _,line in MetDiseDF.iterrows():
    if (( float(line['Mean_metabolite']) > 0) & ( float(line['Mean_disease']) > 0)):
        sign_general.append("+")
    if (( float(line['Mean_metabolite']) < 0) & ( float(line['Mean_disease']) < 0)):
        sign_general.append("+")
    if (( float(line['Mean_metabolite']) > 0) & ( float(line['Mean_disease']) < 0)):
        sign_general.append("-")
    if (( float(line['Mean_metabolite']) < 0) & ( float(line['Mean_disease']) > 0)):
        sign_general.append("-")

# add the sign to the df
MetDiseDF['sign_general'] = sign_general
#%% 
MetDiseDF[['Mean_metabolite', 'Mean_disease']] = MetDiseDF[['Mean_metabolite', 'Mean_disease']].apply(pd.to_numeric)
#%%
#%% K-means clustering of the data
singArray, signIndex = pd.factorize(MetDiseDF['sign'])
MetDiseDF['sign'] = singArray.tolist()
X = MetDiseDF[['Mean_metabolite','Mean_disease', 'sign']].to_numpy()
y_pred = GaussianMixture(n_components=4, init_params='kmeans', random_state=64).fit_predict(X)
MetDiseDF['sign'] = signIndex[singArray]
MetDiseDF['clusters'] = y_pred.astype('str')
#%% 
# calculate the robust regression
def robust_regression(df):
    X = df['Mean_metabolite'].values.reshape(-1, 1)
    y = df['Mean_disease'].values
    model = TheilSenRegressor()
    model.fit(X, y)
    # get the slope and intercept
    slope = model.coef_
    return (slope)

observed_slope = MetDiseDF.groupby(['sign','clusters']).apply(robust_regression)
#%%
# facet the plot according to the sign
# free x and y axis
g = sns.FacetGrid(MetDiseDF, col='sign_general', sharex=False, sharey=True, hue='sign_general', 
                  hue_order=['+','-'], col_order=['+','-'] )
# then make a line plot that links same gene
g.map(sns.scatterplot, 'Mean_metabolite', 'Mean_disease', alpha=0.5)
g.map(sns.regplot, 'Mean_metabolite', 'Mean_disease', scatter=False, robust=True)
plt.show()
#%% 
g2 = sns.FacetGrid(MetDiseDF, col='sign', sharex=False, sharey=True, hue='sign', 
                   hue_order=['++','--','+-','-+'],  col_order=['++','--','+-','-+'] )
# plot the dots and then add the regression using the betas calculated in the observed_slope
g2.map(sns.scatterplot, 'Mean_metabolite', 'Mean_disease', alpha=0.5)
g2.map(sns.regplot, 'Mean_metabolite', 'Mean_disease', robust=True) # regression values
# g2.set(xlim=(-1.5, 1)) # set x axis limits
plt.show()

#%% 
# create random pairs
def random_pairs(df, n_pairs=100):
    # shuffle the df
    df = df.sample(n=n_pairs*2).reset_index(drop=True)
    # get the first half of the df
    df1 = df.sample(n=n_pairs)
    # get the second half of the df
    df2 = df[~df.isin(df1)].dropna(how='all').sample(n=n_pairs)
    df3 = pd.DataFrame({'Mean_x': df1['mean'].reset_index(drop=True), 
                        'Mean_y': df2['mean'].reset_index(drop=True)})
    # define the sings
    sign = []
    for _,line in df3.iterrows():
        if (( float(line['Mean_x']) > 0) & ( float(line['Mean_y']) > 0)):
            sign.append("++")
        if (( float(line['Mean_x']) < 0) & ( float(line['Mean_y']) < 0)):
            sign.append("--")
        if (( float(line['Mean_x']) > 0) & ( float(line['Mean_y']) < 0)):
            sign.append("+-")
        if (( float(line['Mean_x']) < 0) & ( float(line['Mean_y']) > 0)):
            sign.append("-+")
    # add the sign to the df
    df3['sign'] = sign
    return df3

#%% Perform the bootstrap
# function for robust regression

def robust_regression(df):
    X = df['Mean_x'].values.reshape(-1, 1)
    y = df['Mean_y'].values
    model = TheilSenRegressor()
    model.fit(X, y)
    slope = model.coef_
    return (slope) 

n_bootstraps = 10000
bootstrap_slopes_pp = np.empty(n_bootstraps)
bootstrap_slopes_nn = np.empty(n_bootstraps)
bootstrap_slopes_np = np.empty(n_bootstraps)
bootstrap_slopes_pn = np.empty(n_bootstraps)

for i in range(n_bootstraps):
    bootstrap_sample = random_pairs(df = df, n_pairs= MetDiseDF.shape[0])
    bootstrap_sample['Mean_x'] = bootstrap_sample['Mean_x'].astype('float')
    bootstrap_sample['Mean_y'] = bootstrap_sample['Mean_y'].astype('float')   
    try:
        bootstrap_model = bootstrap_sample.groupby('sign',).apply(robust_regression)
        bootstrap_slopes_pp[i] = bootstrap_model[0]
        bootstrap_slopes_pn[i] = bootstrap_model[1]
        bootstrap_slopes_np[i] = bootstrap_model[2]
        bootstrap_slopes_nn[i] = bootstrap_model[3]
    except:
        continue

#%% Calculate the p-values between observed and bootstrap
observed_slope_pp = observed_slope[0]
observed_slope_pn = observed_slope[1]
observed_slope_np = observed_slope[2]
observed_slope_nn = observed_slope[3]

# round to the 4th decimal
p_value_pp = (np.sum(bootstrap_slopes_pp > observed_slope_pp) + 1/n_bootstraps) / n_bootstraps
p_value_pn = (np.sum(bootstrap_slopes_pn < observed_slope_pn) + 1/n_bootstraps) / n_bootstraps
p_value_np = (np.sum(bootstrap_slopes_np < observed_slope_np) + 1/n_bootstraps) / n_bootstraps
p_value_nn = (np.sum(bootstrap_slopes_nn > observed_slope_nn) + 1/n_bootstraps) / n_bootstraps
#%% 
# plot bootstrap
width_px = 8000
height_px = 1700

# DPI for your system. For most systems, it's 96, 100, or 120. You may need to adjust this value.
dpi = 600

# Calculate size in inches
width_in = width_px / dpi
height_in = height_px / dpi

bootDF = pd.DataFrame({'pp': bootstrap_slopes_pp, 'nn': bootstrap_slopes_nn, 
                       'pn': bootstrap_slopes_pn, 'np': bootstrap_slopes_np})
# create subplot with 4 axes (4 columns) and plot the 4 distributions, one per axis
# the size of the plot should be 4 times the size of a single plot
fig, axs = plt.subplots(1, 4, figsize=(width_in, height_in))
sns.histplot(data=bootDF, x = 'pp', fill=True, ax=axs[0], stat='probability')
# add pvalue estimate onto the plot
axs[0].axvline(x=observed_slope_pp, color='red')
axs[0].set_title('sign = ++')
axs[0].set(ylim = (0, 0.12))
sns.histplot(data=bootDF, x = 'nn', fill=True, ax=axs[1], stat='probability')
axs[1].axvline(x=observed_slope_nn, color='red')
axs[1].set_title('sign = --')
axs[1].set(ylabel='', ylim = (0, 0.12))
sns.histplot(data=bootDF, x = 'pn', fill=True, ax=axs[2], stat='probability')
axs[2].axvline(x=observed_slope_pn, color='red')
axs[2].set_title('sign = +-')
axs[2].set(ylabel='', ylim = (0, 0.12))
sns.histplot(data=bootDF, x = 'np', fill=True, ax=axs[3], stat='probability')
axs[3].axvline(x=observed_slope_np, color='red')
axs[3].set_title('sign = -+')
axs[3].set(ylabel='', ylim = (0, 0.12))
plt.show()
