#%% 
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

def treshold_df(df, treshold1 = 0.7, threshold2 = 0.5):
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

dir = '/container/AS_stats_out_20240325' # this is without AZ
files = os.listdir(dir)
files = [os.path.join(dir, i) for i in files]
# read the files in python
dfs = [pd.read_csv(i) for i in files]
# Concatenate the dataframes
#%% 
df = pd.concat(dfs)
#%%
gen = df[df['Unnamed: 1'] == 'slope']
gen = gen[treshold_df(gen, treshold1 = 0.80, threshold2 = 0.55)]



