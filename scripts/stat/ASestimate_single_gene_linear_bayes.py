#%% 
import pandas as pd
from sys import argv
import sys
sys.path.append('/container/scripts/modules_py')
import moduleAS as mAS
import json
import numpy as np
import pprint
from cmdstanpy import CmdStanModel
import pickle

#%% 
h1, gene, pq_dataset = argv[1], argv[2], argv[3]

#%% 
# parse the pd dataset paths
pq_dataset_ls = pq_dataset.replace("[|]", "").split(',')

# read all the datasets
pq_list = [pd.read_parquet(pq_dataset, engine='fastparquet') for pq_dataset in pq_dataset_ls]

# if it is a list concatenate the diffferent entries of the list
df = pd.concat(pq_list)


# read pickle file
with open('/container/Resources/coloc_ref_dict.pkl', 'rb') as f:
    rename_dict = pickle.load(f)

def filter_duplecate_lead(df, lead_vars_effect_col, lead_vars_id_col, rename_dict ):
    """
    Filter duplicate lead variants from the lead_vars dataframe.
    Some duplicated variants are still in the dataframe because the lead in different GWAS may be different
    To overcome this one can use the coloc GWAS from OT. It has a list of lead variants for each GWAS and their colocalising
    """
    df[lead_vars_id_col] = df[lead_vars_id_col].replace(rename_dict)
    df.sort_values(by=lead_vars_effect_col, inplace = True)
    df.drop_duplicates(subset=lead_vars_id_col, keep='last', inplace=True)
    
    return df

# Group by the source of evidences
# and remove eventual duplicated entried due to lead variant called in different names
df = df\
        .groupby(['as_disease','GsourceLab','GqtlLab'])\
        .apply(lambda x: filter_duplecate_lead(df=x, lead_vars_effect_col='yc',lead_vars_id_col='variant', rename_dict=rename_dict))\
        .reset_index(drop=True)

###############################################################################
#%% functions
def clean_posteriorForAs(serie):
    return {'mean':np.mean(serie),
            'median':np.median(serie),
            '1%':np.percentile(serie, 1),
            '2.5%':np.percentile(serie, 2.5),
            '5%':np.percentile(serie, 5),
            '10%':np.percentile(serie, 10),
            '25%':np.percentile(serie, 25),
            '40%':np.percentile(serie, 40),
            '50%':np.percentile(serie, 50),
            '60%':np.percentile(serie, 60),
            '75%':np.percentile(serie, 75),
            '90%':np.percentile(serie, 90),
            '95%':np.percentile(serie, 95),
            '97.5%':np.percentile(serie, 97.5),
            '99%':np.percentile(serie, 99),
            'PP_slope>0':(serie > 0).mean(),
            'PP_slope<0':(serie < 0).mean()}

#%% 
###############################################################################
# Single variants estimates
def AS_singleVars(df, gene):
    df_dict = {
        'h1': float(h1),
        'N': len(df.variant),
        'numG1':df.GsourceLab.values[0],
        'numG2':df.GqtlLab.values[0],
        'xc':df.xc.values[0], 
        'xcse':df.xcse.values[0], 
        'yOR':df.yc.values[0], 
        'yORse':df.ycse.values[0], 
        'as_blosum62':df.as_blosum62.values[0],
        'as_conservation':df.as_conservation.values[0], 
        'as_sift':df.as_sift.values[0], 
        'as_polyphen':df.as_polyphen.values[0],
        'as_clinicalSignificance':df.as_clinicalSignificance.values[0], 
        'as_cadd':df.as_cadd.values[0],
        'as_alphamissense' : df.as_alphamissense.values[0],
        'as_primateai': df.as_primateai.values[0]
    }

    fit = model2.variational(
        data=df_dict,
        seed=412,
        algorithm='fullrank',
        grad_samples=20,
        draws=1000,
        require_converged=False,
        show_console=False,
        refresh=1000
        )
    
    ###############################################################################
    # save the summary file
    ###############################################################################
    # Get the ouput I want from the CmdStanVB model
    slope_posteriors = fit.stan_variable('slope', mean=False)

    # Return a pandas dataframe with posteriors
    return pd.DataFrame({
                        'gene':gene,
                        'mean':np.mean(slope_posteriors),
                        'median':np.median(slope_posteriors),
                        '1%':np.percentile(slope_posteriors, 1),
                        '2.5%':np.percentile(slope_posteriors, 2.5),
                        '5%':np.percentile(slope_posteriors, 5),
                        '10%':np.percentile(slope_posteriors, 10),
                        '25%':np.percentile(slope_posteriors, 25),
                        '40%':np.percentile(slope_posteriors, 40),
                        '50%':np.percentile(slope_posteriors, 50),
                        '60%':np.percentile(slope_posteriors, 60),
                        '75%':np.percentile(slope_posteriors, 75),
                        '90%':np.percentile(slope_posteriors, 90),
                        '95%':np.percentile(slope_posteriors, 95),
                        '97.5%':np.percentile(slope_posteriors, 97.5),
                        '99%':np.percentile(slope_posteriors, 99),
                        'PP_slope>0':(slope_posteriors > 0).mean(),
                        'PP_slope<0':(slope_posteriors < 0).mean(),
                        'n_variants':len(df.variant.unique()),
                        'source': df_dict['numG1'],
                        'qtl': df_dict['numG2'],
                        'model': 'single_variant'
    }, index=[0])

#%% 
def AS_multiVars(df, gene):
    df_dict = {
    'h1': float(h1),
    'nu': len(df.variant.unique()) - 1,
    'N': len(df.variant),
    'numG1':df.GsourceLab.tolist(),
    'numG2':df.GqtlLab.tolist(),
    'xc':df.xc.tolist(), 
    'xcse':df.xcse.tolist(), 
    'yOR':df.yc.tolist(), 
    'yORse':df.ycse.tolist(), 
    'bO':df.bO.tolist(), 
    'bOse':df.bOse.tolist(), 
    'as_blosum62':df.as_blosum62.tolist(),
    'as_conservation':df.as_conservation.tolist(), 
    'as_sift':df.as_sift.tolist(), 
    'as_polyphen':df.as_polyphen.tolist(),
    'as_clinicalSignificance':df.as_clinicalSignificance.tolist(), 
    'as_cadd':df.as_cadd.tolist(),
    'as_alphamissense' : df.as_alphamissense.tolist(),
    'as_primateai': df.as_primateai.tolist(),
    'as_revel': df.as_revel.tolist()
    }
    
    if df_dict['nu'] == 0:
        df_dict['nu'] = 1
    
    fit = model.variational(
        data=df_dict,
        seed=412,
        algorithm='fullrank',
        grad_samples=20,
        draws=1000,
        require_converged=False,
        show_console=False,
        refresh=1000
        )
    
    ###############################################################################
    # save the summary file
    ###############################################################################
    # sources of the data
    sources = [list(set(df_dict['numG1']))] # set is to return only the unique values
    # QTL type(s) of the data
    qtls = [list(set(df_dict['numG2']))]

    combination_observed = [list(pair) for pair in df[['GsourceLab', 'GqtlLab']].drop_duplicates().itertuples(index=False)]
    combination_slope = {1:[0,0], 2:[0,1], 3:[3,2], 4:[1,2], 5:[2,2]}
    slope_to_meta_analyse = []

    def get_key_for_value(dict, value_to_find):
        for key, value in dict.items():
            if value == value_to_find:
                return key
        return None

    for combination in combination_observed:
        key = get_key_for_value(combination_slope, combination)
        if key is not None:
            slope_to_meta_analyse.append(key)
    
    list_posteriors = []
    for slope in slope_to_meta_analyse:
        list_posteriors += fit.variational_sample_pd[f'slope_random[{slope}]'].tolist()
    metaPP = pd.DataFrame(clean_posteriorForAs(pd.Series(list_posteriors)),index=['meta_slope'])

    # Get the ouput from the CmdStanVB model
    posteriors = fit.variational_sample_pd
    #  Create a list of column names based on the indices
    colNamesKeep = posteriors.columns[posteriors.columns.str.startswith(('slope', 'slope_random[', 'intercept'))]
    # subset posterior to get just the slopes
    pp = pd.DataFrame(posteriors.loc[:, colNamesKeep])
    parsed_pp = pp.apply(lambda x:clean_posteriorForAs(x))
    pp_indexes = parsed_pp.index
    output = pd.json_normalize(pd.DataFrame(parsed_pp)[0]).set_index(pp_indexes)
    output.loc[len(output.index)] = metaPP.loc['meta_slope'] = metaPP.iloc[0].to_list()
    output = output.rename(index={6:'metaPP'})
    output['gene'] = gene
    output['n_variants'] = len(df.variant.unique())
    output['source'] = sources*len(output) 
    output['qtl'] = qtls*len(output)
    output['model'] = 'multiple_variant'
    return output

#%% 
def as_error_report(df, gene):
    return pd.DataFrame({
                        'gene':gene,
                        'n_variants':len(df.variant.unique()),
                        'model': 'error_report',
                        'source': [list(set(df['GsourceLab']))], # set is to return only the unique values
                        'qtl': [list(set(df['GqtlLab']))]
    }, index=[0])

#%%
def fitASmodels(df, gene):
    if len(df['variant']) > 1:
        try:
            return AS_multiVars(df, gene)
        except:
            return as_error_report(df, gene)
    else:
        try:
            return AS_singleVars(df, gene)
        except:
            return as_error_report(df, gene)

#%% 
###############################################################################

# Model multi-variant estimates
path_to_model = '/container/stan_models/VIDRA.stan'
model = CmdStanModel(stan_file=path_to_model,cpp_options={'STAN_THREADS':'true'})

# Model single-variant estimates
path_to_model2 = '/container/stan_models/VIDRA_single_variant.stan'
model2 = CmdStanModel(stan_file=path_to_model2,cpp_options={'STAN_THREADS':'true'})

# Stan parameters
chains = 4
parallel_chains = chains
threads_per_chain = 4 

#%% 
df = df.groupby(['as_disease', 'GsourceLab', 'GqtlLab']).filter(lambda x: False if ( (len(x['variant']) == 1 )&(x['GsourceLab']==3).all() ) else True )
result = df \
    .groupby('as_disease')\
    .apply(lambda x: fitASmodels(x, gene) )

#%% 
result.to_csv(f'./{gene}_Bayes_hierarchical.csv', index=True)