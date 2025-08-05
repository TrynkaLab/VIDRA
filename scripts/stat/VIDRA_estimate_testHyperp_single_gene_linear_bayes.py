import pandas as pd
from sys import argv
import sys
sys.path.append('/container/scripts/modules_py')
import moduleAS as mAS
import json
import numpy as np
import pprint
from cmdstanpy import CmdStanModel

# cmdstanpy.set_cmdstan_path('.cmdstan/cmdstan-2.33.0/')
h1, gene, phenotype, pq_dataset = argv[1], argv[2], argv[3], argv[4]

# parse the pd dataset paths
pq_dataset_ls = pq_dataset.replace("[|]", "").split(',')
# read all the datasets
pq_list = [pd.read_parquet(pq_dataset) for pq_dataset in pq_dataset_ls]
# if it is a list concatenate the diffferent entries of the list
df = pd.concat(pq_list)

######### Check if there is data for this gene. It make sense to run only if enough variants

if df.empty:
    print('No data for this gene')
    exit()
if len(df.variantId.unique()) <= 1:
    print('Only less than two independent variants for this gene-phenotype combo')
    exit()
# if df.shape[0] <= 2:
#     print('Only less than two variant for this gene')
#     exit()
print(f'More than two variants for this gene, i.e. {df.shape[0]}')

###############################################################################
path_to_model = '/container/stan_models/VIDRA.stan'
model = CmdStanModel(stan_file=path_to_model,cpp_options={'STAN_THREADS':'true'})

chains = 4
parallel_chains = chains
threads_per_chain = 4 # floor( ( multiprocessing.cpu_count() - 2 ) / chains)

df_dict = {
    'h1': float(h1),
    'nu': len(df.variantId.unique()) - 1,
    'N': len(df.variantId),
    'numG1':df.GsourceLab.tolist(),
    'numG2':df.GqtlLab.tolist(),
    'xc':df.xc.tolist(), 
    'xcse':df.xcse.tolist(), 
    'yOR':df.yc.tolist(), 
    'yORse':df.ycse.tolist(), 
    'bO':df.bO.tolist(), 
    'bOse':df.bOse.tolist(), 
    'as_blosum62':df.as_blosum62.tolist(),
    'as_foldx':df.as_foldx.tolist(), 
    'as_plddt':df.as_plddt.tolist(), 
    'as_conservation':df.as_conservation.tolist(), 
    'as_sift':df.as_sift.tolist(), 
    'as_polyphen':df.as_polyphen.tolist(),
    'as_clinicalSignificance':df.as_clinicalSignificance.tolist(), 
    'as_cadd':df.as_cadd.tolist(),
    'as_alphamissense' : df.as_alphamissense.tolist(),
    'as_consequence': df.as_consequence.tolist(),
    'as_primateai': df.as_primateai.tolist()
}

try:
    optimized = model.optimize(data=df_dict)
    print("Optimization was successful")
    fit = model.sample(
        data=df_dict,
        chains=chains, 
        parallel_chains=parallel_chains,
        threads_per_chain=threads_per_chain,
        thin=10,
        max_treedepth=15,
        adapt_step_size=500,
        adapt_metric_window=500,
        adapt_init_phase=500,
        adapt_delta=0.85,
        iter_warmup=1000,
        iter_sampling=1500,
        seed=412,
        inits=optimized.optimized_params_dict,
        show_progress=False
        )
except RuntimeError:
    print("Optimization was not successful")
    fit = model.sample(
        data=df_dict,
        chains=chains, 
        parallel_chains=parallel_chains,
        threads_per_chain=threads_per_chain,
        thin=10,
        max_treedepth=15,
        adapt_step_size=500,
        adapt_metric_window=500,
        adapt_init_phase=500,
        adapt_delta=0.85,
        iter_warmup=1000,
        iter_sampling=1500,
        seed=412,
        show_progress=False
        )

print(f'divergences:\n{fit.divergences}\niterations at max_treedepth:\n{fit.max_treedepths}')
print(fit.diagnose())

# Return a pandas dataframe with posteriors
summary_df = fit.summary(percentiles=(1, 2.5, 5, 10, 25, 40, 50, 60, 75, 90, 95, 97.5, 99), sig_figs=3) 
# add the phenotype and gene columns
summary_df['phenotype'] = phenotype 
summary_df['gene'] =  gene 
summary_df['h1'] =  h1
summary_df['n_variants'] = len(df.variantId.unique())


# save the summary file
summary_df.to_csv(f'./{phenotype}_{gene}_hype{h1}_Bayes_testHyperp.csv')