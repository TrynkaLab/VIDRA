"""

Simulation script to test the AS model

The simulation needs to start from a given initial condition (i.e. true dose response relationship) 
and then split it in 3 different noised datasets to simulate 3 different experimental conditions.

We test the AS model on different slopes, different noise levels and different number of data points (i.e. variants).

"""
#%% Imports
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from cmdstanpy import CmdStanModel
import multiprocessing
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy.stats as stats
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

def generate_random_regression(n_observations, slope, sigma, noise, meaned):
    # Generate independent variable 
    x = np.random.normal(meaned, sigma, n_observations)
    y = 0 + (slope * x + np.random.normal(0, noise, n_observations))
    return x, y

def generate_AS(slope_true, 
                n_common, n_rare, n_categorical, 
                sigma_common, sigma_rare, sigma_categorical,
                noise_common, noise_rare, noise_categorical,
                meaned_common = 1, meaned_rare = .5, meaned_categorical = 0):
    
    common_x, common_y = generate_random_regression(n_common, slope_true, sigma_common, noise_common, meaned_common)
    rare_x, rare_y = generate_random_regression(n_rare, slope_true, sigma_rare, noise_rare, meaned_rare)
    categorical_x, categorical_y = generate_random_regression(n_categorical, slope_true, sigma_categorical, noise_categorical, meaned_categorical)
    # Generate independent variable 
    common_df = pd.DataFrame({'x': common_x, 'y': common_y, 'source': 0})
    rare_df = pd.DataFrame({'x': rare_x, 'y': rare_y, 'source': 1})
    categorical_df = pd.DataFrame({'x': categorical_x, 'y': categorical_y, 'source': 2})
    return pd.concat([common_df, rare_df, categorical_df])

def conver_df_to_ModelSutableDict(df):
    return {
        'N': df.shape[0],
        # nu is the degree of freedom of the t-distribution if nu = 0 then set it to 1
        'nu': df.shape[0] - 1 if (df.shape[0] - 1) > 0 else 1,
        'numG1':df.source.values.tolist(),
        'xc':df.x.values.tolist(), 
        'yOR':df.y.values.tolist(), 
    }

def calculate_posterior(df, modelBs):
    df_dict = conver_df_to_ModelSutableDict(df)
    # Fit the model
    fit = modelBs.variational(
        data=df_dict,
        seed=412,
        algorithm='fullrank',
        grad_samples=20,
        draws=1000,
        require_converged=False,
        show_console=False,
        refresh=1000
        )
    # Save the posterior samples
    posteriors = fit.variational_sample_pd
    slope_posteriors = posteriors.iloc[:,posteriors.columns.str.startswith('slope')]
    return slope_posteriors.apply(lambda x: np.mean(x), axis = 0)

def process_data(slope, iteration, n_common, n_rare, n_categorical, 
                 sigma_common, sigma_rare, sigma_categorical, 
                 noise_common, noise_rare, noise_categorical):
    # Calculate the posterior
    df = generate_AS(slope, n_common, n_rare, n_categorical, 
                     sigma_common, sigma_rare, sigma_categorical, 
                     noise_common, noise_rare, noise_categorical)
    post = calculate_posterior(df, model).to_dict()
    post['iteration'] = iteration
    post['slope_true'] = slope
    post['n_sources'] = sum([n_common > 0, n_rare > 0, n_categorical > 0])
    post['n_common'] = n_common
    post['n_rare'] = n_rare
    post['n_categorical'] = n_categorical
    post['sigma_common'] = sigma_common
    post['sigma_rare'] = sigma_rare
    post['sigma_categorical'] = sigma_categorical
    post['noise_common'] = noise_common
    post['noise_rare'] = noise_rare
    post['noise_categorical'] = noise_categorical
    return post


def generate_AS_ss(slope_true, 
                n_common, 
                sigma_common, 
                noise_common,
                meaned_common = 1):
    
    common_x, common_y = generate_random_regression(n_common, slope_true, sigma_common, noise_common, meaned_common)
    # Generate independent variable 
    common_df = pd.DataFrame({'x': common_x, 'y': common_y, 'source': 0})
    return common_df

def process_data_ss(slope, iteration, n_common, sigma_common, noise_common):
    # Calculate the posterior
    df = generate_AS_ss(slope, n_common,
                     sigma_common,
                     noise_common)
    post = calculate_posterior(df, model_ss).to_dict()
    post['iteration'] = iteration
    post['slope_true'] = slope
    post['n_common'] = n_common
    post['sigma_common'] = sigma_common
    post['noise_common'] = noise_common
    return post

def process_data_ss_frequentist(slope, iteration, n_common, sigma_common, noise_common):
    # Calculate the posterior
    df = generate_AS_ss(slope, n_common,
                     sigma_common,
                     noise_common)
    try:
        model = smf.ols(formula='y ~ x', data=df).fit()
        post = model.params.to_dict()
        post['iteration'] = iteration
        post['slope_true'] = slope
        post['n_common'] = n_common
        post['sigma_common'] = sigma_common
        post['noise_common'] = noise_common
        return post
    except:
        post = {'slope': None,
                'iteration': iteration,
                'slope_true': slope,
                'n_common': n_common,
                'sigma_common': sigma_common,
                'noise_common': noise_common}
        return post

#%% Load the model for the simulation
# Model for the simulation
path_to_model = '/container/stan_models/simulation_hierarchical_mix_models_robust_beta_Huber.stan'
model = CmdStanModel(stan_file=path_to_model,cpp_options={'STAN_THREADS':'true'})
# Single source model 
path_to_model_ss = '/container/stan_models/simulation_single_source.stan'
model_ss = CmdStanModel(stan_file=path_to_model_ss,cpp_options={'STAN_THREADS':'true'})

# Stan parameters
chains = 4
parallel_chains = chains
threads_per_chain = 4 
#%% Simulation on multiple sources

nv_common = 3
nv_rare = 3
nv_categorical = 3
iterations = 100
max_sigma_common = 0.5
max_sigma_rare = 0.5
max_sigma_categorical = 1
max_noise_common = 0.4
max_noise_rare = 0.4
max_noise_categorical = 0.4
list_post = []  
# Number of processes to run in parallel
num_processes = multiprocessing.cpu_count()
# Create a pool of processes
pool = multiprocessing.Pool(processes=num_processes)
# Generate the list of arguments for each process
arguments = [ (slope, iteration, n_common, n_rare, n_categorical, \
                sigma_common, sigma_rare, sigma_categorical, \
                noise_common, noise_rare, noise_categorical) \
                for slope in np.arange(-2, 3, step=0.3) \
                    for iteration in range(iterations) \
                        for n_common, n_rare, n_categorical, \
                            sigma_common, sigma_rare, sigma_categorical, \
                            noise_common, noise_rare, noise_categorical in \
                            [(np.random.randint(1, nv_common),
                                np.random.randint(1, nv_rare),
                                np.random.randint(1, nv_categorical), 
                                np.random.uniform(0, max_sigma_common),
                                np.random.uniform(0, max_sigma_rare),
                                np.random.uniform(0, max_sigma_categorical),
                                np.random.uniform(0, max_noise_common),
                                np.random.uniform(0, max_noise_rare),
                                np.random.uniform(0, max_noise_categorical)) ] ]
# Run the processes in parallel
results = pool.starmap(process_data, arguments)
# Close the pool of processes
pool.close()
# Combine the results into a single list
list_post = results
# Create the DataFrame
df_simulation = pd.DataFrame(list_post)

# Remove the entrie where n_sources == 0
df_simulation = df_simulation[df_simulation['n_sources'] > 0]
df_simulation.to_csv('/container/Results/simulation_AS_model.csv', index = False)

#%% Scatterplot for the simulations
sns.scatterplot( data = df_simulation, x='slope',y='slope_true', hue = 'n_sources', alpha = 0.3, palette='viridis')
plt.savefig('/container/Results/simulation_slope_number_of_sources.pdf', bbox_inches='tight')
#%% 
df_simulation['n_variants_tot'] = df_simulation['n_common'] + df_simulation['n_rare'] + df_simulation['n_categorical']
df_simulation['n_variants_bin'] = pd.cut(df_simulation['n_variants_tot'], bins=10)
sns.scatterplot( data = df_simulation, x='slope',y='slope_true', hue = 'n_variants_bin', alpha = 0.3, palette='viridis')
plt.savefig('/container/Results/simulation_slope_number_of_variants.pdf', bbox_inches='tight')
#%% 
df_simulation['product_noise_bin'] = pd.cut(df_simulation['noise_categorical'], bins=10)
sns.scatterplot( data = df_simulation, x='slope',y='slope_true', hue = 'product_noise_bin', alpha = 0.3, palette='viridis')
plt.savefig('/container/Results/simulation_slope_noise.pdf', bbox_inches='tight')
#%% 
df_simulation['sigmas'] = df_simulation['sigma_common'] + df_simulation['sigma_rare'] + df_simulation['sigma_categorical']
df_simulation['product_sigma_categorical'] = pd.cut(df_simulation['sigmas'], bins=10)
sns.scatterplot( data = df_simulation, x='slope',y='slope_true', hue = 'product_sigma_categorical', alpha = 0.3, palette='viridis')
plt.savefig('/container/Results/simulation_slope_sigma.pdf', bbox_inches='tight')

#%% 
dep_r2_noise = df_simulation\
    .groupby('product_noise_bin')\
    .apply(lambda x: pd.Series({'r2': r2_score(x['slope_true'], x['slope'])}))\
    .reset_index()
#%% 
df_simulation['sigma_categorical_bin'] = pd.cut(df_simulation['sigma_categorical'], bins=10)
dep_r2_sigma = df_simulation\
    .groupby('sigma_categorical_bin')\
    .apply(lambda x: pd.Series({'r2': r2_score(x['slope_true'], x['slope'])}))\
    .reset_index()

dep_r2_sigma.sigma_categorical_bin = dep_r2_sigma.sigma_categorical_bin.astype(str)
sns.scatterplot(data = dep_r2_sigma, x = 'sigma_categorical_bin', y = 'r2')
plt.xticks(rotation=45)
plt.savefig('/container/Results/simulation_sigma_R2dependency.pdf', bbox_inches='tight')


############################################################################################
#%%  Single source estimate
############################################################################################

nv_common = 3
iterations = 100                          
max_sigma_common = 1
max_noise_common = 0.4
list_post = []  
# Number of processes to run in parallel
num_processes = multiprocessing.cpu_count()
# Create a pool of processes
pool = multiprocessing.Pool(processes=num_processes)
# Generate the list of arguments for each process
arguments = [ (slope, iteration, n_common, \
                sigma_common, \
                noise_common,) \
                for slope in np.arange(-2, 3, step=0.3) \
                    for iteration in range(iterations) \
                        for n_common, sigma_common, noise_common in \
                            [(np.random.randint(1, nv_common),
                                np.random.uniform(0, max_sigma_common),
                                np.random.uniform(0, max_noise_common)) ] ]
# Run the processes in parallel
results_ss = pool.starmap(process_data_ss, arguments)
# Close the pool of processes
pool.close()
# Combine the results into a single list
list_post_ss = results_ss
# Create the DataFrame
df_simulation_ss = pd.DataFrame(list_post_ss)
df_simulation_ss = df_simulation_ss[df_simulation_ss['n_common'] > 0]

##############################################################################
#%% Frequentist
##############################################################################


nv_common = 3
iterations = 100                          
max_sigma_common = 1
max_noise_common = 0.4
list_post = []  
# Number of processes to run in parallel
num_processes = multiprocessing.cpu_count()
# Create a pool of processes
pool = multiprocessing.Pool(processes=num_processes)
# Generate the list of arguments for each process
arguments = [ (slope, iteration, n_common, \
                sigma_common, \
                noise_common,) \
                for slope in np.arange(-2, 3, step=0.3) \
                    for iteration in range(iterations) \
                        for n_common, sigma_common, noise_common in \
                            [(np.random.randint(0, nv_common),
                                np.random.uniform(0, max_sigma_common),
                                np.random.uniform(0, max_noise_common)) ] ]
# Run the processes in parallel
results_ss_frequentist = pool.starmap(process_data_ss_frequentist, arguments)
# Close the pool of processes
pool.close()
# Combine the results into a single list
list_post_ss_frequentist = results_ss_frequentist
# Create the DataFrame
df_simulation_ss_frequentist = pd.DataFrame(list_post_ss_frequentist)
df_simulation_ss_frequentist.dropna(inplace=True, subset=['x', 'slope_true'])

#%% Scatterplots

fig, axs = plt.subplots(1, 3, figsize=(15, 3), sharex=False, sharey=True)

# sns.histplot( data = df_simulation_ss_frequentist, x='x', ax=axs[0,0])
# sns.histplot( data = df_simulation_ss, x='slope', ax=axs[0,1])
# sns.histplot( data = df_simulation, x='slope', ax=axs[0,2])
sns.scatterplot( data = df_simulation_ss_frequentist.sample(500), x='x',y='slope_true', ax=axs[0], alpha = 0.07)
sns.scatterplot( data = df_simulation_ss.sample(500), x='slope',y='slope_true', ax=axs[1], alpha = 0.07)
sns.scatterplot( data = df_simulation.sample(500), x='slope',y='slope_true', ax=axs[2], alpha = 0.07)
# plot x=y line on all the three plots
for ax in axs:
    ax.plot([-2, 3], [-2, 3], ls="--", c=".3") # x = y line
    ax.set_xlim(-3, 3)

plt.savefig('/container/Results/simulation_scatterplots.pdf', bbox_inches='tight', dpi=700)

#%% 
sns.scatterplot( data = df_simulation_ss_frequentist, x='x',y='slope_true', alpha = 0.07)
plt.savefig('/container/Results/simulation_frequentist_full_plot.pdf', bbox_inches='tight', dpi=700)

#%% Pearson correlation

stats.pearsonr(df_simulation_ss_frequentist['x'], df_simulation_ss_frequentist['slope_true'])
stats.pearsonr(df_simulation_ss['slope'], df_simulation_ss['slope_true'])
stats.pearsonr(df_simulation['slope'], df_simulation['slope_true'])

#%% Spearman correlation

stats.spearmanr(df_simulation_ss_frequentist['x'], df_simulation_ss_frequentist['slope_true'])
stats.spearmanr(df_simulation_ss['slope'], df_simulation_ss['slope_true'])
stats.spearmanr(df_simulation['slope'], df_simulation['slope_true'])

#%% get the accuracy metrics

df_simulation_Acc_meas = df_simulation\
    .groupby('slope_true')\
    .apply(lambda x: pd.Series({
        'mse': mean_squared_error(x['slope_true'], x['slope']),
        'rmse': np.sqrt(mean_squared_error(x['slope_true'], x['slope'])),
        'mae': mean_absolute_error(x['slope_true'], x['slope'])
    }))


#%% How many times goes in the wrong direction?

def wrong_dir(df, col1, col2):
    return sum( (df[col1] * df[col2]) < 0) / df.shape[0]

wrong_dir(df_simulation_ss_frequentist, 'x', 'slope_true')
wrong_dir(df_simulation_ss, 'slope', 'slope_true')
wrong_dir(df_simulation, 'slope', 'slope_true')