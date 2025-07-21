#%% Import the libraries
import pandas as pd
import os 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.discrete.count_model import ZeroInflatedPoisson as zip
from scipy.stats.contingency import odds_ratio
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint
from scipy.stats import fisher_exact, norm
from sklearn.metrics import accuracy_score, roc_auc_score
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)
# models for the priority score
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier, BaggingClassifier
from sklearn.svm import LinearSVC, NuSVC
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import StackingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.gaussian_process import GaussianProcessClassifier

from sklearn.inspection import permutation_importance
from sklearn.inspection import PartialDependenceDisplay

#################################################################################################
# filter dinamically the probability - use .5 for all but qtl type 2 and .7 for qtl type 2
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

#%% Import disease domain
disease_domain = pd.read_csv('/container/Resources/group_EFO_by_therapeuthicaAreas.csv')
cancer_pheno = pd.read_csv('/container/Resources/cancer_phenotypes.csv')
###################################################################################################
#%%  Load drug information on genes and diseases
###################################################################################################

moleculeInfo = pd.read_parquet('/container/Resources/molecule/')
targets = pd.json_normalize(moleculeInfo['linkedTargets'])
targets.columns = ['targetId', 'targetcount']
disease = pd.json_normalize(moleculeInfo['linkedDiseases'])
disease.columns = ['diseaseId', 'diseasecount']
molInfo = moleculeInfo[['id', 'yearOfFirstApproval', 'maximumClinicalTrialPhase', 'isApproved']]
molInfo = pd.concat([molInfo, targets, disease], axis=1).dropna(subset=['targetId', 'diseaseId'])
molInfo = molInfo.explode('targetId')
molInfo = molInfo.explode('diseaseId')
molInfo.dropna(subset=['targetId', 'diseaseId'], inplace=True)
# group by gene and reatain only the row with the maximum clinical trial phase
molecules_G = molInfo.sort_values(by='maximumClinicalTrialPhase', ascending=False)\
    .drop_duplicates(subset='targetId', keep='first')\
    .dropna(subset=['targetId'])

#############################################################################################################
# Read in the AS data
#############################################################################################################
#%% import AS estimates
dir = '/container/AS_stats_out_20240325' 
files = os.listdir(dir)
files = [os.path.join(dir, i) for i in files]
# read the files in python
dfs = [pd.read_csv(i) for i in files]
# Concatenate the dataframes
#%% concatenate the dataframes 
df = pd.concat(dfs)
#%% remove not-relevant lines and filter to keep more robust evidences
df = df[df['Unnamed: 1'] == 'slope']
as_estimes = df[treshold_df(df, 0.8, 0.55)]

############################################################################################################
#%%  Number of phenotypes per gene
############################################################################################################
# Group by 'gene' and calculate the number of unique 'as_disease' entries per gene
unique_counts = as_estimes.groupby('gene')['as_disease'].nunique()
# Calculate the median of the unique counts
median_unique_counts = unique_counts.median()
iqr = unique_counts.quantile(0.75) - unique_counts.quantile(0.25)

print(f"Median number of phenotypes per gene: {median_unique_counts}")
print(f"Interquartile range: {iqr}")

############################################################################################################
#%% Calculate the priority score
################################################################################################################

# Include genes that are not in the chembl database - it is a negative control - should score the lowest
not_targets = as_estimes[~as_estimes.gene.isin(molecules_G.targetId)].sample(100, random_state=1)
# Add information on the Drug to the AS
as_estimes = as_estimes.merge(
     molecules_G[['isApproved','targetId','diseaseId', 'maximumClinicalTrialPhase','yearOfFirstApproval']], 
     left_on=['gene'], 
     right_on=['targetId'], 
     how='inner')
# Attache the gene without drugs to the main df
as_estimes = pd.concat([as_estimes,not_targets])
# Remove eventual duplicates
as_estimes = as_estimes.drop_duplicates(keep='first')

# Define the Approved status
as_estimes['isApproved'] = as_estimes.isApproved.fillna(False) # Initialise with False status
as_estimes['isApproved'] = as_estimes.maximumClinicalTrialPhase >=3 # Code to True is it is 3 or 4

# Fill in the missing values
as_estimes['maximumClinicalTrialPhase'] = as_estimes.maximumClinicalTrialPhase.fillna(0)

# Rename the targetID to gene
as_estimes['targetId'] = as_estimes.gene

#%% ############################################################################################################
#%% Train the priority score model
# X are independent data use in the to train the model. i.e. Mean AS slope, intercept, std confidences
# y is the dependent dat, in this case is / not is drug
X = as_estimes[['median', 'n_variants', 'source','qtl', 'PP_slope<0']]
X = pd.get_dummies(X, columns=['source', 'qtl']) # convert categorical variables to dummies
y = as_estimes[['isApproved','maximumClinicalTrialPhase', 'yearOfFirstApproval']].replace({True: 1, False: 0}) 

# Define the parameter distribution
param_dist = {"max_depth": [3, None],
              "max_features": randint(1, 4),
              "min_samples_split": randint(2, 11),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=0)

# Use an aggregated model to predict the priority score
estimators = [
     ('gbc', GradientBoostingClassifier(random_state=42, n_estimators=50) ),
     ('svc', make_pipeline(StandardScaler(), LinearSVC(random_state=42))),
     ('gpc', make_pipeline(StandardScaler(), GaussianProcessClassifier(random_state=42)) )
     ]

# Stack the model's classifier prediciton
clf = StackingClassifier( estimators=estimators, final_estimator=LogisticRegression() )

# Fit the priority score model to predict approved drugs
clf.fit(X_train, y_train.isApproved)

# Pedict using isApproved model on the test set
predictions = clf.predict_proba(X_test)

# Get a bit of quick summary statistics checks
accuracy = accuracy_score(y_test.isApproved, clf.predict(X_test) )                 
roc_auc = roc_auc_score(y_test.isApproved, clf.predict(X_test))
print(f"Accuracy: {accuracy}, 'ROC AUC: {roc_auc}")

############################################################################################################
#%%  Plot the distribution of the priority score
############################################################################################################

# Use the trained model to predict the probability labels of the test set
X_test['predictions'] = clf.predict_proba(X_test).T[1]
X_test[['isApproved','maximumClinicalTrialPhase', 'yearOfFirstApproval']] = y_test # Attach the test labels to the predictions
# Plot the distribution 
sns.histplot(X_test, x='predictions', hue='isApproved', kde=True, bins=40, multiple='dodge', stat='percent')
# x axis label 
plt.xlabel('Priority score')
# Save the plot
plt.savefig('/container/Results/priority_score_distribution_percent_yesNo.pdf', dpi=600, bbox_inches='tight' )

##################################################################################
#%% Expand the prediciton to whole AS df
##################################################################################

# These are the independent prediction features
as_estimesPrediction = as_estimes[['median', 'n_variants', 'source','qtl', 'PP_slope<0']]
as_estimesPrediction = pd.get_dummies(as_estimesPrediction, columns=['source', 'qtl']) # convert categorical variables to dummies

# Predict the probability of being a drug
as_estimesPrediction['prediction'] = clf.predict_proba(as_estimesPrediction).T[1]

# Fill in the information bits of the prediction df that were not used in the prediction data
as_estimesPrediction['gene'] = as_estimes.gene
as_estimesPrediction['as_disease'] = as_estimes.as_disease
as_estimesPrediction['isApproved'] = as_estimesPrediction.gene.isin( molInfo[molInfo['maximumClinicalTrialPhase'] >= 3].targetId )
as_estimesPrediction['isApproved'] = as_estimesPrediction['isApproved'].replace({False:0,True:1})

##################################################################################################################################
#%%  Plot the new prediction score
##################################################################################################################################

sns.histplot(as_estimesPrediction, x='prediction')
# Add vertical line on the threshold
plt.axvline(x=0.75, color='black', linestyle='--')
# Add arrow pointing to TYK2 gene
plt.annotate('TYK2', 
             xy=(as_estimesPrediction[as_estimesPrediction.gene == 'ENSG00000105397'].prediction.max(), 0), 
             xytext=(as_estimesPrediction[as_estimesPrediction.gene == 'ENSG00000105397'].prediction.max(), 5), 
             arrowprops=dict(facecolor='blue', shrink=0.1),
             textcoords='offset points')
# Add arrow pointing to PCSK9 gene
plt.annotate('PCSK9', 
             xy=(as_estimesPrediction[as_estimesPrediction.gene == 'ENSG00000169174'].prediction.max(), 0), 
             xytext=(as_estimesPrediction[as_estimesPrediction.gene == 'ENSG00000169174'].prediction.max(), 5), 
             arrowprops=dict(facecolor='red', shrink=0.1),
             textcoords='offset points')

# Save plot
plt.savefig('/container/Results/priority_Score_extended_to_all_genes.pdf', bbox_inches='tight')

############################################################################################################
#%% Fisher enrichment for genetic evidences based on priority score threshold
############################################################################################################

# Join the information on the AS and priority score to the molecular info
DrugPrediction = as_estimesPrediction[['prediction','gene', 'median', 'PP_slope<0', 'as_disease']].merge(molInfo[['maximumClinicalTrialPhase','isApproved','targetId']], left_on='gene', right_on='targetId', how='inner')
DrugPrediction = DrugPrediction.sort_values(by='maximumClinicalTrialPhase', ascending=False)\
    .drop_duplicates(subset='gene', keep='first')
DrugPrediction.dropna(inplace=True)
DrugPrediction['isApproved'] = DrugPrediction.isApproved.replace({True:1,False:0})

listOut = []
# With Fisher exact test get the quantile for the priority score
# get the quantile for the priority score
for treshold in np.arange(.401, 0.9, 0.1):
    # Create contingency table
    # it has to be by row
    contingency_table = [
        # With genetic evidence
        [sum(DrugPrediction[DrugPrediction.prediction >= treshold].isApproved) + 0.5, 
        sum(DrugPrediction[DrugPrediction.prediction >= treshold].isApproved == 0) + 0.5], # What become drug - What did not become Drug
        # With out genetic evidence
        [sum(moleculeInfo[moleculeInfo.linkedTargets.isna()].isApproved.dropna()) + 0.5, 
        sum(moleculeInfo[moleculeInfo.linkedTargets.isna()].isApproved == 0) + 0.5]] # What become drug - What did not become Drug
    # Calculate the odds ratio and pvalues with fisher exact
    odds_ratio, p_value = fisher_exact(contingency_table)
    # Calculate the standard error and confidence intervals
    se = np.sqrt((1/contingency_table[0][0]) + (1/contingency_table[0][1]) + (1/contingency_table[1][0]) + (1/contingency_table[1][1]))
    conf_interval = [np.exp(np.log(odds_ratio) - norm.ppf(0.975) * se), np.exp(np.log(odds_ratio) + norm.ppf(0.975) * se)]
    # attach to output
    listOut.append({"Threshold": treshold, 
                    'OR': odds_ratio,
                    'ci_U': conf_interval[0],
                    'ci_L': conf_interval[1],
                    'pval': p_value, 
                    'TP': sum(DrugPrediction[DrugPrediction.prediction >= treshold].isApproved), 
                    'TN': sum(moleculeInfo[moleculeInfo.linkedTargets.isna()].isApproved == 0), 
                    'FP': sum(DrugPrediction[DrugPrediction.prediction >= treshold].isApproved == 0), 
                    'FN': sum(moleculeInfo[moleculeInfo.linkedTargets.isna()].isApproved.dropna())})

orDF = pd.DataFrame(listOut)

############################################################################################################
#%% Forest plot the results
############################################################################################################
plt.clf()
# Create a pointplot - forsest plot in python
sns.pointplot(data=orDF, x='Threshold', y='OR', join=False)
# Add confidence interval
plt.vlines(x=orDF['Threshold'].astype('str'), ymin=orDF['ci_L'], ymax=orDF['ci_U'], alpha=0.7)
# Add a vertical line at x=0
plt.axhline(y=1, color='black', linestyle='--')
plt.xticks(rotation=20)
plt.savefig('/container/Results/odds_ratio_enrichment_priority_score.pdf', dpi=600, bbox_inches='tight')