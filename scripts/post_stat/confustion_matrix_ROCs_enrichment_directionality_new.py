#%% 
import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
from scipy.stats import t, norm, f
import scipy.stats as stats
from sklearn.preprocessing import LabelBinarizer
from itertools import cycle
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.metrics import roc_curve, auc
from scipy.stats import chi2_contingency
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from scipy.stats import binomtest

#################################################################################################
# function to filter the drug genetics dataframe
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

#################################################################################################
#%% Import disease domain
disease_domain = pd.read_csv('/container/Resources/group_EFO_by_therapeuthicaAreas.csv')
cancer_pheno = pd.read_csv('/container/Resources/cancer_phenotypes.csv')

# Load the data
dir = '/container/AS_stats_out_20240325/' 
files = os.listdir(dir)
files = [os.path.join(dir, i) for i in files]

# read the files in python
dfs = [pd.read_csv(i) for i in files]

# Concatenate the dataframes
#%% 
df = pd.concat(dfs)

df['PP_slope<0'] = df['PP_slope<0'].where(~df['as_disease'].isin(cancer_pheno.efo), 1 - df['PP_slope<0'])
df['PP_slope>0'] = df['PP_slope>0'].where(~df['as_disease'].isin(cancer_pheno.efo), 1 - df['PP_slope>0'])
df['mean'] = df['mean'].where(~df['as_disease'].isin(cancer_pheno.efo), -1 * df['mean'])

#%% 
df = df[df['Unnamed: 1'] == 'slope']
df = df[treshold_df(df, treshold1 = 0.80, threshold2 = 0.55)]
df['MoA_AS'] = df['mean'].apply(lambda x: 'inhibitor' if x > 0 else 'activator')

#%% 
# get credible intervals consistent with the Mean directionality
# find all the columns that contain '%' and return the latest the return the positive ones
def find_credible_interval(line):
    colOut = 0
    for col in df.filter(like='%').columns:
        if (line['mean'] * line[col] >= 0):
            colOut += 1
    return(colOut)

#%% 
df['confidence'] = df.apply(lambda line: find_credible_interval(line), axis=1)
df['confidence_norm'] = df.apply(lambda line: find_credible_interval(line), axis=1) / df.filter(like='%').columns.size

#%% 
# ----------------------------------------------------------------------------------
# Drug information
# ----------------------------------------------------------------------------------
################################################################################################
drugInfo = pd.read_parquet('/container/Resources/mechanismOfAction/')
# remove reference info which I won't use in here
drugInfo.drop(columns=['references'], inplace=True)

# Keep only drugs that target a single protein (i.e. target)
removeTypes = ['protein family', 'protein complex', 'protein complex group', 'chimeric protein', 'nucleic-acid', 'selectivity group']
drugInfo = drugInfo[~drugInfo.targetType.isin(removeTypes)]
# pass from column with list to a one entry per row in a long format
drugInfo = drugInfo.explode('targets')
drugInfo = drugInfo.explode('chemblIds')

# remove rows with no target or chemblId
drugInfo.dropna(subset=['targets', 'chemblIds'], inplace=True)

# info on the molecule
moleculeInfo = pd.read_parquet('/container/Resources/molecule/')
targets = pd.json_normalize(moleculeInfo['linkedTargets'])
targets.columns = ['targetId', 'targetcount']
disease = pd.json_normalize(moleculeInfo['linkedDiseases'])
disease.columns = ['diseaseId', 'diseasecount']
molInfo = moleculeInfo[['id', 'yearOfFirstApproval', 'maximumClinicalTrialPhase']]
molInfo = pd.concat([molInfo, targets, disease], axis=1).dropna(subset=['targetId', 'diseaseId'])
molInfo = molInfo.explode('targetId')
molInfo = molInfo.explode('diseaseId')
molInfo.dropna(subset=['targetId', 'diseaseId'], inplace=True)

# split targets column, which has lists of targets, into multiple columns
molInfoMoA = molInfo.merge(drugInfo, left_on=['id', 'targetId'], right_on=['chemblIds', 'targets'], how='inner')
drugInfoHI = molInfoMoA

# categorise the mode of action in acrivators, inhibitors and unknown
# activators id the description contains 'agonist' or 'activator' or 'positive modulator' or 'stabiliser'
# inhibitors if the description contains 'antagonist' or 'inhibitor' or 'degrader'
# unknown if the description does not contain any of the above
inhibitors = [
    "RNA INHIBITOR",
    "NEGATIVE MODULATOR",
    "NEGATIVE ALLOSTERIC MODULATOR",
    "ANTAGONIST",
    "ANTISENSE INHIBITOR",
    "BLOCKER",
    "INHIBITOR",
    "DEGRADER",
    "INVERSE AGONIST",
    "ALLOSTERIC ANTAGONIST",
    "DISRUPTING AGENT",
]

activators = [
    "PARTIAL AGONIST",
    "ACTIVATOR",
    "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR",
    "AGONIST",
    "SEQUESTERING AGENT",
    "STABILISER", 
]

MoA = []

for index, row in drugInfoHI.iterrows():
    if any(keyword in row.actionType.upper() for keyword in activators):
        MoA.append('activator')
    elif any(keyword in row.actionType.upper() for keyword in inhibitors):
        MoA.append('inhibitor')
    else:
        MoA.append('other')

drugInfoHI['MoA'] = MoA

#%% 
################################################################################################
# combine the drug information with the allelic series information
consistency_DF = drugInfoHI[['diseaseId', 'id', 'mechanismOfAction', 'maximumClinicalTrialPhase', 'targetId', 'MoA']] \
                    .merge( df, 
                            left_on=['diseaseId', 'targetId'], 
                            right_on=['as_disease', 'gene'],
                            how = 'inner' ) 

# ----------------------------------------------------------------------------------
# Roc curve for classifier absed on the mean of the allelic series estimates
# ----------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 6))
classes_of_interest = ['activator','inhibitor']
colors = ["darkorange", "cornflowerblue"]
for class_of_interest, color in zip(classes_of_interest, colors):
    consistencyMoAdummy = pd.get_dummies(consistency_DF.MoA)
    pp = consistency_DF['PP_slope>0'] if class_of_interest == 'inhibitor' else consistency_DF['PP_slope<0'] 
    display = RocCurveDisplay.from_predictions(
        consistencyMoAdummy[[class_of_interest]],
        pp,
        name=class_of_interest,
        color=color,
        ax=ax,
        plot_chance_level=(class_of_interest == 'inhibitor')
    )
    _ = ax.set(
    xlabel="False Positive Rate",
    ylabel="True Positive Rate"
    )

fig.savefig('/container/Results/auroc_mechanismofaction.pdf', dpi = 600, bbox_inches='tight')

################################################################################################
#%% confusion matrix
################################################################################################

cm = confusion_matrix(consistency_DF[consistency_DF.MoA != 'other'].MoA, 
                      consistency_DF[consistency_DF.MoA != 'other'].MoA_AS, 
                      labels=consistency_DF[consistency_DF.MoA != 'other'].MoA.unique() )
disp = ConfusionMatrixDisplay(confusion_matrix=cm, 
                              display_labels=consistency_DF[consistency_DF.MoA != 'other'].MoA.unique())
disp.plot(cmap='Blues')

################################################################################################
#%%  To calculate the statistical significance (p-value) of an AUC (Area Under the ROC Curve) score I used a permutation test.
################################################################################################

### INHIBITOR PERMUTATION TEST
# Number of permutations
n_permutations = 10000
# Assume y_true is your true labels and y_score is your predicted scores
y_true = consistency_DF.MoA.replace({'activator': 0, 'inhibitor': 1, 'other': 0})
y_score = consistency_DF.MoA_AS.replace({'activator': 0, 'inhibitor': 1, 'other': 0})
# Calculate the observed AUC
observed_auc = roc_auc_score(y_true, y_score)
# Initialize a counter for permuted AUCs that are greater than or equal to the observed AUC
count = 0
auc_distributionINH = []
# Perform the permutation test
for _ in range(n_permutations):
    # Permute the true labels
    y_true_permuted = np.random.permutation(y_true)
    # Calculate the permuted AUC
    permuted_auc = roc_auc_score(y_true_permuted, y_score, multi_class='ovr')
    auc_distributionINH.append(permuted_auc)
    # If the permuted AUC is greater than or equal to the observed AUC, increment the counter
    if permuted_auc >= observed_auc:
        count += 1
# Calculate the p-value
# Addition of a small number to avoid division by zero - Laplace smoothing
p_value = (count + 0.00000001) / (n_permutations + 0.00000001)
print(f'p-value: {p_value}')

#### ACTIVATOR PERMUTATION TEST
# Assume y_true is your true labels and y_score is your predicted scores
y_true = consistency_DF.MoA.replace({'activator': 1, 'inhibitor': 0, 'other': 0})
y_score = consistency_DF.MoA_AS.replace({'activator': 1, 'inhibitor': 0, 'other': 0})
# Calculate the observed AUC
observed_auc = roc_auc_score(y_true, y_score)
# Initialize a counter for permuted AUCs that are greater than or equal to the observed AUC
count = 0
auc_distributionACT = []
# Perform the permutation test
for _ in range(n_permutations):
    # Permute the true labels
    y_true_permuted = np.random.permutation(y_true)
    # Calculate the permuted AUC
    permuted_auc = roc_auc_score(y_true_permuted, y_score)
    auc_distributionACT.append(permuted_auc)
    # If the permuted AUC is greater than or equal to the observed AUC, increment the counter
    if permuted_auc >= observed_auc:
        count += 1
# Calculate the p-value + Laplace smoothing
p_value = (count + 0.00000001) / (n_permutations + 0.00000001)
print(f'p-value: {p_value}')

#%%  plot the distribution of the AUCs in the same plot, 1 column 2 rows
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(5,5))
sns.histplot(ax=axes[0], data = pd.DataFrame(auc_distributionINH, columns=['inhibitor']), stat = 'probability' ) 
axes[0].axvline(0.93, color = 'red') # vertical line with inhibitor ROC curve value
sns.histplot(ax=axes[1], data = pd.DataFrame(auc_distributionACT, columns=["activator"]), stat = 'probability' ) 
axes[1].axvline(0.91, color = 'red') # vertical line with activator ROC curve value
fig.savefig('/container/Results/auroc_permutation_test.pdf', dpi = 600, bbox_inches='tight')
################################################################################################
#%% 
# ----------------------------------------------------------------------------------
# Roc curve for classifier absed on the mean of the allelic series estimates at increasing thresholds of confidence
# ----------------------------------------------------------------------------------

listOut = []
for cred in consistency_DF.confidence.unique():
    # select only confident ones in both DF
    tmp = consistency_DF[ consistency_DF.confidence >= cred]
    plt.clf()
    fig, ax = plt.subplots(figsize=(6, 6))
    classes_of_interest = ['activator','inhibitor']
    colors = ["darkorange", "cornflowerblue"]
    for class_of_interest, color in zip(classes_of_interest, colors):
        consistencyMoAdummy = pd.get_dummies(tmp.MoA)
        pp = tmp['PP_slope>0'] if class_of_interest == 'inhibitor' else tmp['PP_slope<0'] 
        display = RocCurveDisplay.from_predictions(
            consistencyMoAdummy[[class_of_interest]],
            pp,
            name=class_of_interest,
            color=color,
            ax=ax,
            plot_chance_level=(class_of_interest == 'inhibitor')
        )
        _ = ax.set(
        xlabel="False Positive Rate",
        ylabel="True Positive Rate"
        )
    plt.savefig(f'/container/ROC_curve_test_{cred}.png', dpi = 600)
    tmp = tmp[tmp['MoA'] != 'other']
    cm = confusion_matrix(tmp.MoA, tmp.MoA_AS, labels=tmp.MoA.unique() )
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=tmp.MoA.unique())
    disp.plot(cmap='Blues').figure_.savefig(f'/container/AS_MoA_confusion_matrix{cred}.png', dpi = 600)

    # Calculate the Chi-square test of independence
    chi2, p, dof, expected = chi2_contingency(cm)

    # Assume cm is your confusion matrix
    # cm[0, 0] = TN, cm[0, 1] = FP, cm[1, 0] = FN, cm[1, 1] = TP
    TN, FP, FN, TP = cm.ravel()
    # Calculate the odds ratio
    OR = ( (TP + 0.1) * (TN + 0.1) ) / ( (FP + 0.1) * (FN + 0.1) ) # to overcome eventual division by 0 (i.e. Laplace smoothing)
    dictOut = {'cred': [cred], 'n_drugs': tmp.shape[0],'chi2': chi2, 'p_value': p, 'dof': dof, 'OR': OR, 'TP': TP, 'TN': TN, 'FP': FP, 'FN': FN}
    listOut.append(pd.DataFrame(dictOut))

pd.concat(listOut).to_csv('/container/chi2_test_AS_MoA_directionality_stats.csv', index=False)


################################################################################################
#%% Predicted genes directionalities
################################################################################################
predicted_genes = pd.read_csv('/container/Resources/genelist_with_priorityScore_greater_than07.csv', skiprows=1, header=None)
df['predicted'] = df.gene.isin(predicted_genes[0])
# return only the lines with the max absolute of the mean
idx = df.groupby('gene')['median'].apply(lambda x: x.abs().idxmax())
max_abs_medians = df.loc[idx]
sns.histplot( data = max_abs_medians[max_abs_medians.predicted], x = 'median', hue = 'MoA_AS' )
#%% 
# drugInfoHI.groupby('MoA').count()
# df.groupby('MoA_AS').count()
df_barplot = pd.DataFrame( {'source': ['drug','drug','drug', 'gentic', 'gentic', 'gentic'],
                'type': ['activator','inhibitor','other', 'activator', 'inhibitor', 'other'],
                'n': [14181, 28019, 3318,4919,2936,0]} )

# calculate the proportion of n given the source
df_barplot = df_barplot.groupby('source')\
    .apply(lambda x: x.n / sum(x.n))\
        .reset_index(drop=False)\
            .rename(columns={'level_1': 'type'})

df_barplot.type = df_barplot.type.replace({0:'activator', 1:'inhibitor', 2:'other', 3:'activator', 4:'inhibitor', 5:'other' })

sns.barplot( data=df_barplot, y = 'n', x = 'type', hue='source')
plt.ylabel('Percentage')
