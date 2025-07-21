from google.oauth2 import service_account
import pandas_gbq 
import pandas as pd
import numpy as np
import sys

# tmp = '/container/work/80/9b7c9e0ec4609599b0f784242bb94b/burden_tests.csv'
# df_unnested_pheno = tmp
genes_file = sys.argv[1] 
# Phenotypes
genes = pd.read_csv(genes_file)

#BigQuery 
## BigQuery credentials
## Note: these are private, so don't have to figure in pucliv repositories

# SQL request
SQL_request = """
  SELECT 
    datasourceId, 
    targetId, 
    allelicRequirements, 
    ancestryId, 
    beta, 
    betaConfidenceIntervalLower, 
    betaConfidenceIntervalUpper, 
    cohortId, 
    diseaseFromSource, 
    diseaseFromSourceMappedId, 
    oddsRatio, 
    oddsRatioConfidenceIntervalLower, 
    oddsRatioConfidenceIntervalUpper, 
    pValueExponent, 
    pValueMantissa, 
    projectId, 
    resourceScore, 
    statisticalMethod, 
    targetFromSourceId, 
    score,
    studyCases,
    studyCasesWithQualifyingVariants,
    studySampleSize
  FROM 
    `open-targets-prod.platform.evidence`
  WHERE datasourceId = "gene_burden"
  """

df = pandas_gbq.read_gbq(
    SQL_request,
    project_id="",
    credentials=""
    ).dropna(axis=1,how='all')

# conver beta to OR
df['oddsRatio'] = df["oddsRatio"].fillna(np.exp(df["beta"]))
df['oddsRatioConfidenceIntervalLower'] = df["oddsRatioConfidenceIntervalLower"]\
  .fillna(
      np.exp(df["beta"]) - 1.99 * 0.5346 
      )
df['oddsRatioConfidenceIntervalUpper'] = df["oddsRatioConfidenceIntervalUpper"]\
  .fillna(
      np.exp(df["beta"]) - 1.99 * 0.5346 
      )
# remove rows where still get nan
# and remove beta measures  - everythong was incorporated in the OR columns
df = df.loc[~df.oddsRatio.isna(),:]\
  .drop(['Unnamed: 0','beta', 
    'betaConfidenceIntervalLower', 
    'betaConfidenceIntervalUpper'], 
    axis=1,
    errors='ignore')

# # Flatten clincal significance
df.to_csv("burden_tests.csv")
