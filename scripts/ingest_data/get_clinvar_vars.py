from google.oauth2 import service_account
import pandas_gbq 
import pandas as pd
import sys

qtl_info = sys.argv[1] 
qtls_genes = pd.read_csv(qtl_info)\
    .iloc[0:100000] # This is to test the script and make it quicker. Is hould go for the proper setup

#BigQuery 
# SQL request
SQL_request = """
  SELECT 
    clinicalSignificances.list,
    confidence,
    diseaseFromSourceId,
    diseaseFromSourceMappedId,
    studyId,
    targetFromSourceId,
    targetId,
    variantFunctionalConsequenceId,
    variantHgvsId,
    variantId,
    variantRsId,
    diseaseId,
    score,
    datasourceId,
  FROM 
    `open-targets-prod.platform.evidence`,
  UNNEST(clinicalSignificances.list) AS clincalsignificanceunnest
  WHERE datasourceId = "eva"
    AND targetId IN UNNEST(@genes)
    AND confidence != "criteria provided, conflicting interpretations"
    AND confidence != "no assertion criteria provided" 
    AND confidence != "no assertion provided" 
    AND variantId IS NOT NULL
  """

# This is to filter the SQL request
# It commented out because it is filtering fro pheontype
# Can be implemented back but is should filter by gene 

query_conf= {
    'query': {
        'parameterMode': 'NAMED',
        'queryParameters': [
            # Fitler 1
            {
                'name': 'genes',
                'parameterType': {'type': 'ARRAY',
                                  'arrayType': {'type': 'STRING'}},
                'parameterValue': {'arrayValues': [{'value': i} for i in qtls_genes.gene_id.dropna().unique().tolist() ]}
            }
        ]
    }
}


df = pandas_gbq.read_gbq(
    SQL_request,
    project_id='',
    credentials='',
    configuration=query_conf
    ).dropna(axis=1,how='all')

# # Flatten clincal significance
df["clinicalSignificances"] = df['list'].apply(lambda x: list(x.flatten()[0].values())[0] )

df.to_csv("clinvar_variants.csv")
