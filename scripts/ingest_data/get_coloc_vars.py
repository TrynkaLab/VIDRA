from google.oauth2 import service_account
import pandas_gbq 
import pandas as pd
from sys import argv
import numpy as np

df_unnested_pheno = argv[1] 
coloc_th = argv[2] # theshold for the colocalisation H4 filter
#Cluster of phenotypes
phenotypes = pd.read_csv(df_unnested_pheno)
pheno_list = (
pd.unique(
    np.array(
        phenotypes.filter(
            regex=r'^pheno', axis=1
            )
            )
            .flatten()
            )
            .tolist()
)
not_na_index = [ x is not None for x in pheno_list ] 

lis = np.array(pheno_list)
filter = np.array(not_na_index)
pheno_list_clean = lis[filter].tolist()

studies = phenotypes.study_id.drop_duplicates().tolist()

#BigQuery 
# the -- is needed to comment SQL scripts
SQL_request = """
  SELECT 
    *
  FROM 
    `open-targets-genetics.genetics.variant_disease_coloc`
  WHERE coloc_h4 > @ct
    AND ( right_study IN UNNEST(@stu) 
          OR left_study IN UNNEST(@stu) )
    AND ( right_bio_feature IS NOT NULL 
          OR right_phenotype IS NOT NULL 
          OR right_gene_id IS NOT NULL )
  """

query_conf= {
    'query': {
        'parameterMode': 'NAMED',
        'queryParameters': [
            # Fitler 1
            {
                'name': 'stu',
                'parameterType': {'type': 'ARRAY',
                                  'arrayType': {'type': 'STRING'}},
                'parameterValue': {'arrayValues': [{'value': i} for i in studies]}
            },
            # Fitler 2
            {
                'name': 'ct',
                'parameterType': {'type': 'FLOAT'},
                'parameterValue': {'value': coloc_th}
            }
        ]
    }
}

df = pandas_gbq.read_gbq(
    SQL_request,
    project_id="",
    credentials="",
    configuration=query_conf
    )
df.to_csv("colocalising_variants.csv")