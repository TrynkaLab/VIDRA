from google.oauth2 import service_account
import pandas_gbq 
import pandas as pd
import numpy as np
import sys

coloc_vars = sys.argv[1]

#Cluster of phenotypes
coloc_vars = pd.read_csv(coloc_vars)

studies = (
          coloc_vars.
            left_study.
            drop_duplicates().
            tolist()
          )

chromosomes = (
    coloc_vars.
      left_chrom.
      drop_duplicates().
      tolist()
)

position = (
    coloc_vars.
      left_pos.
      drop_duplicates().
      tolist()
)

#BigQuery 

# the -- is needed to comment SQL scripts
SQL_request = """
  SELECT 
    study_id, chrom, pos, ref, alt, beta, se, pval, n_total, n_cases, eaf, mac, mac_cases, info, is_cc
  FROM 
    `open-targets-genetics.genetics.sa_gwas`
  GROUP BY study_id, chrom, pos, ref, alt, beta, se, pval, n_total, n_cases, eaf, mac, mac_cases, info, is_cc
  HAVING chrom IN UNNEST(@chromosome)
    AND study_id IN UNNEST(@study)
    AND pos IN UNNEST(@position)
  """

# I do not use ref and alt because some may be flipped and this would cause the line to be filtered
query_conf= {
    'query': {
        'parameterMode': 'NAMED',
        'queryParameters': [
            # Fitler 1
            {
                'name': 'study',
                'parameterType': {'type': 'ARRAY',
                                  'arrayType': {'type': 'STRING'}},
                'parameterValue': {'arrayValues': [{'value': i} for i in studies]}
            },
            {
                'name': 'chromosome',
                'parameterType': {'type': 'ARRAY',
                                  'arrayType': {'type': 'STRING'}},
                'parameterValue': {'arrayValues': [{'value': i} for i in chromosomes]}
            },
            {
                'name': 'position',
                'parameterType': {'type': 'ARRAY',
                                  'arrayType': {'type': 'INTEGER'}},
                'parameterValue': {'arrayValues': [{'value': i} for i in position]}
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
df.to_csv('GWAS_variants.csv')