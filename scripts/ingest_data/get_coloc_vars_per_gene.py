from google.oauth2 import service_account
import pandas_gbq 
import pandas as pd
from sys import argv
import numpy as np

gene_file = argv[1] 
genes = pd.read_csv(gene_file, header=None)

coloc_th = argv[2] # theshold for the colocalisation H4 filter

#BigQuery 
# the -- is needed to comment SQL scripts
SQL_request = """
  SELECT 
    *
  FROM 
    `open-targets-genetics.genetics.variant_disease_coloc`
  WHERE coloc_h4 > @ct
    AND ( right_type != 'gwas' AND right_type != 'sqtl'  )
    AND right_gene_id IN UNNEST(@genes)
  """

query_conf= {
    'query': {
        'parameterMode': 'NAMED',
        'queryParameters': [
            {
                'name': 'genes',
                'parameterType': {'type': 'ARRAY',
                                  'arrayType': {'type': 'STRING'}},
                'parameterValue': {'arrayValues': [{'value': i} for i in genes[0].tolist()]}
            },
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