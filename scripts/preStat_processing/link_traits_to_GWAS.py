from google.oauth2 import service_account
import pandas_gbq 
import pandas as pd
from sys import argv

#BigQuery 
## BigQuery credentials
# n_variants = 0

# the -- is needed to comment SQL scripts
SQL_request = """
  SELECT 
    study_id,
    n_cases,
    n_initial,
    n_replication,
    num_assoc_loci,
    trait_reported,
    trait_efos AS trait_efos
  FROM 
    `open-targets-genetics.genetics.studies`, 
  UNNEST(trait_efos.list) AS tr
  WHERE num_assoc_loci > @number_variants
  """

query_conf= {
    'query': {
        'parameterMode': 'NAMED',
        'queryParameters': [
            # Fitler 1
            {
                'name': 'number_variants',
                'parameterType': {'type': 'INTEGER'},
                'parameterValue': {'value': n_variants}
            }
        ]
    }
}

df = pandas_gbq.read_gbq(
    SQL_request,
    project_id="",
    credentials="",
    configuration=query_conf,
    index_col="study_id"
    )

df.reset_index(inplace=True)

df2 = pd.json_normalize(df["trait_efos"]).list.apply(pd.Series)
df2 = df2.applymap(lambda x: x['element'] if isinstance(x, dict) else np.nan )

df3 = pd.concat([df, df2], axis=1)
df3.drop(columns=["n_cases", "n_initial","n_replication", "num_assoc_loci", "trait_reported", "trait_efos"], inplace=True)

df_long = df3.melt(id_vars='study_id', var_name='groups', value_name='phenotypes')
df_long = df_long[~df_long.phenotypes.isna()]
df_long.reset_index(drop=True, inplace=True)

df.to_csv("gwas_study_info.csv.csv")
df_long.to_csv("trait_to_gwas.csv")