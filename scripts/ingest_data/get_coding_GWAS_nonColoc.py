from google.oauth2 import service_account
import pandas_gbq 
import pandas as pd
import numpy as np
import sys
########################################################################
#BigQuery 
#############################################################################
# read in info and decide what to do
gene_list = sys.argv[1]
# check the entry is present
if gene_list == '':
  # the -- is needed to comment SQL scripts
  SQL_request = """
    SELECT 
      *
    FROM 
      `open-targets-genetics.genetics.variants` 
    WHERE 
      gene_id_prot_coding IS NOT NULL
    AND most_severe_consequence IN ('inframe_insertion',
      'frameshift_variant',
      'stop_gained',
      'splice_donor_variant',
      'coding_sequence_variant',
      'stop_lost',
      'stop_retained_variant',
      'missense_variant',
      'incomplete_terminal_codon_variant',
      'protein_altering_variant',
      'start_lost',
      'synonymous_variant',
      'splice_acceptor_variant',
      'inframe_deletion')
    """
  # Get the list of coding variants presents in Open Targets Genetics
  df_req1 = pandas_gbq.read_gbq(
      SQL_request,
      project_id="",
      credentials=""
      )
else:
  # read the list of genes
  df_genes = pd.read_csv(gene_list, header=None)
  # the -- is needed to comment SQL scripts
  SQL_request = """
    SELECT 
      *
    FROM 
      `open-targets-genetics.genetics.variants` 
    WHERE 
      gene_id_prot_coding IN UNNEST(@genes)
    AND most_severe_consequence IN ('inframe_insertion',
      'frameshift_variant',
      'stop_gained',
      'splice_donor_variant',
      'coding_sequence_variant',
      'stop_lost',
      'stop_retained_variant',
      'missense_variant',
      'incomplete_terminal_codon_variant',
      'protein_altering_variant',
      'start_lost',
      'splice_acceptor_variant',
      'inframe_deletion')
    """
  # Get the list of coding variants presents in Open Targets Genetics
  df_req1 = pandas_gbq.read_gbq(
      SQL_request,
      project_id="",
      credentials="",
      configuration={'query': {'parameterMode': 'NAMED','queryParameters': [{'name': 'genes','parameterType': {'type': 'ARRAY','arrayType': {'type': 'STRING'}},'parameterValue': {'arrayValues': [{'value': i} for i in df_genes[0].tolist()]}}]}}
      )

########################################################################
# We got the variants, let's get the GWAS information
print(f"Number of coding variants in OT: {df_req1.shape[0]}")

# Use this variants to extract the GWAS variant
# Step 1: get the variants filtering by position
positionsFilter = df_req1.position.drop_duplicates(inplace=False).tolist()
length = len(positionsFilter)
part1 = positionsFilter[:length//3]
part2 = positionsFilter[length//3:2*length//3]
part3 = positionsFilter[2*length//3:]

SQL_request = """
  SELECT 
    * 
  FROM 
    `open-targets-genetics.genetics.sa_gwas` 
  WHERE 
    pos IN UNNEST(@position)
  AND
    pval <= 0.00000005
"""

#  I need to split the query in 3 parts because the query is too long
dfOut = []
for part in [part1, part2, part3]:
  query_conf= {
      'query': {
          'parameterMode': 'NAMED',
          'queryParameters': [
              # Fitler 1
              {
                  'name': 'position',
                  'parameterType': {'type': 'ARRAY',
                                    'arrayType': {'type': 'INTEGER'}},
                  'parameterValue': {'arrayValues': [{'value': i} for i in part]}
              }
          ]
      }
  }
  dfTmp = pandas_gbq.read_gbq(
      SQL_request,
      project_id="",
      credentials="",
      configuration=query_conf
      )
  dfOut.append(dfTmp)

df_req2 = pd.concat(dfOut, axis=0, ignore_index=True)
df_req2.rename(columns={'chrom': 'chr_id', 'pos':'position', 'ref': 'ref_allele', 'alt': 'alt_allele'}, inplace=True)
# Step 2: get the variants match usgin left join to the original dataframe
# join on position and chromosome, also ref and alt alleles
df_req1 = df_req1.set_index(['chr_id', 'position', 'ref_allele', 'alt_allele'])
df_req2 = df_req2.set_index(['chr_id', 'position', 'ref_allele', 'alt_allele'])
df_vars = df_req1.join(
                    df_req2, 
                    on=['chr_id', 'position', 'ref_allele', 'alt_allele'], 
                    how='left',
                    validate='many_to_many') 
df_vars = df_vars[~df_vars['study_id'].isna()]         

df_req1.to_csv('GWAS_coding_variants_fromOT.csv')
df_vars.to_csv('GWAS_coding_variants_NoColocalising.csv')