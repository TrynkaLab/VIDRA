import pandas as pd
import sys

df_unnested_pheno = sys.argv[1] 
pheno_conv_table = sys.argv[2] 
az_data = sys.argv[3] 
# Phenotypes
phenotypes = pd.read_csv(df_unnested_pheno)
# conv_table
conv_table = pd.read_csv(pheno_conv_table)
# AZ
az_ds = pd.read_parquet(az_data)

# conv_table_link = ( conv_table[["diseaseFromSource","diseaseFromSourceMappedId"]] )
conv_table_link = ( conv_table[["Study tag", "Mapped trait"]] )
# This limits the AZ variants to the phenotypes we are interested in
# I commented this out because I think it is interesting to see all the phenotypes
  # [ lambda x: x.diseaseFromSourceMappedId.isin( phenotypes['1'].to_list() ) ]
  # )

az_ds = az_ds.set_index( 'Phenotype' )
conv_table_link = conv_table_link.set_index('Study tag')
AZ_vars = pd.merge(
              az_ds, 
              conv_table_link, 
              left_index=True, 
              right_index=True, 
              how='left', 
              validate="m:m")
# change column name to beck fit previous nomenclature.
AZ_vars['diseaseFromSourceMappedId'] = AZ_vars['Mapped trait']

AZ_vars[['Variant', 'Varianttype', 'Category', 'Model', 'Consequencetype',
       'Gene', 'Aminoacidchange', 'Exonrank',
       'p-value', 'Oddsratio', 'OddsratioLCI', 'OddsratioUCI',
       'diseaseFromSourceMappedId']].to_csv("az_variants.csv")