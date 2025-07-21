import pandas as pd
import sys

genes_file = sys.argv[1] 
pheno_conv_table = sys.argv[2] 
az_data = sys.argv[3] 
# Phenotypes
genes = pd.read_csv(genes_file, header=None)
# conv_table
conv_table = pd.read_csv(pheno_conv_table)
# AZ
az_ds = pd.read_parquet(az_data)

# limit to the genes of interest
# Remove the ' from the gene names - This is an AZ thing that is carried over from the previous pipeline
filter_list = az_ds.Gene \
                .str.replace("'","") \
                .isin( genes[0].to_list() )
az_ds = az_ds[ filter_list ]
# filter for significant pvalues
az_ds = az_ds[ az_ds['p-value'] < 0.00000005 ]

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
# drop rows where phenotypes is NaN
AZ_vars = AZ_vars.dropna(subset=['diseaseFromSourceMappedId'])
# some phenotype codes are divided by a pipe, so we need to split them and add them as new rows
AZ_vars['diseaseFromSourceMappedId'] = AZ_vars.diseaseFromSourceMappedId.str.split('|')
AZ_vars = AZ_vars.explode('diseaseFromSourceMappedId').reset_index(drop=True)

AZ_vars[['Variant', 'Varianttype', 'Category', 'Model', 'Consequencetype',
       'Gene', 'Aminoacidchange', 'Exonrank',
       'p-value', 'Oddsratio', 'OddsratioLCI', 'OddsratioUCI',
       'diseaseFromSourceMappedId']].to_csv("az_variants.csv")