#%% 
import pandas as pd
import sys
import numpy as np
from math import floor
from sklearn.preprocessing import MinMaxScaler
sys.path.append('/container/scripts/modules_py')
import moduleAS as mAS
import os
import scipy.stats as ss
#%% 
####################################################################################
# read all the files into the memory
gene_name_conversion = pd.read_csv(sys.argv[13]) # Table to convert AZ variants gene symbol to ensembl ID
coloca_variants = pd.read_csv(sys.argv[1]) # Coloc - # This table is used to link GWAS and QTL
GWAS_variants = pd.read_csv(sys.argv[2]) # Read GWAS - information on phenotype effect of common variants
QTL_variants = pd.read_csv(sys.argv[3]) # Read QTLs - proxy for the variant effect on the protein acitivy
pheno_gwas_studies = pd.read_csv(sys.argv[4]) # This is a conversion table between study ID in the GWAS and phenotype (EFO)
AZvariants = pd.read_csv(sys.argv[5], low_memory=False) # AZvariants - Rare variants from UKB
RareVars_AZ_protein_info = pd.read_csv(sys.argv[6]) # RareVars_AZ effec on protein_info
clinvar_variants = pd.read_csv(sys.argv[7]) # clinvar_variants - Rare vriant and rare diseases
RareVars_ClinVar_Protein_info = pd.read_csv(sys.argv[8]) # RareVars_ClinVar effec on protein_info
burdens_AZ = pd.read_csv(sys.argv[10]) # Burden effet 
cdGWAS = pd.read_csv(sys.argv[11]) # coding missense variants identified in GWAS
cdGWASprotein = pd.read_csv(sys.argv[12]) # predivtion via VEP of the effect of GWAS missense (where QTLs are not available)
OlfEFOtoNewTable = pd.read_csv(sys.argv[14]) # This is a conversion table between old EFO IDs and the most updated ones
dir_out_name = sys.argv[9]

#%% 
####################################################################################
# convert the EFO table to dictionary
dict_EFOreplace = pd.Series(OlfEFOtoNewTable.id.values,index=OlfEFOtoNewTable.old_ids).to_dict()

# Remap all the EFO in the dataframes to the most updated version
# Replace the values in the DataFrame column
AZvariants['diseaseFromSourceMappedId'] = AZvariants['diseaseFromSourceMappedId'].replace(dict_EFOreplace)
burdens_AZ['diseaseFromSourceMappedId'] = burdens_AZ['diseaseFromSourceMappedId'].replace(dict_EFOreplace)
clinvar_variants['diseaseFromSourceMappedId'] = clinvar_variants['diseaseFromSourceMappedId'].replace(dict_EFOreplace)
pheno_gwas_studies['phenotypes'] = pheno_gwas_studies['phenotypes'].replace(dict_EFOreplace)

#%% 
####################################################################################
# Protein function info from VEP
# These are the columns that are kept from the VEP annotation
cols_keep = ['variantId','gene_id',
             'primateai','revel','cadd_phred','cadd_raw','pli_gene_value',
             'impact','gene_symbol','conservation','loftool','canonical',
             'consequence_terms','polyphen_prediction','am_pathogenicity','sift_prediction',
             'sift_score','am_class','blosum62','polyphen_score']

####################################################################################
# Work on the common GWAS QTL variants
####################################################################################
# These are scaffold on the coloc ones
# slim coloc file to retain only the info that are required to assemble GWAS and QTL
coloca_variants  = coloca_variants[['left_type', 'left_study', 'left_chrom','left_pos', 'left_ref', 'left_alt', 
                                    'right_type', 'right_study','right_bio_feature', 'right_phenotype', 'right_chrom', 'right_pos',
                                    'right_ref', 'right_alt', 'is_flipped', 'right_gene_id']]
# slim GWAS spredsheet
GWAS_variants = GWAS_variants[['study_id','chrom','pos','ref', 'alt', 'beta','se', 'pval']]
# slim QTL spreadsheet
QTL_variants = QTL_variants[['type_id', 'study_id', 'bio_feature', 'phenotype_id',
                               'gene_id', 'chrom', 'pos', 'ref', 'alt', 'beta', 'se', 'pval','is_cc']]

# slim pheno_gwas_studies
# Keep only the study ID and the phenotype translation to EFO
phenoTab = pheno_gwas_studies[['study_id','phenotypes']].drop_duplicates()
phenoTab = phenoTab.rename(columns={"phenotypes": "diseaseFromSourceMappedId"}) # Most tables use diseaseFromSourceMappedId

# Link GWAS to colocalising variants
# Assembling on coloc variants
gwas = pd.merge(
                coloca_variants,
                GWAS_variants,
                left_on=['left_study', 'left_chrom','left_pos', 'left_ref', 'left_alt'], 
                right_on=['study_id','chrom','pos','ref', 'alt'],
                how="inner"
                )

# Add suffix '_qtl' to the columns in QTL_variants
QTL_variants.columns = [str(col) + '_qtl' for col in QTL_variants.columns]

# Link QTL effects to colocalising variants
common_var = pd.merge(
                gwas,
                QTL_variants,
                left_on=['right_type','right_study','right_bio_feature','right_phenotype',
                         'right_gene_id',
                         'right_chrom','right_pos','right_ref','right_alt'], 
                right_on=['type_id_qtl','study_id_qtl','bio_feature_qtl','phenotype_id_qtl',
                          'gene_id_qtl',
                          'chrom_qtl', 'pos_qtl', 'ref_qtl', 'alt_qtl'],
                how="inner"
                )

# Attache phenotype to the GWAS
common_var = pd.merge(
    common_var,
    phenoTab,
    left_on="left_study",
    right_on="study_id",
    )

##################################
# Retain only the variants that have association in the selected tissues
cells_to_keepBLOOD = ['UBERON_0000178','BLOOD', 'UBERON_0001969', 'TREG_NAIVE', 'TREG_MEMORY', 'TH2_MEMORY', 'TH1_MEMORY', 'TH17_MEMORY',
                 'TH1-17_MEMORY','TFH_MEMORY','NK-CELL_NAIVE','MONOCYTE_NAIVE','MONOCYTE_CD16_NAIVE','CD8_T-CELL_NAIVE','CD8_T-CELL_ANTI-CD3-CD28',
                 'CD4_T-CELL_NAIVE','CD4_T-CELL_ANTI-CD3-CD28', 'B-CELL_NAIVE', 'MONOCYTE_R848', 'MONOCYTE_PAM3CSK4', 'MONOCYTE_LPS',
                 'MONOCYTE_IAV','MACROPHAGE_SALMONELLA','MACROPHAGE_NAIVE','MACROPHAGE_LISTERIA','NEUTROPHIL_CD16','T-CELL_CD8','T-CELL_CD4','NEUTROPHIL',
                 'MONOCYTE','MACROPHAGE_IFNG','MACROPHAGE_IFNG+SALMONELLA',
                 'T-CELL','MONOCYTE_LPS24','MONOCYTE_LPS2','MONOCYTE_IFN24','B-CELL_CD19','PLATELET','NEUTROPHIL_CD15','MONOCYTE_CD14']

common_var = common_var[ common_var.right_bio_feature.isin(cells_to_keepBLOOD) ]

# Keep only coloc where left is GWAS
common_var = common_var[common_var['left_type'] == 'gwas' ]

# create variant IDs
common_var['variantId'] = common_var \
                                .left_chrom.astype('int').astype("str") + \
                                "_" + common_var.left_pos.astype('int').astype("str") + \
                                "_" + common_var.left_ref + \
                                "_" + common_var.left_alt

# Relevant columns for GWAS and QTL
keep_gwas_qtl = [
                # GWAS
                'variantId', 'left_type', 'beta', 'se', 'pval', 'left_type', 'diseaseFromSourceMappedId',
                # QTL
                'right_type', 'right_study', 'right_bio_feature', 'right_phenotype', 'right_gene_id', 'beta_qtl', 
                'se_qtl', 'pval_qtl'
                ]

common_var = common_var[keep_gwas_qtl]
# group by variant, gene, disease and see how many tissues are associated to the variant
# if more than one of the remaining tissue/study, then take only the one with hte lowest GWAS and qtl pvalue

# Define the function to return the lowest pvalue row
def best_qtl(grouped_df):
    condition = grouped_df['pval_qtl'] == np.min(grouped_df['pval_qtl'])
    return grouped_df[condition] 

def best_gwas(grouped_df):
    condition = grouped_df['pval'] == np.min(grouped_df['pval'])
    return grouped_df[condition] 

common_var = common_var.groupby(['variantId','right_gene_id','diseaseFromSourceMappedId']) \
                    .apply(best_gwas) \
                    .reset_index(drop=True) \
                    .groupby(['variantId','right_gene_id','diseaseFromSourceMappedId']) \
                    .apply(best_qtl) \
                    .reset_index(drop=True)\
                    .drop_duplicates()\
                    .reset_index(drop=True)

# Add cols that will be helpful for the AS
common_var['var_type'] = "common_NC"

#%% 
####################################################################################
## AZ VARIANTS
####################################################################################
# Add Ensemble ID - astra zeneca returns them with gene symbol, i prefer to use ENSEMBL IDs
AZvariants['Gene'] = AZvariants.Gene.str.replace("\'","") # Gene symbols are quoted - remove '
AZvariants = AZvariants.merge(gene_name_conversion, left_on='Gene', right_on='symbol', how='left')

AZvariants = AZvariants[['Variant','gene_id','Model',"p-value", "Oddsratio", "OddsratioLCI", "OddsratioUCI", "diseaseFromSourceMappedId"]] \
    .drop_duplicates() \
    .dropna(subset=['diseaseFromSourceMappedId']) \
    .reset_index()

# There are multiple sep in the variants - this is to harmonise the ids
AZvariants['Variant'] = AZvariants.Variant.str.replace("-"," ")
AZvariants['Variant'] = AZvariants.Variant.str.replace(" ","_")

# Select only the allelic mode of inheritance to be consisted with GWAS
model_to_select = ['allelic']

# Remove the others MOI from the df
AZvariants = AZvariants[AZvariants.Model.isin(model_to_select)]

# Keep only decided protein info
RareVars_AZ_protein_info = RareVars_AZ_protein_info.iloc[ :, RareVars_AZ_protein_info.columns.isin(cols_keep)]

# The IDs in rare variants ends with _. - remove it
RareVars_AZ_protein_info['variantId'] = RareVars_AZ_protein_info.variantId.str.replace("_.","")

# Add info on variant effects on protein
AZ_df = pd.merge(
   AZvariants,
   RareVars_AZ_protein_info, 
   left_on='Variant',
   right_on="variantId",
   how='left' # left is to keep all the variants - some may be filtered in protein annotation and want to cehck those 
)

# convert Odds ratio to betas
AZ_df['beta_gwas'] = np.log(AZ_df.Oddsratio)
AZ_df['se_gwas'] = np.sqrt(50) * ( (AZ_df.OddsratioUCI - AZ_df.OddsratioLCI)/3.92 )

# Clean the AZ DF                    
AZ_df.drop(['Model', 'Oddsratio', 'OddsratioLCI', 'OddsratioUCI'],
                axis=1,
                inplace=True,
                errors='ignore')

# Attach AZ label
AZ_df['var_type'] = "rare_AZ"

# change name of variants id cols
AZ_df = AZ_df.rename(columns={"variantId":"Variant", \
                              "gene_id_x":"targetId", \
                              "gene_id_y":"targetId", \
                              "diseaseFromSourceMappedId": "diseaseId"})

#%% 
####################################################################################
# CLINVAR VARIANTS
####################################################################################
##clean cols
clinvar_variants = clinvar_variants[['variantId', 'targetId', 'diseaseId', 'score', 'clinicalSignificances']]

# Keep only decided protein info
RareVars_ClinVar_Protein_info = RareVars_ClinVar_Protein_info.iloc[ :, RareVars_ClinVar_Protein_info.columns.isin(cols_keep)]

# Combine rare ClinVar variants
clinvar_df = pd.merge(
   clinvar_variants,
   RareVars_ClinVar_Protein_info, 
   on='variantId',
   how='left'
)

# Remove unwanted columns 
clinvar_df.drop(['score','gene_id','gene_symbol'],
                axis=1,
                inplace=True,
                errors='ignore')

# Add info on source of the information
clinvar_df['var_type'] = "rare_CV"

# change name of variants id cols
clinvar_df = clinvar_df.rename(columns={"diseaseFromSourceMappedId": "diseaseId"}) 

#%% 
####################################################################################
# Coding missense variants identified in GWAS
####################################################################################
# Create variant IDs
cdGWAS['variantId'] = cdGWAS.chr_id.astype('str') + \
    "_"+cdGWAS.position.astype('str') + \
    "_"+cdGWAS.ref_allele + \
    "_"+cdGWAS.alt_allele

# Add the phenotype to the GWAS
cdGWAS = pd.merge(
    cdGWAS,
    phenoTab, 
    on="study_id",
    validate="m:m"
    )

# Keep the protine infomration from VEP decided at the beginning in cols keep variable
cdGWASprotein = cdGWASprotein.iloc[ :, cdGWASprotein.columns.isin(cols_keep)]

# Format attach gwas to protein info
cdGWAS_df = pd.merge(
   cdGWAS,
   cdGWASprotein, 
   on='variantId',
   how = 'left'
)

# Type ID is used to defint QTL type, sho should remove from this df
cdGWAS_df.drop(['type_id', 'position_b37', 'rs_id', 'chr_id', \
                'position', 'alt_allele', 'ref_allele', \
                'chr_id_b37', 'gene_id_any', 'info', \
                'gene_id', 'exon', 'inputStr', 'disgenet', \
                'lof', 'lof_info', 'nmd', 'lof_flags', \
                'rf_score', 'ada_score', 'study_id', \
                'gene_id_prot_coding_distance', 'gene_symbol', \
                'n_total', 'n_cases', 'eaf', 'mac', 'mac_cases', 'is_cc', \
                'gene_id_any_distance', 'cadd', 'af' ],
                axis=1,
                inplace=True,
                errors='ignore')

# group by variant, gene, disease and see how many studies are associated to the variant
# if multiple, retain only the one with the lowest pvalue
cdGWAS_df = cdGWAS_df.groupby(['variantId','gene_id_prot_coding','diseaseFromSourceMappedId']) \
                    .apply(best_gwas) \
                    .reset_index(drop=True)\
                    .drop_duplicates()

# rename gene_id column and diseaseFromSourceMappedId so that theyare consistent with the other df
cdGWAS_df = cdGWAS_df \
    .rename(
        columns={"gene_id_prot_coding": "gene_id", 
                 "diseaseFromSourceMappedId": "diseaseId",
                 "beta": "beta_gwas",
                 "se": "se_gwas",
                 "pval": "pval_gwas"})

cdGWAS_df['var_type'] = "common_CD"

#%% 
####################################################################################
# BURDEN TEST
####################################################################################

# change name of variants id cols
statistical_methods_to_keep = ['pLoF','ptv','ptvraredmg']
cohort_to_kepp = ['UK Biobank 450k']
project_to_keep = ['AstraZeneca PheWAS Portal']
burdens_AZ = burdens_AZ.rename(columns={"targetId": "gene_id"}) 

# constrain studies to wanted ones
# and wanted statistical methods
burdens_AZ = burdens_AZ[burdens_AZ.cohortId.isin(cohort_to_kepp) & \
                        burdens_AZ.projectId.isin(project_to_keep) & \
                        burdens_AZ.statisticalMethod.isin(statistical_methods_to_keep) ]

# Function to select the lowest pvalue once grouped by gene and disease
def best_burden(grouped_df):
    condition = grouped_df['pValueExponent'].min() == grouped_df['pValueExponent']
    return grouped_df[condition][0:1] # This is in the case of multiple variants haveing the same p-values. It returns the first one

# If there are multiple burden for same gene-disease pair, keep only the one with the lowest pvalue
# Since I decided to onlyuse one source of burden, this shouldn't make any difference - unless the same gene is identified in ptv and prvraredmg
burdens_AZ = burdens_AZ \
    .groupby(['gene_id', 'diseaseFromSourceMappedId']) \
    .apply(best_burden) \
    .reset_index(drop=True) 

burdens_AZ['sd'] = np.sqrt(burdens_AZ.studyCasesWithQualifyingVariants) * ( (burdens_AZ.oddsRatioConfidenceIntervalUpper - burdens_AZ.oddsRatioConfidenceIntervalLower)/3.92 )
burdens_AZ['beta'] = np.log(burdens_AZ.oddsRatio)

# remove unwanted columns
burdens_AZ.drop(['Unnamed: 0', 'cohortId', 'score',\
                 'allelicRequirements', 'datasourceId', \
                 'oddsRatio', \
                 'projectId', 'statisticalMethod', 'resourceScore', \
                 'pValueMantissa', 'oddsRatioConfidenceIntervalLower', \
                 'oddsRatioConfidenceIntervalUpper', \
                 'pValueExponent', 'studyCases', \
                 'studyCasesWithQualifyingVariants', \
                 'studySampleSize', \
                 'diseaseFromSource', 'targetFromSourceId', 'ancestryId'],
                 axis=1,
                 inplace=True, 
                 errors='ignore')

# Attach burden to every column name
burdens_AZ.columns = [col + '_burden' for col in burdens_AZ.columns]

# rename gene_id column and diseaseFromSourceMappedId so that theyare consistent with the other df
burdens_AZ = burdens_AZ.rename(columns={"gene_id_burden": "gene_id", "diseaseFromSourceMappedId_burden": "diseaseId"})

#%% 
####################################################################################
# Assemble all variants in same table 
####################################################################################
# Assemble all the dataframes in one  - crate a list of dataframes to concatenate
df_to_concat = [common_var, AZ_df, clinvar_df, cdGWAS_df]

# Rename the relavant columns to make them consistent across all df
for i in range(len(df_to_concat)):
    df_to_concat[i] = df_to_concat[i].rename(columns= 
                            {'Variant': 'variantId', 
                            'targetId': 'gene_id', 
                            'diseaseFromSourceMappedId': 'diseaseId', 
                            'right_gene_id': 'gene_id',
                            'beta': 'beta_gwas',
                            'se': 'se_gwas',
                            'se_qtl': 'se_qtl',
                            })

# Concatenate the dataframes along the rows
finDF = pd.concat(df_to_concat, axis=0)

# Attach burden measures
finDF = pd.merge(
                finDF,
                burdens_AZ,
                on = ['diseaseId','gene_id'],
                how = 'left'
                )

# Drop AS not relevant columns
finDF.drop(['left_type', 'pval', 'pval_gwas','right_study',  'right_bio_feature', 
            'right_phenotype', 'pval_qtl', 'gene_symbol',
            'index','p-value'],
            errors='ignore', axis=1, inplace=True)

####################################################################################
# NORMALISING THE DATA
####################################################################################
# normalising functions 

scaler = MinMaxScaler()
# Fill in na per with values that comes from other predictors in the same row.

# Prediction on protein activity
# REVEL has the right distribution. it only needs to be inverted (i.e. 1 - revel)
finDF['revel'] = 1 - finDF['revel'] # In revel pathogenic variants are close to 1. I need them to be close to 0

# CADD goes from -4 to 12, with highest being more pathongenic. 
# scale between 0 and 1 and invert as 1-cadd
finDF['cadd_phred'] = 1 - scaler.fit_transform(finDF['cadd_phred'].to_numpy().reshape(-1, 1)).flatten()

# conservation is ok, it only needs scaling between 0 and 1 and flipped
finDF['conservation'] = 1 - scaler.fit_transform(finDF['conservation'].to_numpy().reshape(-1, 1)).flatten()

# polyphen_score needs to be inverted
finDF['polyphen_score'] = 1 - finDF['polyphen_score'] 

####################################################################################
# CALCULATE THE STD DEVIATION TO USE IN BAYESIAN MODELLING
####################################################################################
# Iterate across the columns 
for col in finDF.columns:
    # if the column is not numeric, skip it
    if not np.issubdtype(finDF[col].dtype, np.number):
        continue
    # calculate the standard deviation of the column
    print( col, np.std(finDF[col]) )
    # This is output in the stdout and then pasted in the bayesian model

#%% 
####################################################################################
# FILL IN THE MISSING VALUES WITH MEAN OF THE COLUMN
####################################################################################
colz = ['revel', 'cadd_phred', 'conservation', 'sift_score', 'polyphen_score']
# The fill is it is better to be done per gene 
# - so that a gene that has higher values doesn't skew the mean
# and at the same time it consider the implicit relevance of the gene
finDFg = finDF[colz + ['gene_id']].groupby('gene_id', group_keys=False)
finDFg_colz = finDFg\
                    .apply(lambda x: x.fillna(value=np.mean(x)) )\
                    .reset_index(drop=True)

# Attach the filled columns to the main dataframe
finDF = pd.concat([finDF.drop(colz, axis=1), finDFg_colz], axis=1)

# Manually hardcode 0 to the premature stop codons or start lost
finDF.loc[ finDF.most_severe_consequence.isin(['start_lost','stop_gained']), colz ] = 0

# Drop the most severe consequence column
finDF.drop('most_severe_consequence', axis=1, inplace=True)

#%% 
####################################################################################
# INTERPRET "NUMERICALLY" clinical significance 
####################################################################################
# In this case clinical significance and primate AI will need to be: pathogenic variants close to 1 and benign close to 0
# In these scores the higher the score the more pathogenic the variant is. Which is ok for the model.
order_list_CV = ['not provided','association not found', 'other', 'benign', 'likely benign',
              'low penetrance', 'confers sensitivity', 'uncertain risk allele', 
              'drug response',
              'uncertain significance', 'association', 'affects', 
              'likely risk allele', 'risk factor',
              'established risk allele', 'likely pathogenic', 'pathogenic']

finDF['clinicalSignificances'] = finDF.clinicalSignificances.replace(np.nan, 'not provided')
finDF['clinicalSignificances'] = pd.Categorical(finDF['clinicalSignificances'], categories=order_list_CV, ordered=True)
finDF['clinicalSignificances'] = finDF['clinicalSignificances'].cat.codes
finDF['clinicalSignificances'] = scaler.fit_transform(finDF['clinicalSignificances'].to_numpy().reshape(-1, 1)).flatten()

# Impute missing values with the mean of the column
colz2 = ['primateai', 'clinicalSignificances']
finDF[colz2] = finDF[colz2]\
                    .apply(
                    lambda x: x\
                        .fillna(value=pd.to_numeric(x,errors='coerce')\
                        .dropna()\
                        .mean()), 
                    axis=1)

# remove rows with missing values on disease or gene
finDF = finDF[~finDF.gene_id.isna()]
finDF = finDF[~finDF.diseaseId.isna()]

# harmonise the way variants are named
finDF.variantId = finDF.variantId.str.replace(" ","_")

####################################################################################
# HARD CODE IN SOME VALUE MANUALLY AND STRUCTURE THE OUTPUT
####################################################################################
values = {"am_pathogenicity": 1, "primateai": 0, "revel": 1, "cadd_phred": 1, 
          "clinicalSignificances": 0, "conservation": 1, 
          "sift_score": 1, "polyphen_score": 1,
          "blosum62": 1, "loftool": 1, 
          'sd_burden': 2, 'beta_burden': 0, 
          'se_gwas': .14, 'beta_gwas': 0,
          'beta_qtl': 0, 'se_qtl': .1}

finDF.fillna(value=values, inplace=True)

# Rename the columns to be consistent with the model nomenclature
finDF.rename(columns={'beta_gwas': 'yc',
                        'se_gwas': 'ycse',
                        'beta_qtl': 'xc',
                        'se_qtl': 'xcse',
                        'beta_burden': 'bO',
                        'sd_burden': 'bOse',
                        'blosum62': 'as_blosum62',
                        'as_foldx' 
                        'as_plddt' 
                        'loftool': 'as_loftool',
                        'conservation': 'as_conservation',
                        'sift_score': 'as_sift',
                        'polyphen_score': 'as_polyphen',
                        'most_severe_consequence': 'as_consequence',
                        'clinicalSignificances': 'as_clinicalSignificance',
                        'cadd_phred': 'as_cadd',
                        'revel':'as_revel',
                        'primateai': 'as_primateai',
                        'am_pathogenicity': 'as_alphamissense',
                        'diseaseId': 'as_disease',
                        'gene_id': 'as_gene',
                        'variantId': 'variant'
                }, inplace=True)

#%% 
# Source
order_list_source = ['common_NC', 'rare_AZ', 'rare_CV', 'common_CD']
finDF['GsourceLab'] = pd.Categorical(finDF.var_type, categories=order_list_source, ordered=True)
finDF['GsourceLab'] = finDF['GsourceLab'].cat.codes

# QTLs
order_list_qtl = ['eqtl', 'pqtl','no_qtl']
finDF.right_type.fillna('no_qtl', inplace=True)
finDF['GqtlLab'] = pd.Categorical(finDF.right_type, categories=order_list_qtl, ordered=True)
finDF['GqtlLab'] = finDF['GqtlLab'].cat.codes
finDF.drop(['var_type','right_type'], axis=1, inplace=True)

#%% 
####################################################################################
# OUTPUT THE DATAFRAME PARTITIONS
####################################################################################
# If scaling up to all genes and diseases this creates too many partitions and causes memory issues
# Therefore it is more conveninet to partition by gene only - this is because the model is gene centric
# This reflect in the bayes script too

# The output dir shuold be /container/pre_AS/output_all_genes_20240320/bayes/preASwrangle_bayesDF/
paritions = len(finDF.as_gene.unique())

# sort by partition column
finDF = finDF.sort_values('as_gene')

# output the dataframes
finDF.to_parquet(
            path = dir_out_name + "_bayesDF",
            compression = 'snappy',
            partition_cols = 'as_gene',
            engine='fastparquet', # if use pyarrow it will have problem with some partitions
            max_partitions = paritions
        )