import sys
import pandas as pd
import re

cluster_file = sys.argv[1] 
subset_cluster = sys.argv[2] 
subset_cluster = pd.to_numeric(subset_cluster.strip('][').split(","))
#Cluster of phenotypes
cluster_pheno = pd.read_csv(cluster_file)

if pd.isnull(subset_cluster).any():
    ids = (
        cluster_pheno
        .ID
        .str.split("-", n=2, expand=True)[1]
        .dropna()
        .drop_duplicates()
        .reset_index() 
    )
else:
    ids = (
        cluster_pheno[cluster_pheno['cl_08'].isin(subset_cluster)]
        .ID
        .str.split("-", n=2, expand=True)[1]
        .dropna()
        .drop_duplicates()
        .reset_index()
    )

ids[1].to_csv("phenotypes.csv", index=False)
