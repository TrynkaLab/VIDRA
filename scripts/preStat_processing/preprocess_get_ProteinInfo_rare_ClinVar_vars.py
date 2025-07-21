import pandas as pd
import sys
import numpy as np
from math import ceil

# Output name and path. It will be a parquet file with the ProtVar info.
input_name, col_vars, separator = sys.argv[1], sys.argv[2], sys.argv[3]
# Read in rare variants
vars = pd.read_csv(input_name) 

# This function should convert any variant ID to a ProtVar ID.
# Protvar ID structure for genotypes is "chr pos ref alt", e.g. "1 10393 A T"
# AZ uses '-' while clinvar uses "_"
input_vars = vars[col_vars]\
    .drop_duplicates()\
    .str.split(separator,expand=True)
input_vars.columns = ["chrom","pos","ref","alt"]
partition_size = 1000
length_df = len(input_vars.index)

split_df = np.array_split(input_vars, ceil(length_df/partition_size))

# Write to csv
for counter, partition in enumerate(split_df):
    partition['id'] = "."
    partition = partition[["chrom", "pos", "id", "ref", "alt"]]
    partition.to_csv(
        f"./rare_vars-{counter}.csv", 
        sep = " ",
        header=False,
        index=False)
