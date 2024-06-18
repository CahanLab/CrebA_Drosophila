import numpy as np 
import pandas as pd 
import os

flybase_list = []
gene_list = []

with open('../input/reference_genome/dmel-all-r6.33.gtf', 'r') as file:
#with open(snakemake.input[0], 'r') as file:
    # Read and print each line
    for line in file:
        flybase_id = line.strip().split("\t")[8].split(";")[0]
        gene_id = line.strip().split("\t")[8].split(";")[1]

        gene_id = gene_id.removeprefix(" gene_symbol ").replace("\"", "")
        flybase_id = flybase_id.removeprefix("gene_id ").replace("\"", "")
        flybase_list.append(flybase_id)
        gene_list.append(gene_id)

convert_pd = pd.DataFrame()
convert_pd['gene_names'] = gene_list
convert_pd['flybase'] = flybase_list

convert_pd = convert_pd.drop_duplicates()

if os.path.isdir("../input/flybase_gene_conversion/") == False: 
    os.makedirs("../input/flybase_gene_conversion/")

convert_pd.to_csv("../input/flybase_gene_conversion/conversion_tab.csv", index=False)
#convert_pd.to_csv(snakemake.output[0], index=False)