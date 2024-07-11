import pandas as pd
import numpy as np 
import os
# the identification of tss regions is from this paper 
# https://academic.oup.com/bioinformatics/article/38/20/4806/6674500

def get_max_length(file_path, chr_prefix = False):
    output_dict = dict()
    with open(file_path, 'r') as file: 
        for line in file: 
            if line.startswith(">"):
                if chr_prefix == False:
                    output_dict[line.split(" ")[0][1:]] = line.split(" ")[6].removeprefix('length=').removesuffix(";")
                else: 
                    output_dict['chr' + line.split(" ")[0][1:]] = line.split(" ")[6].removeprefix('length=').removesuffix(";")
    return output_dict

gtf_table = pd.read_csv("../input/reference_genome/dmel-all-r6.33.gtf", sep = '\t', header = None)
gtf_table = gtf_table.loc[gtf_table[2] == 'gene', :]

# read in the fasta file to get the length 

max_dict = get_max_length('../input/reference_genome/dmel-all-chromosome-r6.33.fasta')

# the 1000 upstream and 300 downstream for the promoter region was found from 
# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2005-6-8-r72
# I think the above paper was using 700 - 300 

# I think traditionally people use 1kb upstream of tss 
# https://www.genome.gov/Multimedia/Slides/ENCODE2015-ResearchAppsUsers/23_ENCODE-Factorbook_Purcaro.pdf
# https://academic.oup.com/nar/article/35/2/406/2400679?login=false
w = 1000
w_end = 300 # this is just arbritary 
tss_table = pd.DataFrame(data = None, index = gtf_table.index, columns = ['chrom', 'strand', 'promoter_start', 'promoter_end', 'tss', 'gene_id', 'gene_name', 'gene_length'])
for tmp_index in gtf_table.index: 
    chrom = gtf_table.loc[tmp_index, :][0]
    strand = gtf_table.loc[tmp_index, :][6]
    if strand == "+":
        tss = gtf_table.loc[tmp_index, :][3] - 1
        promoter_start = np.max([tss - w, 0])
        promoter_end = np.min([tss + w_end, int(max_dict[chrom])])
    else:
        tss = gtf_table.loc[tmp_index, :][4] - 1
        promoter_start = np.min([tss + w, int(max_dict[chrom])])
        promoter_end = np.max([tss - w_end, 0])
        
    gene_strings = gtf_table.loc[tmp_index, :][8]
    gene_strings = gene_strings.split('"')
    gene_id = gene_strings[1]
    gene_name = gene_strings[3]
    gene_length = abs(gtf_table.loc[tmp_index, :][4] - gtf_table.loc[tmp_index, :][3])
    row_insert = [chrom, strand, promoter_start, promoter_end, tss, gene_id, gene_name, gene_length]
    tss_table.loc[tmp_index, :] = row_insert

if os.path.isdir("../output/tss_table") == False: 
    os.makedirs("../output/tss_table")
tss_table.to_csv("../output/tss_table/tss_table.txt", index = False)

for i in tss_table.index: 
    print(i)
    if tss_table.loc[i, 'promoter_end'] < tss_table.loc[i, 'promoter_start']:
        temp_end = tss_table.loc[i, 'promoter_start']
        tss_table.loc[i, 'promoter_start'] = tss_table.loc[i, 'promoter_end']
        tss_table.loc[i, 'promoter_end'] = temp_end
tss_table.to_csv("../output/tss_table/promoter_regions.bed", index = False, header = False, sep = '\t')