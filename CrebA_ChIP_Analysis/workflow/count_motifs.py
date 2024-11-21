import numpy as np 
import pandas as pd
import argparse

# 
parser = argparse.ArgumentParser(description='check core motif location')
parser.add_argument('--motif', type=str, help = 'motif seq')
parser.add_argument('--fasta', type=str, help = 'the path for fasta file')
parser.add_argument('--out', type=str, help = 'the output path')
args = parser.parse_args()

motif = args.motif
fasta_file_path = args.fasta
output_path = args.out

##### make all the helper functions ######
comp_dict = dict()
comp_dict['A'] = 'T'
comp_dict['C'] = 'G'
comp_dict['G'] = 'C'
comp_dict['T'] = 'A'

def get_complement(motif):
    comp_motif = ''
    for nuc in motif: 
        comp_motif = comp_motif + comp_dict[nuc]
    return comp_motif

def get_motif_loc(motif, seq):
    locations = list()
    for i in range(len(seq)):
        if seq[i:i+len(motif)] == motif: 
            locations.append(str(i))
    return ",".join(locations)

def parse_fasta_to_dict(file_path):
    fasta_dict = {}
    with open(file_path, "r") as file:
        header = None
        sequence = None
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if header:  # Save the previous header and sequence
                    fasta_dict[header] = sequence
                header = line[1:]  # Remove the ">" character
                sequence = None  # Reset sequence collector
            else:  # Sequence line
                sequence = line
        if header:  # Save the last sequence
            fasta_dict[header] = sequence
    return fasta_dict

rev_motif = motif[::-1]
comp_motif = get_complement(motif)
comp_rev_motif = get_complement(rev_motif)

fasta_dict = parse_fasta_to_dict(fasta_file_path)
output_df = pd.DataFrame(columns = ['peaks', 'motif', 'rev_motif', 'comp_motif', 'comp_rev_motif'])

for tmp_peak in fasta_dict.keys():
    motif_loc = get_motif_loc(motif, fasta_dict[tmp_peak])
    rev_motif_loc = get_motif_loc(rev_motif, fasta_dict[tmp_peak])
    comp_motif_loc = get_motif_loc(comp_motif, fasta_dict[tmp_peak])
    comp_rev_motif_loc = get_motif_loc(comp_rev_motif, fasta_dict[tmp_peak])
    output_df.loc[len(output_df)] = [tmp_peak, motif_loc, rev_motif_loc, comp_motif_loc, comp_rev_motif_loc]

output_df.to_csv(output_path, index=False)