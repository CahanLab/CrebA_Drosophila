import pandas as pd 
import numpy as np 

creba_motifs = open("../input/previous_CrebA_motifs/AllCrebAMotifs_probabilitymatricies.txt", 'r')
Lines = creba_motifs.readlines()

# store all the probability matrix in a dictionary 
motif_dict = dict()
for line in Lines:
    if line[0] == ">":
        cur_motif = line.strip()
        motif_dict[cur_motif] = []
    else: 
        if len(line.strip()) == 0:
            continue
        motif_dict[cur_motif].append(line.strip().split("\t"))

with open('../input/previous_CrebA_motifs/CrebAMotifs_meme.meme', 'w') as file:
    file.write('MEME version 5.4.1 \n') 
    file.write("\n")
    file.write("ALPHABET= ACGT \n")
    file.write("\n")
    file.write("strands: + - \n")
    file.write("\n")
    file.write("Background letter frequencies (from uniform background): \n")
    file.write("A 0.25000 C 0.25000 G 0.25000 T 0.25000 \n")
    file.write("\n")
    for temp_motif in motif_dict.keys():
        print(temp_motif)    
        temp_prob = motif_dict[temp_motif]
        prob_df = pd.DataFrame(temp_prob).T
        prob_df.columns = prob_df.iloc[0, :]
        prob_df = prob_df.iloc[1:, :].copy()
        prob_df = prob_df.astype(float)
        prob_df = prob_df.loc[:, ['A', 'C', 'G', 'T']]
        row_sums = prob_df.sum(axis=1)

        # Calculate percentage
        prob_df = prob_df.div(row_sums, axis=0)
        prob_df = prob_df.round(6)
        prob_df = prob_df.astype(str)
        temp_motif = temp_motif.strip(">")
        file.write("MOTIF " + temp_motif + " " + temp_motif + "\n")
        file.write("\n")
        file.write("letter-probability matrix: alength= 4 w= " + str(prob_df.shape[0]) +"\n")
        for temp_row in prob_df.to_numpy():
            file.write("\t".join(temp_row) + " \n")
        file.write("\n")
        

