import pandas as pd
import numpy as np 
import os

###### to make the working directory ###### 
output_path = "../output/match_manual_automated_peaks"
if os.path.isdir(output_path) == False: 
    os.makedirs(output_path)    

##### this is the function to reproduce data ##### 
def find_bound_genes(peak_path, DE_path, closest_gene = False, dist_limit = 5000):
    narrow_peaks = pd.read_csv(peak_path, index_col=0)
    narrow_peaks = narrow_peaks.loc[narrow_peaks['nearest_gene'].isna() == False, :]

    if closest_gene == True: 
        bound_genes = np.unique(list(narrow_peaks['nearest_gene']))
    else: 
        bound_genes = list(narrow_peaks.loc[narrow_peaks['nearest_distance_1'] <= dist_limit, "nearest_gene_1"]) + list(narrow_peaks.loc[narrow_peaks['nearest_distance_2'] <= dist_limit, 'nearest_gene_2'])
        bound_genes = np.unique(bound_genes)
    
    output_df = pd.DataFrame(index = ['manual_bound_all', 'manual_unbound_all', 'manual_bound_SG', 'manual_unbound_SG'], columns = ['down_genes', 'up_genes', 'static_genes'])
    for i in range(0, 3): 
        if i == 0: 
            type_gene = 'down_genes'
        elif i == 1: 
            type_gene = 'up_genes'
        elif i == 2: 
            type_gene = 'static_genes'

        DE_df = pd.read_excel(DE_path, sheet_name = i)
        DE_df['bound_status'] = 0
        DE_df.loc[DE_df['target_genes'].isin(bound_genes), 'bound_status'] = 1

        sub_DE = DE_df.loc[[True if 'Unbound' in x else False for x in DE_df['Bound Status']], :]
        output_df.loc['manual_unbound_all', type_gene] = np.sum(sub_DE['bound_status']) / sub_DE.shape[0]

        sub_DE = sub_DE.loc[sub_DE['Salivary Gland'] != 'other', :]
        output_df.loc['manual_unbound_SG', type_gene] = np.sum(sub_DE['bound_status']) / sub_DE.shape[0]

        sub_DE = DE_df.loc[[True if 'Bound' in x else False for x in DE_df['Bound Status']], :]
        output_df.loc['manual_bound_all', type_gene] = np.sum(sub_DE['bound_status']) / sub_DE.shape[0]

        sub_DE = sub_DE.loc[sub_DE['Salivary Gland'] != 'other', :]
        output_df.loc['manual_bound_SG', type_gene] = np.sum(sub_DE['bound_status']) / sub_DE.shape[0]

    return output_df 

###### bidirection unique genes ######
peak_path = "../output/match_nearest_gene/oregon_fkh_sage_unique_nearest_genes.csv"
DE_path = "../input/manual_curated_DE_genes/manual_curated_DE.xlsx"
closest_gene = False
dist_limit = 5000
output_df = find_bound_genes(peak_path, DE_path, closest_gene, dist_limit)
print(output_df)

narrow_peaks = pd.read_csv(peak_path, index_col=0)
narrow_peaks = narrow_peaks.loc[narrow_peaks['nearest_gene'].isna() == False, :]
if closest_gene == True: 
    bound_genes = np.unique(list(narrow_peaks['nearest_fly_id']))
else: 
    bound_genes = list(narrow_peaks.loc[narrow_peaks['nearest_distance_1'] <= dist_limit, "nearest_fly_id_1"]) + list(narrow_peaks.loc[narrow_peaks['nearest_distance_2'] <= dist_limit, 'nearest_fly_id_2'])
    bound_genes = np.unique(bound_genes)
print(len(bound_genes))

spcg_df = pd.read_excel('../../SPCG_files/SPCG List.xlsx')
np.sum(spcg_df['Drosophila FBgn'].isin(bound_genes)) / spcg_df.shape[0]

###### bidirection intersecting genes ######
peak_path = "../output/match_nearest_gene/oregon_fkh_sage_intersect_nearest_genes.csv"
DE_path = "../input/manual_curated_DE_genes/manual_curated_DE.xlsx"
closest_gene = False
dist_limit = 5000
output_df = find_bound_genes(peak_path, DE_path, closest_gene, dist_limit)
print(output_df)

narrow_peaks = pd.read_csv(peak_path, index_col=0)
narrow_peaks = narrow_peaks.loc[narrow_peaks['nearest_gene'].isna() == False, :]
if closest_gene == True: 
    bound_genes = np.unique(list(narrow_peaks['nearest_fly_id']))
else: 
    bound_genes = list(narrow_peaks.loc[narrow_peaks['nearest_distance_1'] <= dist_limit, "nearest_fly_id_1"]) + list(narrow_peaks.loc[narrow_peaks['nearest_distance_2'] <= dist_limit, 'nearest_fly_id_2'])
    bound_genes = np.unique(bound_genes)
print(len(bound_genes))

spcg_df = pd.read_excel('../../SPCG_files/SPCG List.xlsx')
np.sum(spcg_df['Drosophila FBgn'].isin(bound_genes)) / spcg_df.shape[0]

###### this is just to explore data for unique peaks ######
narrow_peaks = pd.read_csv("../output/match_nearest_gene/oregon_fkh_sage_unique_nearest_genes.csv", index_col=0)
narrow_peaks = narrow_peaks.loc[narrow_peaks['nearest_gene'].isna() == False, :]

down_df = pd.read_excel("../input/manual_curated_DE_genes/manual_curated_DE.xlsx", sheet_name = 0)
down_df['bound_status'] = 0
down_df.loc[down_df['target_genes'].isin(narrow_peaks['nearest_gene']), 'bound_status'] = 1
down_df.to_csv(os.path.join(output_path, 'oregon_fkh_sage_unique_bound_down.csv'))
sub_down = down_df.loc[[True if 'Unbound' in x else False for x in down_df['Bound Status']], :]
np.sum(sub_down['bound_status']) / sub_down.shape[0]

sub_down = sub_down.loc[sub_down['Salivary Gland'] != 'other', :]
np.sum(sub_down['bound_status']) / sub_down.shape[0]

sub_down = down_df.loc[[True if 'Bound' in x else False for x in down_df['Bound Status']], :]
np.sum(sub_down['bound_status']) / sub_down.shape[0]

sub_down = sub_down.loc[sub_down['Salivary Gland'] != 'other', :]
np.sum(sub_down['bound_status']) / sub_down.shape[0]

up_df = pd.read_excel("../input/manual_curated_DE_genes/manual_curated_DE.xlsx", sheet_name = 1)
up_df['bound_status'] = 0
up_df.loc[up_df['target_genes'].isin(narrow_peaks['nearest_gene']), 'bound_status'] = 1
up_df.to_csv(os.path.join(output_path, 'oregon_fkh_sage_unique_bound_up.csv'))
sub_up = up_df.loc[[True if 'Unbound' in x else False for x in up_df['Bound Status']], :]
np.sum(sub_up['bound_status']) / sub_up.shape[0]

sub_up = sub_up.loc[sub_up['Salivary Gland'] != 'other', :]
np.sum(sub_up['bound_status']) / sub_up.shape[0]

sub_up = up_df.loc[[True if 'Bound' in x else False for x in up_df['Bound Status']], :]
np.sum(sub_up['bound_status']) / sub_up.shape[0]

sub_up = sub_up.loc[sub_up['Salivary Gland'] != 'other', :]
np.sum(sub_up['bound_status']) / sub_up.shape[0]

static_df = pd.read_excel("../input/manual_curated_DE_genes/manual_curated_DE.xlsx", sheet_name = 2)
static_df['bound_status'] = 0
static_df.loc[static_df['target_genes'].isin(narrow_peaks['nearest_gene']), 'bound_status'] = 1
static_df.to_csv(os.path.join(output_path, 'oregon_fkh_sage_unique_bound_static.csv'))

sub_static = static_df.loc[[True if 'Unbound' in x else False for x in static_df['Bound Status']], :]
np.sum(sub_static['bound_status']) / sub_static.shape[0]

sub_static = sub_static.loc[sub_static['Salivary Gland'] != 'other', :]
np.sum(sub_static['bound_status']) / sub_static.shape[0]

sub_static = static_df.loc[[True if 'Bound' in x else False for x in static_df['Bound Status']], :]
np.sum(sub_static['bound_status']) / sub_static.shape[0]

sub_static = sub_static.loc[sub_static['Salivary Gland'] != 'other', :]
np.sum(sub_static['bound_status']) / sub_static.shape[0]

##### SPCG ###### 
spcg_df = pd.read_excel('../../SPCG_files/SPCG List.xlsx')
np.sum(spcg_df['Drosophila FBgn'].isin(narrow_peaks['nearest_fly_id']))

###### this is just to explore data for intersect peaks ######
narrow_peaks = pd.read_csv("../output/match_nearest_gene/oregon_fkh_sage_intersect_nearest_genes.csv", index_col=0)
narrow_peaks = narrow_peaks.loc[narrow_peaks['nearest_gene'].isna() == False, :]

down_df = pd.read_excel("../input/manual_curated_DE_genes/manual_curated_DE.xlsx", sheet_name = 0)
down_df['bound_status'] = 0
down_df.loc[down_df['target_genes'].isin(narrow_peaks['nearest_gene']), 'bound_status'] = 1
down_df.to_csv(os.path.join(output_path, 'oregon_fkh_sage_intersect_bound_down.csv'))
sub_down = down_df.loc[[True if 'Unbound' in x else False for x in down_df['Bound Status']], :]
np.sum(sub_down['bound_status']) / sub_down.shape[0]

sub_down = sub_down.loc[sub_down['Salivary Gland'] != 'other', :]
np.sum(sub_down['bound_status']) / sub_down.shape[0]

sub_down = down_df.loc[[True if 'Bound' in x else False for x in down_df['Bound Status']], :]
np.sum(sub_down['bound_status']) / sub_down.shape[0]

sub_down = sub_down.loc[sub_down['Salivary Gland'] != 'other', :]
np.sum(sub_down['bound_status']) / sub_down.shape[0]

up_df = pd.read_excel("../input/manual_curated_DE_genes/manual_curated_DE.xlsx", sheet_name = 1)
up_df['bound_status'] = 0
up_df.loc[up_df['target_genes'].isin(narrow_peaks['nearest_gene']), 'bound_status'] = 1
up_df.to_csv(os.path.join(output_path, 'oregon_fkh_sage_intersect_bound_up.csv'))
sub_up = up_df.loc[[True if 'Unbound' in x else False for x in up_df['Bound Status']], :]
np.sum(sub_up['bound_status']) / sub_up.shape[0]

sub_up = sub_up.loc[sub_up['Salivary Gland'] != 'other', :]
np.sum(sub_up['bound_status']) / sub_up.shape[0]

sub_up = up_df.loc[[True if 'Bound' in x else False for x in up_df['Bound Status']], :]
np.sum(sub_up['bound_status']) / sub_up.shape[0]

sub_up = sub_up.loc[sub_up['Salivary Gland'] != 'other', :]
np.sum(sub_up['bound_status']) / sub_up.shape[0]

static_df = pd.read_excel("../input/manual_curated_DE_genes/manual_curated_DE.xlsx", sheet_name = 2)
static_df['bound_status'] = 0
static_df.loc[static_df['target_genes'].isin(narrow_peaks['nearest_gene']), 'bound_status'] = 1
static_df.to_csv(os.path.join(output_path, 'oregon_fkh_sage_intersect_bound_static.csv'))

sub_static = static_df.loc[[True if 'Unbound' in x else False for x in static_df['Bound Status']], :]
np.sum(sub_static['bound_status']) / sub_static.shape[0]

sub_static = sub_static.loc[sub_static['Salivary Gland'] != 'other', :]
np.sum(sub_static['bound_status']) / sub_static.shape[0]

sub_static = static_df.loc[[True if 'Bound' in x else False for x in static_df['Bound Status']], :]
np.sum(sub_static['bound_status']) / sub_static.shape[0]

sub_static = sub_static.loc[sub_static['Salivary Gland'] != 'other', :]
np.sum(sub_static['bound_status']) / sub_static.shape[0]

##### SPCG ###### 
spcg_df = pd.read_excel('../../SPCG_files/SPCG List.xlsx')
np.sum(spcg_df['Drosophila FBgn'].isin(narrow_peaks['nearest_fly_id']))