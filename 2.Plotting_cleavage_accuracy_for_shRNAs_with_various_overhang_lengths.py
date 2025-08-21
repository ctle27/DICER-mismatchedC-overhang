# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 14:38:29 2025

@author: congt
"""

path = 'D:/HKUST_Research/01.Research_projects/HT-DICER-processing-projects/DICER-lower stem project/ngs_analysis/'
import pandas as pd
processed_df_wt = pd.read_csv(path+'processed_df_wt_tlr1627_2reps.bed',sep='\t')
processed_df_wt['Symmetric_structure'] = processed_df_wt['concrete_struct'].apply(lambda strc: 'Yes' if 'A' not in strc and 'B' not in strc else 'No')
processed_df_wt.drop_duplicates(subset=['Variant','shRNA_sequence','5p-3p-alternative'], keep='first', inplace=True)

'''
Note:
It is important to add all cleavage sites for each shRNA sequence, even though the value is 0.
For example:
    shRNA-a has 3 clv sites: 21-21, 22-22, and other. Sum of accuracy is 1.
    Need to add other cleavage site: 19-19, 20-20, 23-23 with accuracy score = 0, so that the number of cleavage sites is similar for all variants.
    It affects the calculate of median value in later steps if the number of cleavage sites are not consistent among variants.
    Example: in DICER-WT there are clv site at 21-21 detected in 4000 variants, mean of 21-21 accuracy is 0.6, median is 0.5
    if DICERdeldsRBD there are 1000 variants with clv site at 21-21 detected in only 1000 variants, mean is 0.1 (expected),
    but median was 0.5 (due to too little variants with < 0.5 DC21 score)
    --> need to add 0 value for DC21 score in other 3000 variants as in DICER-WT --> median value now is 0.1 --> reflect correct pattern 
    
Do similarly for cleavage efficiency but add a pseudo efficiency of -50 (indicating very low efficiency)
'''
def adding_missing_cleavage_site(df_input, metrics, fillin_value):
    df = df_input.copy()
    df = df[['shRNA_sequence', 'concrete_struct', 'New_define_structure_2', 'Symmetric_structure', '5p-3p-alternative', metrics]]
    df.drop_duplicates(subset=['shRNA_sequence','5p-3p-alternative'], keep='first', inplace=True)
    df = df.pivot(index=['shRNA_sequence', 'concrete_struct', 'New_define_structure_2', 'Symmetric_structure'], columns='5p-3p-alternative', values=metrics)
    df.fillna(fillin_value, inplace=True)
    df.reset_index(inplace=True)
    df = pd.melt(df, id_vars=['shRNA_sequence', 'concrete_struct', 'New_define_structure_2', 'Symmetric_structure'], value_vars=['19-19', '20-20', '21-21', '22-22', '23-23', 'other'],
                        var_name='5p-3p-alternative', value_name=metrics)
    df.sort_values(['shRNA_sequence', '5p-3p-alternative'], ascending=[True, True], inplace=True)
    df.reset_index(inplace=True, drop=True)
    
    return (df)

df_wt = adding_missing_cleavage_site(processed_df_wt, 'Mean_Cleavage_accuracy_of_alternative_5p_3p', 0).merge(
    adding_missing_cleavage_site(processed_df_wt, 'Mean_Positional_efficiency_of_alternative_5p_3p', -50),
    on=['shRNA_sequence', 'concrete_struct', 'New_define_structure_2', 'Symmetric_structure', '5p-3p-alternative'],
    how='inner')

df_wt['5p_overhang'] = df_wt['New_define_structure_2'].str.count('F')

#%%checking cleavage accuracy of structures with different overhang length (box plot)
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
path = 'D:/HKUST_Research/01.Research_projects/HT-DICER-processing-projects/DICER-lower stem project/ngs_analysis/'
df = df_wt.copy()
sample = 'DICERWT'
clv_site = '21-21'
checking = 'Mean_Positional_efficiency_of_alternative_5p_3p' #Cleavage_accuracy or Positional_efficiency
if 'accuracy' in checking:
    title = 'accuracy'
if 'efficiency' in checking:
    title = 'efficiency'
    
df = df[(df['Symmetric_structure'] == 'Yes') & (df['5p-3p-alternative'] == clv_site)]
strc_list = ['0F; 23-SSSSSSS SSSSSSS-23; 36-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-36; 2T',
             '1F; 22-SSSSSSS SSSSSSS-22; 35-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-35; 3T',
             '2F; 21-SSSSSSS SSSSSSS-21; 34-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-34; 4T',
             '3F; 20-SSSSSSS SSSSSSS-20; 33-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-33; 5T',
             '4F; 19-SSSSSSS SSSSSSS-19; 32-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-32; 6T',
             '5F; 18-SSSSSSS SSSSSSS-18; 31-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-31; 7T',
             '6F; 17-SSSSSSS SSSSSSS-17; 30-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-30; 8T']

tri_nu_comb = []
for nu1 in ['A','T','G','C']:
     for nu2 in ['A','T','G','C']:
         for nu3 in ['A','T','G','C']:
             tri_nu_comb.append(nu1+nu2+nu3)
ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
seq = ''
seq_list = []
sub_list = []
tri5p = []
tri3p = []

for i in range(1,13):
    for comb1 in tri_nu_comb:
        for comb2 in tri_nu_comb:
            seq = ini_seq[:i] + comb1 + ini_seq[i+3:99-i] + comb2 + ini_seq[102-i:]
            seq_list.append(seq.replace('N'*32,''))
            sub_list.append(i)
            tri5p.append(comb1)
            tri3p.append(comb2)
df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
df_subgroup['Subgroup'] = sub_list
df_subgroup['shRNA_sequence'] = seq_list
df_subgroup['tri5p'] = tri5p
df_subgroup['tri3p'] = tri3p
df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner")
df_subgroup = df_subgroup[df_subgroup['Subgroup'].isin([1,2,3])]
df_subgroup = df_subgroup[df_subgroup['concrete_struct'].isin(strc_list)]
df_subgroup.reset_index(inplace = True, drop = True)


for sub in range(1,4):
    df_plot = df_subgroup[df_subgroup['Subgroup'] == sub]
    order = [item for item in strc_list if item in df_plot['concrete_struct'].unique()]
    if sub == 1:
        figsize = (5,3)
    if sub == 2:
        figsize = (4,3)
    if sub == 3:
        figsize = (3,3)  
    
    ax = plt.figure(figsize=figsize)
    mpl.rcParams['axes.linewidth'] = 1 #set the value globally
    mpl.rcParams['axes.spines.right'] = True
    mpl.rcParams['axes.spines.top'] = True
    ax = sns.boxplot(data= df_plot, x='concrete_struct', y=checking, zorder=10, order=order, color='royalblue', showfliers=False)
    ax = sns.stripplot(data=df_plot, x='concrete_struct', y=checking, color='royalblue', order=order,
                       size=2, dodge=True, edgecolor=".2", zorder=20)

    # plt.ylim(-0.1,1.1)
    # plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
    
    plt.ylim(-6.5,6.5)
    plt.yticks([-6,-3,0,3,6])
    
    
    ax.tick_params(axis='y', width = 1, length=8)
    ax.tick_params(axis='x', width = 1, length=8)
    plt.grid(axis = 'y', color = 'black', linestyle = '--', linewidth = 0.5, zorder=1)

    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(visible=False,size=8, rotation = 90)
    plt.yticks(visible=False)
    plt.savefig(path+f'plot/box_plot_{title}_subgroup{str(sub)}_{sample}_{clv_site}_different_overhang_length.png', dpi=150, bbox_inches='tight')
    plt.show()

#%%checking cleavage accuracy of each structures in each subgroup (heatmap)
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

path = 'D:/HKUST_Research/01.Research_projects/HT-DICER-processing-projects/DICER-lower stem project/ngs_analysis/'
df = df_wt.copy()

sample = 'DicerWT'
option = 'symmetric' #or asymmetric or symmetric

if option == 'symmetric':
    df = df[(df['Symmetric_structure'] == 'Yes')]
    folder = '23L-structures'

if option == 'asymmetric':
    df = df[~((df['Symmetric_structure'] == 'Yes'))]
    folder = 'other_structures'

checking = 'Mean_Cleavage_accuracy_of_alternative_5p_3p' #Cleavage_accuracy or Positional_efficiency
if 'accuracy' in checking:
    title = 'accuracy'
if 'efficiency' in checking:
    title = 'effiicency'

tri_nu_comb = []
for nu1 in ['A','T','G','C']:
     for nu2 in ['A','T','G','C']:
         for nu3 in ['A','T','G','C']:
             tri_nu_comb.append(nu1+nu2+nu3)
ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
seq = ''
seq_list = []
sub_list = []
tri5p = []
tri3p = []

for i in range(1,13):
    for comb1 in tri_nu_comb:
        for comb2 in tri_nu_comb:
            seq = ini_seq[:i] + comb1 + ini_seq[i+3:99-i] + comb2 + ini_seq[102-i:]
            seq_list.append(seq.replace('N'*32,''))
            sub_list.append(('TLR'+str(i+15)))
            tri5p.append(comb1)
            tri3p.append(comb2)
df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
df_subgroup['Subgroup'] = sub_list
df_subgroup['shRNA_sequence'] = seq_list
df_subgroup['tri5p'] = tri5p
df_subgroup['tri3p'] = tri3p
df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner")

strc_list = list(set(df_subgroup['concrete_struct'].tolist())) 

for tlr in ['TLR16','TLR17','TLR18','TLR19','TLR20','TLR21','TLR22','TLR23','TLR24','TLR25','TLR26','TLR27']:
    df_check = df_subgroup[df_subgroup['Subgroup'] == tlr]
    for strc in strc_list:
        df_plot = df_check[df_check['concrete_struct'] == strc]
        df_plot = df_plot.copy()
        df_plot.drop_duplicates(subset=['shRNA_sequence','5p-3p-alternative'],keep='first',inplace=True)
        l = len(set(df_plot['shRNA_sequence'].tolist()))
        if l > 29:
            df_plot = df_plot.pivot(index='shRNA_sequence',columns='5p-3p-alternative',values=checking)
            
            df_plot = df_plot.reindex(columns=['19-19','20-20','21-21','22-22','23-23','other'], fill_value=0)
            df_plot.fillna(0, inplace=True)
            df_plot.sort_values(['21-21'],ascending=True,inplace=True)
            
            ax = plt.figure(figsize=(2,1))
            mpl.rcParams['axes.linewidth'] = 0.25 #set the value globally
            mpl.rcParams['axes.spines.right'] = True
            mpl.rcParams['axes.spines.top'] = True
            my_color = sns.color_palette("Blues", as_cmap=True)
            # #for accuracy
            ax = sns.heatmap(data=df_plot,cmap=my_color,vmax=1,vmin=0, cbar=False)
            #for efficiency
            # ax = sns.heatmap(data=df_plot,cmap=my_color,vmax=5,vmin=-10, cbar=False)
            ax.tick_params(axis='y', width = 0, length=0)
            ax.tick_params(axis='x', width = 0, length=0)
            for _, spine in ax.spines.items():
                spine.set_visible(True)
            plt.xlabel('')
            plt.ylabel('')
            plt.yticks([])
            plt.xticks([])
            strc_short = strc.replace('L','')

            plt.savefig(path+f'heatmap/{sample}/{folder}/heatmap_{title}_{tlr}_{sample}_{strc_short}_no. of variant = {str(l)}.png', dpi=150, bbox_inches='tight')
            plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    