# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 18:14:01 2024

@author: congt
"""
'''
Obtain data of SNP frequency, its effect on changing DICER cutting site from running programs in hpc3 server
'''

import pandas as pd
path = 'D:/HKUST_Research/isomirdb/ngs_analysis/'

df1 = pd.read_csv(path + 'Cleavage_altering_score_of_mutated_premiRNAs.bed', sep='\t')
df1 = df1[df1['Mutated_position'] != 'Multiple_location']

df1['Mutated_position'] = pd.to_numeric(df1['Mutated_position'], errors='coerce')
df1['Position_from_5p_end_of_3pmiRNA'] = df1['Mutated_position'] - df1['3pmiRNA_starting'] + 1

'''
The column 'Distance from 5p end' displays the distance (in nt) from the 1st nt of 3p miRNA (miRBase annotation)
Example: nt at position 2 on the 3p miRNA having the distance of 1 nt.

Note: only keep SNPs/mutations that are not mapped to other WT homolog pre-miRNAs
    SOME reads of hsa-mir-92b_rs1350664926_T 4p miRNA were mapped to both miR-92a-3p (from hsa-mir-92a-1, hsa-mir-92a-2) and mir-92b_rs1350664926_T-3p
    --> they were not assigned to hsa-mir-92b_rs1350664926_T
    --> Remove these two SNPs in calculation
    
    Reads mapped to hsa-let-7a-1/7a-3-3p were also mapped to 7f-2-88G>A mutation --> remove this mutation
    Reads mapped to hsa-mir-92b-81C>T also mapped to 92a-1, 92a-2 --> remove this mutation
    Reads mapped to hsa-mir-1269a 85G>A also mapped to -1269b --> remove this mutation
    
'''

'''
Obtain new define structures of WT and mutant pre-miRNAs from hpc3 programs
Only keep mutation that retain the stem length of the pre-miRNA for further analsysis
'''

df2 = pd.read_csv(path + '3.WT_mutated_pre_miRNA_new_define_structure.bed', sep='\t')
'''
replace 'S, T, F' by 'M' and compare new_define_structure1. expect the structures to be the same for ref and alt pre-mirna
'''
for i,struct_ref in enumerate(df2['New_define_structure1_Ref']):
    struct_alt = df2['New_define_structure1_Alt'][i]
    
    struct_ref = struct_ref.replace('S','M')
    struct_ref = struct_ref.replace('F','M')
    struct_ref = struct_ref.replace('T','M')
    
    struct_alt = struct_alt.replace('S','M')
    struct_alt = struct_alt.replace('F','M')
    struct_alt = struct_alt.replace('T','M')
    
    if struct_ref == struct_alt:
        df2.loc[i,'Structure_change'] = 'No'
    if struct_ref != struct_alt:
        df2.loc[i,'Structure_change'] = 'Yes'
    
    
df2 = df2[df2['Structure_change'] == 'No']
df2 = df2[['MiRBase_ID', 'Mutation_profile', 'New_define_structure2_Ref', 'New_define_structure2_Alt', 'New_define_structure_wobble_Ref', 'New_define_structure_wobble_Alt']]

df3 = df1.merge(df2, on=['MiRBase_ID', 'Mutation_profile'], how='inner')

'''
Remove some mutations as described above
'''
df3['Mutation_ID'] = df3['MiRBase_ID'] + '-' + df3['Mutation_profile']
df3 = df3[~(df3['Mutation_ID'].isin(['hsa-let-7f-2-88G>A', 'hsa-mir-92b-81C>T', 'hsa-mir-1269a-85G>A']))]
df3 = df3[df3['Position_from_5p_end_of_3pmiRNA'] >= 4]

df3.reset_index(inplace=True, drop=True)
#%%
df_check = df3[df3['Cleavage_altering_score'] >= 25]
#%%plot the position of SNPs


import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

ax = plt.figure(figsize=(8,4))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

ax = sns.stripplot(data=df3, x='Position_from_5p_end_of_3pmiRNA', y='Cleavage_altering_score', size=8,  edgecolor="black", zorder=1, legend=False,
                   color='darksalmon', linewidth=0.5)

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)

position_list = list(set(df3['Position_from_5p_end_of_3pmiRNA'].tolist()))
plt.xticks(range(len(position_list)), [int(pos) for pos in position_list])
plt.tight_layout()

# plt.title(primiRNA)
plt.ylim(-5,130)
plt.yticks([0,25,50,75,100,125])

# plt.ylabel('Cleavage altering score (min 0, max 200)')
# plt.xlabel('Position on 3p miRNA (miRBase annotation)')
plt.ylabel('')
plt.xlabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)

path_save_fig = 'D:/HKUST_Research/isomirdb/ngs_analysis/plot/'
plt.savefig(path_save_fig+'Cleavage_altering_score_of_mutant_premiRNAs.png', dpi=150, bbox_inches='tight')
plt.show()

#%%plot for the type of structure change by the mutations
df_plot = df3.copy()
for i,struct_ref in enumerate(df_plot['New_define_structure_wobble_Ref']):
    struct_alt = df_plot['New_define_structure_wobble_Alt'][i]
    
    mutated_position = df_plot['Mutated_position'][i] - 30
    ref_nt = struct_ref[mutated_position - 1].upper()
    alt_nt = struct_alt[mutated_position - 1].upper()
    
    if ref_nt == alt_nt:
        df_plot.loc[i, 'Structure_change'] = 'No-change'
    elif ref_nt + '>' + alt_nt in ['M>T','W>T','S>T']:
        df_plot.loc[i, 'Structure_change'] = 'Overhang_lengthening'
    elif ref_nt + '>' + alt_nt in ['T>M','T>W','T>S']:
        df_plot.loc[i, 'Structure_change'] = 'Overhang_shortening'
    else:
        df_plot.loc[i, 'Structure_change'] = ref_nt + '>' + alt_nt

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

ax = plt.figure(figsize=(7,4))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True


order = ['M>S', 'M>W', 'S>M', 'S>W', 'W>M', 'W>S', 'Overhang_lengthening', 'Overhang_shortening', 'No-change']
ax = sns.stripplot(data=df_plot, x='Cleavage_altering_score', y='Structure_change', size=8,  edgecolor="black", zorder=1, legend=False,
                   color='darksalmon', linewidth=0.5, order=order)

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)


# plt.title(primiRNA)
plt.xlim(-5,105)
plt.xticks([0,25,50,75,100])

plt.ylabel('')
plt.xlabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)

path_save_fig = 'D:/HKUST_Research/isomirdb/ngs_analysis/plot/'
plt.savefig(path_save_fig+'Cleavage_altering_score_of_mutant_premiRNAs_structure_change.png', dpi=150, bbox_inches='tight')
plt.show()

#report frequency of structure change
frequency = df_plot['Structure_change'].value_counts()
print (frequency)

#%%plot the count of SNP events at each position
#set minimum of Cleavage altering score to filter SNPs
'''
	Position	Count_SNPs
7	3.0	1
10	4.0	4

'''
threshold = 10

def count_SNPs(df_input, threshold):
    df1 = df_input.copy()
    df2 = df1[df1['Cleavage_altering_score'] >= threshold]
    df2 = df2[['Position', 'SNP_ID']]
    df2['Count_SNPs'] = df2.groupby(['Position'])['SNP_ID'].transform('count')
    del df2['SNP_ID']
    df2.drop_duplicates(subset=['Position'], keep='first', inplace=True)
    pos_list = df2['Position'].tolist()
    add_pos_list = [pos for pos in range(2,22) if pos not in pos_list]
    
    df3 = pd.DataFrame({'Position': add_pos_list, 'Count_SNPs': 0})
    df3 = pd.concat([df2,df3], ignore_index=True)
    df3 = df3.sort_values('Position').reset_index(drop=True)
    return (df3)
df2 = count_SNPs(df1, 10)
df3 = count_SNPs(df1, 20)

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

ax = plt.figure()
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

ax = sns.barplot(data=df2, x='Position', y='Count_SNPs', color='royalblue')
ax = sns.barplot(data=df3, x='Position', y='Count_SNPs', color='darksalmon')

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)

# plt.title(primiRNA)
plt.ylim(0,6)
plt.yticks([0,1,2,3,4,5,6])

position_list = list(set(df2['Position'].tolist()))
plt.xticks(range(len(position_list)), [int(pos) for pos in position_list])

plt.ylabel('Number of SNPs')
# plt.xlabel('Position on 3p miRNA (miRBase annotation)')


# path_save_fig = 'D:/HKUST_Research/isomirdb/isomiR_percentage/WT_vs_SNP/'
# plt.savefig(path_save_fig+f'/{primiRNA}_WT_vs_SNP_isomiR_frequency.png', dpi=150, bbox_inches='tight')
plt.show()



#%%check metadata of SNPs
import pandas as pd
path = 'D:/HKUST_Research/01.Research_projects/HT-DICER-processing-projects/DICER-lower stem project/isomirdb/ngs_analysis/'

df1 = pd.read_csv(path + '8.Metadata_of_mutations.bed', sep='\t')
df2 = pd.read_csv(path + '7.pattern_and_abundance_of_selected_isomiR.bed', sep='\t', usecols=['MiRBase_ID', 'Mutation_profile', 'Sample', 'Percentage_of_miRNA_type_in_the_sample'])

miRNA = 'hsa-mir-760'
mutation = '85G>A'

df1 = df1[(df1['MiRBase_ID'] == miRNA) & (df1['Mutation_profile'] == mutation)]
df2 = df2[(df2['MiRBase_ID'] == miRNA) & (df2['Mutation_profile'] == mutation)]
df2.drop_duplicates(subset=['MiRBase_ID', 'Mutation_profile', 'Sample'], keep='first', inplace=True)
#%%plot for percentage of mutant miRNA reads
df3 = df1.merge(df2, on=['MiRBase_ID','Mutation_profile', 'Sample'],how='inner')

median_values = df3.groupby('isomiRdb_tissue')['Percentage_of_miRNA_type_in_the_sample'].median().sort_values(ascending=False).index

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

ax = plt.figure(figsize=(6,3))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

ax = sns.boxplot(data=df3, x='isomiRdb_tissue', y='Percentage_of_miRNA_type_in_the_sample', color='royalblue', order=median_values, showfliers=False)
ax = sns.stripplot(data=df3, x='isomiRdb_tissue', y='Percentage_of_miRNA_type_in_the_sample', color='royalblue', order=median_values, linewidth=2)

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'y', color = 'black', linestyle = '--', linewidth = 0.2, zorder=0)

# plt.title(primiRNA)
plt.ylim(0,60)
plt.yticks([0,20,40,60])

plt.xticks(visible=True,size=8, rotation = 90)
plt.yticks(visible=False)

plt.ylabel('')
# plt.xlabel('')

mutation_name = mutation.replace('>','to')
path_save_fig = 'D:/HKUST_Research/Dicer-TLR16-27/isomirdb/ngs_analysis/plot/'
# plt.savefig(path_save_fig+f'{miRNA}_{mutation_name}_percentage_of_mutation_reads.png', dpi=150, bbox_inches='tight')
plt.show()

#%%plot for number of samples for each tissue

sample_count = df3['isomiRdb_tissue'].value_counts().to_frame()
sample_count.reset_index(inplace=True)

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

ax = plt.figure(figsize=(6,3))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

ax = sns.barplot(data=sample_count, x='isomiRdb_tissue', y='count', color='royalblue', zorder=10)

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'y', color = 'black', linestyle = '--', linewidth = 0.2, zorder=0)

# plt.title(primiRNA)
plt.ylim(0,20)
plt.yticks([0,5,10,15,20])

plt.xticks(visible=True,size=8,rotation=90)
plt.yticks(visible=False)

plt.ylabel('')
plt.xlabel('')

mutation_name = mutation.replace('>','to')
path_save_fig = 'D:/HKUST_Research/Dicer-TLR16-27/isomirdb/ngs_analysis/plot/'
# plt.savefig(path_save_fig+f'{miRNA}_{mutation_name}_sample_count_each_tissue.png', dpi=150, bbox_inches='tight')
plt.show()

import pandas as pd
df1 = pd.read_csv('7.pattern_and_abundance_of_selected_isomiR.bed', sep='\t')
'''
df1
['MiRBase_ID', '3pmiRNA_starting', 'Mapping_position', 'Mutation_profile', 'Sample', 'Sum_count_isomiR', 'Sum_count_miRNA', 'Percentage_of_cutting_site', 'Percentage_of_miRNA_type_in_the_sample', 'Count_sample_of_miRNA_type', 'Percentage_of_miRNA_type_in_wholepool']
'''

df2 = df1.copy()
primiRNA_list = list(set(df2['MiRBase_ID'].tolist()))

'''
create a dataframe containing information of the number of samples containing wildtype and each mutant pri-miRNAs.
'''
df_sample = df2[['MiRBase_ID', 'Mutation_profile', 'Count_sample_of_miRNA_type']]
df_sample = df_sample.copy()
df_sample.drop_duplicates(subset=['MiRBase_ID', 'Mutation_profile'], keep='first', inplace=True)

df_wt = df_sample[df_sample['Mutation_profile'] == 'Wildtype']
df_wt = df_wt.copy()
df_wt.rename(columns={'Count_sample_of_miRNA_type': 'Number_of_samples_for_Wildtype_primiRNA'}, inplace=True)
del df_wt['Mutation_profile']

df_mut = df_sample[df_sample['Mutation_profile'] != 'Wildtype']
df_mut = df_mut.copy()
df_mut.rename(columns={'Count_sample_of_miRNA_type': 'Number_of_samples_for_mutant_primiRNA'}, inplace=True)

df_sample = df_mut.merge(df_wt, on=['MiRBase_ID'], how='inner')

#%%calculate Cleavage_altering_score
'''
For each mutation:
1) Determine DICER cutting site matrix for the mutated pre-miRNA: Position - Percentage
2) Determine DICER cutting site matrix for the corresponding wildtype pre-miRNA: Position - Percentage
3) Calculate Cleavage_altering_score:
    sum(ABS(%CLi,SNP - %CLi,WT)), %CLi is the percentage of cutting site at position i
'''

df_combine = pd.DataFrame()
for primiRNA in primiRNA_list:
    df_temp2 = df2[df2['MiRBase_ID'] == primiRNA]
    cutting_site = df_temp2['Mapping_position'].unique() #all DICER cutting site on WT and SNPed pre-miRNAs
    
    '''
    add missing cutting sites
    For examples: WT pri-miRNA was cleaved at 2 sites, while SNP-ed pri-miRNA was cleaved at 4 sites.
    Need to add 2 missing cutting site in pri-miRNA with the Percentage of 0 %.
    '''
    
    df_temp2 = df_temp2[['MiRBase_ID', '3pmiRNA_starting', 'Mutation_profile', 'Sample', 'Mapping_position', 'Percentage_of_cutting_site']]
    df_temp2 = df_temp2.pivot(index=['MiRBase_ID', '3pmiRNA_starting', 'Mutation_profile', 'Sample'], columns='Mapping_position', values='Percentage_of_cutting_site')
    df_temp2.fillna(0, inplace=True)
    
    df_temp2.reset_index(inplace=True)
    df_temp2 = pd.melt(df_temp2, id_vars=['MiRBase_ID', '3pmiRNA_starting', 'Mutation_profile', 'Sample'], value_vars=cutting_site,
                        var_name='Mapping_position', value_name='Percentage_of_cutting_site')
    
    df_temp2['Mean_percentage'] = df_temp2.groupby(['MiRBase_ID', 'Mutation_profile', 'Mapping_position'])['Percentage_of_cutting_site'].transform('mean')
    df_temp2 = df_temp2.drop(['Sample', 'Percentage_of_cutting_site'], axis=1).drop_duplicates(subset=['MiRBase_ID', 'Mutation_profile', 'Mapping_position'])
    
    df_mut = df_temp2[df_temp2['Mutation_profile'] != 'Wildtype']
    df_wt = df_temp2[df_temp2['Mutation_profile'] == 'Wildtype']
    
    df_mut = df_mut.copy()
    df_wt = df_wt.copy()
    
    df_mut.rename(columns={'Mean_percentage':'Mean_percentage_mutant'}, inplace=True)
    df_wt = df_wt.rename(columns={'Mean_percentage': 'Mean_percentage_WT'}).drop('Mutation_profile', axis=1)
    
    if not df_mut.empty and not df_wt.empty:

        df_temp2 = df_mut.merge(df_wt, on=['MiRBase_ID', '3pmiRNA_starting', 'Mapping_position'], how='inner')
        df_temp2['Delta_percentage'] = abs(df_temp2['Mean_percentage_WT'] - df_temp2['Mean_percentage_mutant'])
        df_temp2['Cleavage_altering_score'] = df_temp2.groupby(['Mutation_profile'])['Delta_percentage'].transform('sum')
        
        df_temp2 = df_temp2.drop(['Mapping_position', 'Mean_percentage_WT', 'Mean_percentage_mutant', 'Delta_percentage'],
                                 axis=1).drop_duplicates(subset=['MiRBase_ID', '3pmiRNA_starting', 'Mutation_profile'])
        
        df_combine = pd.concat([df_combine, df_temp2], ignore_index=True)
'''
merge with df_sample to obtain information of number of samples containing wildtype and mutant pri-miRNAs
'''
df_combine = df_combine.merge(df_sample, on=['MiRBase_ID', 'Mutation_profile'], how='inner')

#%%
'''
Calculate the distance from mutation sites to the 5p end of the 3p miRNA
The distance measurement was only based on nt counting regardless of the pre-miRNA structure
Distance 1 corresponds to nt at position 2 from the reported miRBase DICER cutting site

If more than one mutations was found in the reads --> mark 'Multiple_location'
'''
df_combine2 = df_combine.copy()
df_combine2['Mutated_position'] = df_combine2['Mutation_profile'].apply(lambda x: int(''.join(filter(str.isdigit, x))) if x.count('>') == 1 else 'Multiple_location').astype('O')

df_combine2.to_csv('Cleavage_altering_score_of_mutated_premiRNAs.bed', sep='\t', index=False)

#%%
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

for primiRNA in primiRNA_list:
    df_plot1 = df2[df2['MiRBase_ID'] == primiRNA]
    
    df_plot1 = df_plot1[['MiRBase_ID', 'Mutation_profile', 'Sample', 'Mapping_position', 'Percentage_of_cutting_site', 'Percentage_of_miRNA_type_in_wholepool']]
    cutting_site = list(set(df_plot1['Mapping_position'].tolist()))
    cutting_site.sort()
    
    '''
    add missing cutting sites
    For examples: WT pri-miRNA was cleaved at 2 sites, while SNP-ed pri-miRNA was cleaved at 4 sites.
    Need to add 2 missing cutting site in pri-miRNA with the Percentage of 0 %.
    '''
    
    df_plot1 = df_plot1.pivot(index=['MiRBase_ID', 'Mutation_profile', 'Sample', 'Percentage_of_miRNA_type_in_wholepool'],
                              columns='Mapping_position', values='Percentage_of_cutting_site')
    df_plot1.fillna(0, inplace=True)
    
    df_plot1.reset_index(inplace=True)
    df_plot1 = pd.melt(df_plot1, id_vars=['MiRBase_ID', 'Mutation_profile', 'Sample', 'Percentage_of_miRNA_type_in_wholepool'], value_vars=cutting_site,
                        var_name='Mapping_position', value_name='Percentage_of_cutting_site')
    mutation_list = df_plot1['Mutation_profile'].unique()
    
    
    for mutation in mutation_list:
        if mutation != 'Wildtype':
            df_temp = df_combine[(df_combine['Mutation_profile'] == mutation) & (df_combine['MiRBase_ID'] == primiRNA)]
            df_temp.reset_index(inplace=True, drop=True)
            if not df_temp.empty and df_temp['Cleavage_altering_score'].iloc[0] >= 20: #only draw mutations if its cleavage_altering_score >= 20 units
                df_plot2 = df_plot1[df_plot1['Mutation_profile'].isin([mutation, 'Wildtype'])]
    
                frequency_mutation = df_plot2[df_plot2['Mutation_profile'] == mutation]['Percentage_of_miRNA_type_in_wholepool'].tolist()[0]
                frequency_mutation_rounded = round(frequency_mutation, 1)
                
                if 'Wildtype' in df_plot2['Mutation_profile'].tolist(): #make sure there are both SNP and WT pri-miRNA for comparison
                    '''
                    only keep samples that found both WT and mutant miRNA reads for comparison
                    '''
                    selected_samples = list(set(df_plot2[df_plot2['Mutation_profile'] == mutation]['Sample'].tolist()))
                    df_plot2 = df_plot2.copy()
                    df_plot2 = df_plot2[df_plot2['Sample'].isin(selected_samples)]
    
                    ax = plt.figure()
                    plt.rcParams['figure.max_open_warning'] = 20
                    mpl.rcParams['axes.linewidth'] = 1 #set the value globally
                    mpl.rcParams['axes.spines.right'] = True
                    mpl.rcParams['axes.spines.top'] = True
                    
                    hue_order = ['Wildtype', mutation]
                    
                    ax = sns.boxplot(data=df_plot2, x='Mapping_position', y='Percentage_of_cutting_site', hue='Mutation_profile', hue_order = hue_order,
                                     palette = ['lightseagreen', 'darksalmon'], zorder=10, showfliers=False)
                    ax = sns.stripplot(data=df_plot2, x='Mapping_position', y='Percentage_of_cutting_site',  hue='Mutation_profile', hue_order = hue_order, size=2,
                                       palette = ['lightseagreen', 'darksalmon'], edgecolor=".2", zorder=1, dodge=True, legend=False)
                    
                    ax.tick_params(axis='y', width = 1, length=8)
                    ax.tick_params(axis='x', width = 1, length=8)
                    plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)
                    
                    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)  # Moves the legend box to the right
                    
                    plt.title(f'{primiRNA}, mutation found in {str(frequency_mutation_rounded)} % of samples')
                    plt.ylim(-5,105)
                    plt.yticks([0,25,50,75,100])
                    
                    plt.ylabel('Percentage of DICER cutting site')
                    plt.xlabel('Position on pri-miRNA (30nt extension)')
                    # plt.xticks(visible=False,size=8)
                    # plt.yticks(visible=False)
                    mutation = mutation.replace('>','_to_')
                    path_save_fig = ''
                    plt.savefig(path_save_fig+f'/{primiRNA}_{mutation}_WT_vs_mutation_isomiR_frequency.png', dpi=150, bbox_inches='tight')
                    #plt.show()
                    plt.close()

        






















