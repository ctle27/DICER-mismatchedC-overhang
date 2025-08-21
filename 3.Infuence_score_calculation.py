# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 20:48:38 2024

@author: congt
"""

path = 'D:/HKUST_Research/Dicer-TLR16-27/ngs_analysis/'
import pandas as pd
processed_df_wt = pd.read_csv(path+'processed_df_wt_tlr1627_2reps.bed',sep='\t')
df = processed_df_wt.copy()

df = df[['shRNA_sequence', 'New_define_structure_1', 'New_define_structure_2', 'concrete_struct', '5p-3p', '5p-3p-alternative', 'Mean_Cleavage_accuracy',
                 'Mean_Cleavage_accuracy_of_alternative_5p_3p']]
df.rename(columns={
    'New_define_structure_1': 'new_define_struct1',
    'New_define_structure_2': 'new_define_struct2'},inplace=True)
#%%combine with data from TLR10-15 and TLR1-9
import pandas as pd
'''
TLR10-15
'''
df1 = pd.read_csv(r'D:\HKUST_Research\Dicer-two-motifs\ngs_analysis\processed_df_wt_trl1015_3reps.bed', sep='\t')
df1_dc = df1[df1['SC_on'] == 'None']

df1_dc = df1_dc[['shRNA_sequence', 'new_define_struct1', 'new_define_struct2', 'concrete_struct', '5p-3p', '5p-3p-alternative', 'Mean_Cleavage_accuracy',
                 'Mean_Cleavage_accuracy_of_alternative_5p_3p']]

'''
TLR1-9
'''

path2 = 'D:/HKUST_Research/Dicer-loop-counting-redefinition/ngs_analysis/TLR1-9-from-TamAnh/Truc-processed/'
df2 = pd.read_csv(path2+'processed_df_wt_tlr19.bed',sep='\t')
df2_dc = df2[df2['Type-general'] == 'DC']
df2_dc = df2_dc.copy()
df2_dc.rename(columns={'Reference sequence': 'shRNA_sequence'}, inplace=True)
df2_dc = df2_dc[['shRNA_sequence', 'new_define_struct1', 'new_define_struct2', 'concrete_struct', '5p-3p', '5p-3p-alternative', 'Mean_Cleavage_accuracy',
                 'Mean_Cleavage_accuracy_of_alternative_5p_3p']]

#combine TLR16-27 with TLR1-9 and TLR10-15
'''
For shRNA variants that appear in more than 1 library, calculate average cleavage accuracy at each cleavage site
'''

df_combine = pd.concat([df, df1_dc, df2_dc], axis=0)
df_combine['Mean_Cleavage_accuracy_3libs'] = df_combine.groupby(['shRNA_sequence','5p-3p'])['Mean_Cleavage_accuracy'].transform('mean')
df_combine['Mean_Cleavage_accuracy_of_alternative_5p_3p_3libs'] = df_combine.groupby(['shRNA_sequence','5p-3p-alternative'])['Mean_Cleavage_accuracy_of_alternative_5p_3p'].transform('mean')

df_combine.drop_duplicates(subset=['shRNA_sequence','5p-3p'], keep='first',inplace=True)
df_combine.drop(['Mean_Cleavage_accuracy','Mean_Cleavage_accuracy_of_alternative_5p_3p'], axis=1, inplace=True)
df_combine.rename(columns={
    'Mean_Cleavage_accuracy_of_alternative_5p_3p_3libs': 'Mean_Cleavage_accuracy_of_alternative_5p_3p',
    'Mean_Cleavage_accuracy_3libs': 'Mean_Cleavage_accuracy'
}, inplace=True)

#save merged dataframe
df_combine.to_csv(path+'processed_dt_wt_merge_TLR1-9_TLR10-15_TLR16-27.bed', sep='\t', index=False)
#%%open dataframe to work
#if only want to check TLR16-27
path = 'D:/HKUST_Research/Dicer-TLR16-27/ngs_analysis/'
import pandas as pd
processed_df_wt = pd.read_csv(path+'processed_df_wt_rep1.bed',sep='\t')
df = processed_df_wt.copy()
df = df[(df['Stem'] == '23L') & (df['Symmetric_structure'] == 'Yes') & (df['5p-3p'].isin(['21-21','22-22']))]# & (df['5p_flanking_length'] == 0)]
df = df[['Variant', 'shRNA_sequence', '5p-3p', 'Mean_Cleavage_accuracy']]
#%%open dataframe to work
#if want to check merged dataframe from 3 libs
path = 'D:/HKUST_Research/Dicer-TLR16-27/ngs_analysis/'
import pandas as pd
processed_df_wt = pd.read_csv(path+'processed_dt_wt_merge_TLR1-9_TLR10-15_TLR16-27.bed',sep='\t')
df = processed_df_wt.copy()
strc_list = list(set(df['concrete_struct'].tolist()))
symmetric_strc = [strc for strc in strc_list if 'A' not in strc and 'B' not in strc and '23-SSSSSSS SSSSSSS-23' in strc]

df = df[(df['concrete_struct'].isin(symmetric_strc)) & (df['5p-3p'].isin(['21-21','22-22']))]

df.reset_index(inplace=True, drop=True)

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
    df = df[['shRNA_sequence', '5p-3p', metrics]]
    df.drop_duplicates(subset=['shRNA_sequence','5p-3p'], keep='first', inplace=True)
    df = df.pivot(index=['shRNA_sequence'], columns='5p-3p', values=metrics)
    df.fillna(fillin_value, inplace=True)
    df.reset_index(inplace=True)
    df = pd.melt(df, id_vars=['shRNA_sequence'], value_vars=['21-21', '22-22'],
                        var_name='5p-3p', value_name=metrics)
    df.sort_values(['shRNA_sequence', '5p-3p'], ascending=[True, True], inplace=True)
    df.reset_index(inplace=True, drop=True)
    
    return (df)

df = adding_missing_cleavage_site(df, 'Mean_Cleavage_accuracy', 0)



#%%check DC21, DC22 for each three-pair at each randomized position
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
position_1 = []

for i in range(1,21):
    for comb1 in tri_nu_comb:
        for comb2 in tri_nu_comb:
            seq = ini_seq[:i] + comb1 + ini_seq[i+3:99-i] + comb2 + ini_seq[102-i:]
            seq_list.append(seq.replace('N'*32,''))
            sub_list.append((i))
            tri5p.append(comb1)
            tri3p.append(comb2)
            position_1.append(comb1 + '-' + comb2[::-1])

df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
df_subgroup['Subgroup'] = sub_list
df_subgroup['shRNA_sequence'] = seq_list
df_subgroup['tri5p'] = tri5p
df_subgroup['tri3p'] = tri3p
df_subgroup['Position 1'] = position_1

df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner")

df1 = pd.melt(df_subgroup, id_vars=['5p-3p', 'Mean_Cleavage_accuracy', 'Subgroup'],
                  value_vars = ['Position 1'], var_name='Position', value_name='Nt_combination')

df1.reset_index(inplace=True, drop=True)

for i,pos in enumerate(df1['Position']):
    sub = int(df1['Subgroup'][i])
    pos_split = pos.split(' ')
    df1.loc[i,'Randomized_position'] = int(pos_split[1]) + sub

df1 = df1[['Nt_combination','Randomized_position', '5p-3p', 'Mean_Cleavage_accuracy']]

def extract_clv_site(df_input, clv_site):
    df = df_input.copy()
    df = df[df['5p-3p'] == clv_site]
    df.rename(columns={'Mean_Cleavage_accuracy':'Mean_Cleavage_accuracy_at_'+clv_site,
                       'Randomized_position': 'Randomized_position_checking_'+clv_site}, inplace=True)
    del df['5p-3p']

    df.reset_index(inplace=True, drop=True)
    return (df)

df_21 = extract_clv_site(df1, '21-21')
df_22 = extract_clv_site(df1, '22-22')

df_21['Shifting_forwarding_position'] = df_21['Randomized_position_checking_21-21'] + 1

df_merge = df_21.merge(df_22, left_on=['Nt_combination', 'Shifting_forwarding_position'], right_on=['Nt_combination','Randomized_position_checking_22-22'], how='inner')

df_merge['Mean_DC21_DC22_accuracy'] = df_merge[['Mean_Cleavage_accuracy_at_21-21', 'Mean_Cleavage_accuracy_at_22-22']].mean(axis=1)

df_merge.sort_values(['Randomized_position_checking_21-21', 'Mean_DC21_DC22_accuracy'], ascending=[True,False], inplace=True)
df_merge.reset_index(inplace=True, drop=True)

df_merge['Randomized_position_checking_21-21'] = df_merge['Randomized_position_checking_21-21'].astype(int)

# df_merge.to_csv(path+'shifting_score_all_position.bed', sep='\t',index=False)
#%%plotting nt frequency revised (checking whole pool frequency & bottom group)
df1 = df_merge.copy()
#average_shifting_score of whole pool
average_shifting_score = df1['Mean_DC21_DC22_accuracy'].mean()

df1['1st_pair'] = df1['Nt_combination'].str[0] + '-' + df1['Nt_combination'].str[4]
df1['2nd_pair'] = df1['Nt_combination'].str[1] + '-' + df1['Nt_combination'].str[5]
df1['3rd_pair'] = df1['Nt_combination'].str[2] + '-' + df1['Nt_combination'].str[6]

#checking position: starting position of the combination
for pos in range(2,20):
    df = df1[df1['Randomized_position_checking_21-21'] == pos]
    df = df.copy()
    df.sort_values(['Randomized_position_checking_21-21', 'Mean_DC21_DC22_accuracy'], ascending=[True, False], inplace=True)
    
    '''
    1) Substract cleavage shifting score of each variant by the average shifting score of the whole pool (all subgroups) --> re_scaled score
    2) Select top and bottom variants with highest/lowest re_scaled score --> average_re_scaled_score
    3) Calculate average cleavage shifting score for each nt in each position on each strand for top and bottom variants.
    4) Calculate frequency of each nt in each position on each strand in the whole pool and top/bottom variants
    5) Divide frequency of each nt in top/bottom group to the frequency in whole window --> re_scaled_frequency
    3) Multiple two values. Then calculate weighted proportion: (average_re_scaled_score * re_scaled_frequency) / sum(average_re_scaled_score * re_scaled_frequency)
    '''
    
    #Step 1: calculate rescaled score for all variants
    def calculate_rescaled_score(df_input, average_shifting_score):
        df = df_input.copy()
        average_shifting_score = average_shifting_score
        # average_shifting_score = df['Mean_DC21_DC22_accuracy'].mean()
        df['Re_scaled_score'] = df['Mean_DC21_DC22_accuracy'] - average_shifting_score
        return df
    
    df = calculate_rescaled_score(df, average_shifting_score)    
    
    #select top/bottom few percents with highest mean accuracy score
    n = int(len(df) * 0.03)
    if n > 30:
        df_top = df.nlargest(n, 'Re_scaled_score')
        df_bottom = df.nsmallest(n, 'Re_scaled_score')
    elif n <= 30:
        df_top = df.nlargest(30, 'Re_scaled_score')
        df_bottom = df.nsmallest(30, 'Re_scaled_score')
    
    #Step 2: calculate average rescaled score for each nt at each position from top/bottom groups
        
    def calculate_average_rescaled_score(df_input, position):
        df = df_input.copy()
        
        all_pairs = []
        for nu1 in ['A','T','G','C']:
            for nu2 in ['A','T','G','C']:
                all_pairs.append(nu1 + '-' + nu2)
                
        df_result = df.groupby(position + '_pair')['Re_scaled_score'].mean().reset_index()
        df_result.rename(columns={position + '_pair': 'Pair', 'Re_scaled_score': position + '_pair'}, inplace=True)
        
        '''
        in top/bottom variants, there might be no A/T/C/G at certain position --> bug in later merging steps.
        Need to add missing pair with rescaled_score = 0
        '''
        existing_pairs = df_result['Pair'].unique()
        missing_pairs = list(set(all_pairs) - set(existing_pairs))
        missing_df = pd.DataFrame({'Pair': missing_pairs, position + '_pair': 0})
        
        df_result = pd.concat([df_result, missing_df], ignore_index=True)
        return df_result
    
    def merge_dataframes(df1, df2, df3):
        df = df1.merge(df2, on=['Pair'], how='inner')
        df = df.merge(df3, on=['Pair'], how='inner')
        df.set_index('Pair', inplace=True)
        df.sort_index(inplace=True)
        return df
    
    df_top_average_rescaled_score = merge_dataframes(calculate_average_rescaled_score(df_top, '1st'), calculate_average_rescaled_score(df_top, '2nd'),
                                                        calculate_average_rescaled_score(df_top, '3rd'))
    df_bottom_average_rescaled_score = merge_dataframes(calculate_average_rescaled_score(df_bottom, '1st'), calculate_average_rescaled_score(df_bottom, '2nd'),
                                                           calculate_average_rescaled_score(df_bottom, '3rd'))
    
    #Step 3: calculate the frequency of each bp at each randomized position on each strand
    '''
    1. calculate frequency of pairs in whole pool and in top/bottom groups separately
    2. In whole pool, if the number of variants containing the checking pair is too low (<10), set the count to 0.
    3. To be fair with the frequency in top/bottom group, the frequency of pair satisfied variant count (>10 variants) is calculated by dividing the pair count by total count of pair
    at a position in the pool (before setting the count of low-count pair to 0)
    '''
    def calculate_frequency_whole_pool(df_input):
        df = df_input.copy()
        df = df[['1st_pair', '2nd_pair', '3rd_pair']]
        frequency_df = df.apply(lambda x: pd.Series(list(x)).value_counts()).fillna(0).astype(int)
        '''
        store sum count of pairs at each position to a list
        '''
        sum_count_pair = frequency_df.sum().tolist()
        '''
        set low-count (<10) to 0
        '''
        frequency_df[frequency_df < 10] = 0
        frequency_df = frequency_df.div(sum_count_pair)

        return frequency_df
    
    df_whole_pool_frequency = calculate_frequency_whole_pool(df)
    if len(df_whole_pool_frequency.index) < 16:
        print ('Warning: not enough 16 pairs in one of the position in the pool')
        
    def calculate_frequency_top_bottom_variants(df_input):
        df = df_input.copy()
        df = df[['1st_pair', '2nd_pair', '3rd_pair']]
        frequency_df = df.apply(lambda x: pd.Series(list(x)).value_counts()).fillna(0).astype(int)
        frequency_df = frequency_df.div(frequency_df.sum())
        frequency_df.sort_index(inplace=True)
        return frequency_df
    
    df_top_frequency = calculate_frequency_top_bottom_variants(df_top)
    df_bottom_frequency = calculate_frequency_top_bottom_variants(df_bottom)

    #Step 4: calculate rescaled frequency of each bp at each randomized position on each strand
    '''
    IMPORTANT NOTE:
        1) In whole pool frequency, low-count pairs was set to 0 count --> frequency of those pairs is 0 in whole pool.
        Those low-count pairs might appear in top/bottom group. Hence, when calculating rescaled frequency by dividing frequency in top/bottom group by frequency
        in whole pool, results return infinite (inf). Need to replace by 0.
        
        2) If a pair not appear in top/bottom group, when divide frequency in top/bottom by frequency in whole pool, results return NaN (even when the frequency in whole pool is 0)
    Replace NaN by 0.
    
    The idea is that is a pair does not appear in top or bottom or its count in whole pool is too low, the impact to shifting activity is 0.
    '''
    def calculate_rescaled_frequency(df_input1, df_input2):
        import numpy as np
        df1 = df_input1.copy() #frequency in top/bottom variant group
        df2 = df_input2.copy() #frequency in whole randomization pool
        df = df1.div(df2)
        df.fillna(0, inplace=True)
        df.replace([np.inf, -np.inf], 0, inplace=True)
        return df
    
    df_top_rescaled_frequency = calculate_rescaled_frequency(df_top_frequency, df_whole_pool_frequency)
    df_bottom_rescaled_frequency = calculate_rescaled_frequency(df_bottom_frequency, df_whole_pool_frequency)
    
    #Step 5: calculate weighted shifting score by multiplying rescaled frequency with rescaled shifting score, then combine top and bottom groups
    def calculate_weighted_proportion(df_input1, df_input2):
        df1 = df_input1.copy() #rescaled_frequency
        df2 = df_input2.copy() #resclaed_score
        df1.fillna(0, inplace=True)
        df2.fillna(0, inplace=True)
        df = df1 * df2
        return df
    
    def combine_top_bottom(df_top, df_bottom):
        df1 = df_top.copy()
        df2 = df_bottom.copy()
        
        df = df1 + df2
        return df
    
    df_weighted = combine_top_bottom(calculate_weighted_proportion(df_top_rescaled_frequency, df_top_average_rescaled_score),
                                        calculate_weighted_proportion(df_bottom_rescaled_frequency, df_bottom_average_rescaled_score))


    #Step 6: plotting
    df = df_weighted.copy()
    df.sort_values(['1st_pair', '2nd_pair', '3rd_pair'], ascending=[True, False, False], inplace=True)
    
    df_1st = df[['1st_pair']]
    df_2nd = df[['2nd_pair']]
    df_3rd = df[['3rd_pair']]
    
    def extract_strand(df_input, pair, strand):
        df = df_input.copy()
        
        #sort values: >0 then ascending, <0 then descending --> expected pattern of the stacks in the bar plot
        neg_val = df[df[pair] < 0].sort_values(pair, ascending=False)
        pos_val = df[df[pair] >= 0].sort_values(pair, ascending=True)
        df = pd.concat([neg_val, pos_val])
        df.reset_index(inplace=True)
        
        if strand == '5p':
            df[strand] = df['index'].apply(lambda x: x.split('-')[0])
        if strand == '3p':
            df[strand] = df['index'].apply(lambda x: x.split('-')[1])
        
        #order of appearance, for the purpose of coloring in the barplot:
            # A, T, C, G --> A1, A2, A3, A4, T1, T2, T3, T4, etc depends on the order of appearance in 16 pairs after sorting
        df[strand+'_order'] = df.groupby(strand).cumcount() + 1
        df[strand] = df[strand] + df[strand+'_order'].astype(str)
        
        df1 = df[[strand, pair]]
        df1.set_index(strand, inplace=True)
    
        return df1
    df_1st_5p = extract_strand(df_1st, '1st_pair', '5p').T
    df_1st_3p = extract_strand(df_1st, '1st_pair', '3p').T
    df_2nd_5p = extract_strand(df_2nd, '2nd_pair', '5p').T
    df_2nd_3p = extract_strand(df_2nd, '2nd_pair', '3p').T
    df_3rd_5p = extract_strand(df_3rd, '3rd_pair', '5p').T
    df_3rd_3p = extract_strand(df_3rd, '3rd_pair', '3p').T
    
    import matplotlib.pyplot as plt
    
    colors = {'T1': 'r', 'A1': 'lightcyan', 'C1': 'lightsalmon', 'G1': 'dodgerblue',
              'T2': 'r', 'A2': 'lightcyan', 'C2': 'lightsalmon', 'G2': 'dodgerblue',
              'T3': 'r', 'A3': 'lightcyan', 'C3': 'lightsalmon', 'G3': 'dodgerblue',
              'T4': 'r', 'A4': 'lightcyan', 'C4': 'lightsalmon', 'G4': 'dodgerblue'}
    
    # Set up the figure with two subplots side by side
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 6, figsize=(2,3))
    
    
    df_1st_5p.plot.bar(stacked=True, color=[colors[col] for col in df_1st_5p.columns], ax=ax1, width=1, linewidth=0.2, edgecolor='black')
    df_1st_3p.plot.bar(stacked=True, color=[colors[col] for col in df_1st_3p.columns], ax=ax2, width=1, linewidth=0.2, edgecolor='black')
    
    df_2nd_5p.plot.bar(stacked=True, color=[colors[col] for col in df_2nd_5p.columns], ax=ax3, width=1, linewidth=0.2, edgecolor='black')
    df_2nd_3p.plot.bar(stacked=True, color=[colors[col] for col in df_2nd_3p.columns], ax=ax4, width=1, linewidth=0.2, edgecolor='black')
    
    df_3rd_5p.plot.bar(stacked=True, color=[colors[col] for col in df_3rd_5p.columns], ax=ax5, width=1, linewidth=0.2, edgecolor='black')
    df_3rd_3p.plot.bar(stacked=True, color=[colors[col] for col in df_3rd_3p.columns], ax=ax6, width=1, linewidth=0.2, edgecolor='black')
    
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.axhline(y=0, color='black', linestyle='dashed', linewidth=1)
        ax.set_ylim(-8, 8)
        ax.set_yticks([-8,-4,0,4,8])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.tick_params(axis='y', width=1, length=5)
        ax.tick_params(axis='x', width=1, length=5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.get_legend().remove()
        if ax in [ax2, ax3, ax4, ax5, ax6]:
            ax.set_yticks([])
            ax.spines['left'].set_visible(False)
        if ax in [ax1, ax3, ax5]:
            ax.set_xlim(-0.75, 0.5)
            ax.spines['right'].set_visible(False)
        if ax in [ax2, ax4, ax6]:
            ax.set_xlim(-0.5, 0.75)
            
    # Adjust the spacing between subplots
    plt.subplots_adjust(wspace=0)
    
    path = 'D:/HKUST_Research/Dicer-TLR16-27/ngs_analysis/'
    plt.savefig(path+f'bp_frequency_to_shifting_score_each_position/bp_no_{pos}_weighted_shifting_score_only_for_symmetric_23L_variants.png', dpi=300, bbox_inches='tight')

#%%4096x21 heatmap for shifting score
tri_nu_comb = []
for nu1 in ['A','T','G','C']:
     for nu2 in ['A','T','G','C']:
         for nu3 in ['A','T','G','C']:
             tri_nu_comb.append(nu1+nu2+nu3)
nt_combination = []

for comb1 in tri_nu_comb:
    for comb2 in tri_nu_comb:
        nt_combination.append(comb1 + '-' + comb2[::-1])

df_plot = pd.DataFrame(columns=range(2, 20), index=nt_combination)
for i,comb in enumerate(df_merge['Nt_combination']):
    pos = df_merge['Randomized_position_checking_21-21'][i]
    val = df_merge['Mean_DC21_DC22_accuracy'][i]
    df_plot.loc[comb, pos] = val

mean_shifting_score = df_merge['Mean_DC21_DC22_accuracy'].mean()

df_plot = df_plot - mean_shifting_score
column_means = df_plot.mean()
df_plot = df_plot.fillna(-999)

mask = df_plot == -999

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
size = (10,5)
ax = plt.figure(figsize=size) #(figsize=(3.75,3)) for SC, (figsize=(2.25,3)) for DC
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

ax = sns.heatmap(data = df_plot, mask = mask, cmap = 'coolwarm', cbar=False, vmax=0.6, vmin=-0.6)
ax.tick_params(axis='y', width = 0, length=0)
ax.tick_params(axis='x', width = 0, length=0)
#set color of nan cells (that were replaced by -999 value)
ax.set_facecolor('xkcd:silver')
for _, spine in ax.spines.items():
    spine.set_visible(True)
# Add vertical lines for each column
for col_index in range(df_plot.shape[1] + 1):
    ax.axvline(x=col_index, color='black', linewidth=0.5)
plt.xlabel('')
plt.ylabel('')
plt.yticks([])
plt.xticks([])

plt.savefig(path+'shifting_score_across_position_human_DICER.png', dpi=300, bbox_inches='tight')
plt.show()
#%%plot boxplot for differentiated shifting score
tri_nu_comb = []
for nu1 in ['A','T','G','C']:
     for nu2 in ['A','T','G','C']:
         for nu3 in ['A','T','G','C']:
             tri_nu_comb.append(nu1+nu2+nu3)
nt_combination = []

for comb1 in tri_nu_comb:
    for comb2 in tri_nu_comb:
        nt_combination.append(comb1 + '-' + comb2[::-1])

df_plot = pd.DataFrame(columns=range(2, 20), index=nt_combination)
for i,comb in enumerate(df_merge['Nt_combination']):
    pos = df_merge['Randomized_position_checking_21-21'][i]
    val = df_merge['Mean_DC21_DC22_accuracy'][i]
    df_plot.loc[comb, pos] = val

mean_shifting_score = df_merge['Mean_DC21_DC22_accuracy'].mean()

df_plot = df_plot - mean_shifting_score
df_plot.reset_index(inplace=True)
df1 = df_plot.melt(id_vars='index', value_vars=[i for i in range(2, 20)], var_name='Position', value_name='Differentaited_shifting_score')
df1.dropna(inplace=True)

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
size = (10,2)
ax = plt.figure(figsize=size) #(figsize=(3.75,3)) for SC, (figsize=(2.25,3)) for DC
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

# ax = sns.boxplot(data = df1, x = 'Position', y='Differentaited_shifting_score', showfliers=False, linewidth=1, color='dodgerblue', zorder=10)
ax = sns.stripplot(data = df1, x = 'Position', y='Differentaited_shifting_score', color='lightblue', edgecolor=".2", zorder=1, size=1)

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)

# ax.get_legend().remove()

plt.ylim(-0.6,0.6)
plt.yticks([-0.5,-0.25,0,0.25,0.5])

plt.ylabel('')
plt.xlabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)


plt.savefig(path+'shifting_score_across_position_human_DICER_stripplot.png', dpi=300, bbox_inches='tight')
plt.show()

#%%logomaker
'''
plot information bit for nt in each strand of the motif
'''
df1 = df_merge.copy()
df1['Nt_combination'] = df1['Nt_combination'].str.replace('T','U')
starting_randomization_position = 9
effect = 'enhancing'
df1 = df1[df1['Randomized_position_checking_21-21'] == starting_randomization_position]
n = int(len(df1) * 0.03)

if effect == 'inhibiting':
    if n > 30:
        df_plot = df1.nsmallest(n, 'Mean_DC21_DC22_accuracy')
    if n <= 30:
        df_plot = df1.nsmallest(30, 'Mean_DC21_DC22_accuracy')
if effect == 'enhancing':
    if n > 30:
        df_plot = df1.nlargest(n, 'Mean_DC21_DC22_accuracy')
    if n <= 30:
        df_plot = df1.nlargest(30, 'Mean_DC21_DC22_accuracy')
        
df_plot["5'-arm (from 5' to 3')"] = df_plot['Nt_combination'].str[:3]
df_plot["3'-arm (from 3' to 5')"] = df_plot['Nt_combination'].str[-3:]

import logomaker as lm
import matplotlib.pyplot as plt

def draw_weblogo (list_sequence, fig_name, save_fig = 'no'):
    path = 'D:/HKUST_Research/Dicer-TLR16-27/ngs_analysis/'
    # counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat['correct_index'] = counts_mat.index.map(lambda x: x+1)
    counts_mat = counts_mat.set_index('correct_index')

    crp_logo  = lm.Logo(counts_mat, 
                        figsize = [0.4*len(counts_mat), 1.2], 
                        color_scheme = {'A': 'lightcyan', 'C': 'lightsalmon',  'G': 'dodgerblue', 'U': 'r', 'T': 'r'},
                        font_name='Arial Rounded MT Bold', zorder = 3)
    
    for _, spine in crp_logo.ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_color('black')
    
    plt.yticks([0,0.75,1.5], fontsize = 0, color = 'white')
    plt.xticks([1,2,3], fontsize = 0, color = 'white')
    plt.tick_params(axis='y', width = 1, length=6)
    plt.tick_params(axis='x', width = 1, length=6)
    plt.grid(axis = 'both', color = 'lightgrey', linestyle = '--', linewidth = 0.5)
    
    if save_fig != 'no':
        # plt.title(save_fig)
        # plt.show()
        plt.savefig(path+f'{fig_name}.png'.format(save_fig), bbox_inches="tight", dpi =1000)
    else:
        plt.show()
    print (counts_mat)
    return()

draw_weblogo(df_plot["5'-arm (from 5' to 3')"].tolist(),  str(starting_randomization_position) + '-' + effect + '_motif_5p_arm', save_fig = 'yes')
draw_weblogo(df_plot["3'-arm (from 3' to 5')"].tolist(),  str(starting_randomization_position) + '-' + effect + '_motif_3p_arm',  save_fig = 'yes')

#%%
df_check = df_merge[df_merge['Randomized_position_checking_21-21'] == 10]
df_check = df_check.copy()

df_check['1st_pair'] = df_check['Nt_combination'].str[0] + '-' + df_check['Nt_combination'].str[4]
df_check['2nd_pair'] = df_check['Nt_combination'].str[1] + '-' + df_check['Nt_combination'].str[5]
df_check['3rd_pair'] = df_check['Nt_combination'].str[2] + '-' + df_check['Nt_combination'].str[6]

# df_check.sort_values(['Mean_Cleavage_accuracy_at_21-21'], ascending=False, inplace=True)
n = int(len(df_check) * 0.03)
if n > 30:
    df_top = df_check.nlargest(n, 'Mean_DC21_DC22_accuracy')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

checking = '1st_pair'
frequency = df_top[checking].value_counts()
print (frequency)
#%%
# df_check = df_merge[(df_merge['Randomized_position_checking_21-21'] == 9) & (df_merge['Nt_combination'].isin(['CAA-GTT','CAA-CTT']))]
df_check = df_merge[(df_merge['Randomized_position_checking_21-21'] == 11) & (df_merge['Nt_combination'].isin(['CTG-CAC','ATG-TAC']))]
# df_check = df_subgroup.copy()
# df_check['shRNA_sequence'] = df_check['shRNA_sequence'].str.replace('AAGCTTGCGCAAGCAA','AAAAAA')
# df_choosen = df_check[df_check['shRNA_sequence'] == 'GGGATATTT CCTCCAGATCAAGAAAAAAAAAAAAAAAACTTGATCTGGACG AAATATTCTTA'.replace(' ','')]


























