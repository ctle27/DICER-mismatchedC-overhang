# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 14:35:08 2025

@author: congt
"""

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

            
import pandas as pd
df_subgroup = pd.DataFrame()
df_subgroup['Subgroup'] = sub_list
df_subgroup['shRNA_sequence'] = seq_list