#!/usr/bin/env python
# coding: utf-8

# In[240]:


#Six amino acid sequene matrix data is searched for multi-fasta file made by H Nakano 20220908.
# Calc argorithm is copied from target_scoreing_fast made by Soki Nakano.
# The version can use previously calculated score_ref  in 20220909
# numpy version should be higher than 0.24
# Bags are fixed by H Nakano 20230613
#score_ref_aa6_matrix.csv file ( big size file) is a kind of dictionary of sum of index to calculation much faster
# 
# You can make it in an line belwo :   new_score_ref = 'Y' # chose new score_reference Y or N
# Once calculated , it is save as score_ref_name, so you can use  by refering the same file next time.
# You can adjust the number of list by changing minimun_score

import pandas as pd
import openpyxl as px
import sys
import csv


# In[241]:


from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import generic_protein


# In[242]:


print(pd.__version__)


# In[243]:


new_score_ref = 'N' # chose new score_reference Y or N
score_ref_name = 'score_ref_aa6_matrix.csv'# name of pre-calc score
score_table = pd.read_csv('tg1.csv') #name of score tabel. 
fasta_in = 'Human_proteins.fasta' # name of  multifasta file for serach.
minimum_score = 9 # sepcify minimun score 
csv_name = 'mt6_output.csv' # output file namme of Q and its surrounding without Q close to the C-terminal

score_table


# In[244]:


score_table_list = score_table.sort_values("AA").values.tolist()


# In[245]:



seq_t=[]
seq_t_t=[]
outfile=[]

columns1 = ["id","desc","posi","score","target","sequence"]

outfile=pd.DataFrame(columns = columns1)
#outfile_1=pd.DataFrame(data=data1,columns = columns1)


# In[246]:


import numpy as np
from collections import OrderedDict
import re


# In[247]:


pep_num_list = [("A", 0), ("C", 1), ("D", 2), ("E", 3), ("F", 4), ("G", 5), ("H", 6),           ("I", 7), ("K", 8), ("L", 9), ("M", 10), ("N", 11), ("P", 12), ("Q", 13),          ("R", 14), ("S", 15), ("T", 16), ("V", 17), ("W", 18), ("Y", 19)]
pep_num_dict=dict(pep_num_list) 


# In[248]:



def pep2num(pep):
    num = 0
    for i in range(6):
        if(pep[i] not in pep_num_dict):
            raise Exception("{} is not in pep_num_dict".format(pep[i]))
        num += pep_num_dict[pep[i]] * 20 **(5-i) #10進法に変換 
    return num


# In[249]:



def num2pep(num):
    base_20 = np.base_repr(num, 20)   # 10進法を一般的20進法に変換
    out = ""
    for char in str(base_20).zfill(6):  # 文字列に変換　4桁に0埋めて１個の文字を取り出す
        base_10 = int(char, 20)         # 10進法に変換
        out += pep_num_list[base_10][0] #配列に変換
    return out


# In[250]:


def pep2score(pep):    #scoreを計算
    score = 0
    for i, item in enumerate(pep):
        num = pep_num_dict[item]
        score += score_table_list[num][i+1]
    return score


# In[251]:


score_ref = []
if new_score_ref == "Y":
    
    
    
    for i in range(20**6):
        pep = num2pep(i)
        score_ref.append([pep, pep2score(pep)])
    with open(score_ref_name, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(score_ref)
else:
    
    df_score = pd.read_csv(score_ref_name)
    score_ref = df_score.to_numpy().tolist()
    


# In[ ]:





# In[252]:


out = []

nn = 0
for record in SeqIO.parse(fasta_in, 'fasta'):
    id_part = record.id
    desc_part = record.description
    seq = record.seq
    #print('id:', id_part)
    #print('desc:', desc_part)
    seq_str = str(seq)
    
    search_result = []
   
 
  
    for i in range (len(seq_str)-5):
        search_result.append(seq_str[i:i+6])
       
    #print (search_result)
    k=0
    for item in search_result:
        
        k += 1
        pep = item
        #print (item)
        try:
            
            num = pep2num(item)
            score = score_ref[num][1]
        except Exception as e:
            
            #print(str(e))
            score = 0
            #print(score)
        if score >= minimum_score:

            seq_t = [id_part, desc_part, k, score, pep, str(seq)]
            #print (seq_s.sum())
            #print (seq_t)
            out.append(seq_t) 


# In[253]:


with open(csv_name, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(out)


# In[ ]:





# In[ ]:





# In[ ]:




