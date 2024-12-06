#!/usr/bin/env python
# coding: utf-8

# In[19]:


from random import randint
import numpy as np


# In[20]:


bound_seqs = []
with open("boundcentered.fasta", 'r') as file:
    to_add = ""
    for line in file:
        if '>' not in line:
            line = line.replace('\n', '')
            to_add += line
            if len(line) == 1:
                bound_seqs.append(to_add)
                to_add = ""
bound_seqs[0]


# In[49]:


nuc_to_num_dict = {'A':0, 'C':1, 'G':2, 'T':3}

def rand_kmer_selection(seqs, k, t):
    motifs = []
    max_pos = len(seqs[0]) - k + 1
    for i in range(t):
        start_pos = randint(0, max_pos)
        motifs.append(seqs[i][start_pos:start_pos+k])
    return motifs

def build_count_matrix(motifs, exclude, k):
    count_matrix = np.zeros((4, k))
    for i in range(len(motifs)):
        if i != exclude:
            for j in range(len(motifs[i])):
                nuc = motifs[i][j]
                row = nuc_to_num_dict[nuc]
                count_matrix[row][j] += 1
    return count_matrix

def include_pseudocounts(count_matrix):
    for i in range(len(count_matrix)):
        for j in range(len(count_matrix[i])):
            count_matrix[i][j] += 1
    return count_matrix

def build_profile_matrix(count_matrix, k):
    divisor = np.sum(count_matrix, axis=0)[0] #sum of each column
    profile_matrix = np.zeros((4, k))
    for i in range(4):
        for j in range(k):
            to_add = (count_matrix[i][j] / divisor)
            profile_matrix[i][j] = to_add
            #print(profile_matrix[i][j])
    return profile_matrix

def best_kmer(profile_matrix, seq, k):
    best_score = 0
    best_kmer = ""
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        score = 0
        for j in range(len(kmer)):
            row = nuc_to_num_dict[kmer[j]]
            score += profile_matrix[row][j]
        if score > best_score:
            best_score = score
            best_kmer = kmer
    return best_kmer

def Score(motifs, profile_matrix):
    total_score = 0
    for motif in motifs:
        for i in range(len(motif)): #should be same as for i in range(k)
            row = nuc_to_num_dict[motif[i]]
            total_score += profile_matrix[row][i]
    return total_score

def GibbsSampler(seqs, k, t, N):
    motifs = rand_kmer_selection(seqs, k, t)
    best_motifs = motifs
    for i in range(N):
        exclude = randint(0, t-1)
        count_matrix = build_count_matrix(motifs, exclude, k)
        count_matrix = include_pseudocounts(count_matrix)
        profile_matrix = build_profile_matrix(count_matrix, k)
        motifs[exclude] = best_kmer(profile_matrix, seqs[exclude], k)
        if Score(motifs, profile_matrix) > Score(best_motifs, profile_matrix):
            best_motifs = motifs
    return best_motifs


# In[62]:


#len(bound_seqs2)


# In[82]:


bound_seqs2 = []
for seq in bound_seqs:
    bound_seqs2.append(seq[75:125])
best_motifs = GibbsSampler(bound_seqs2, 20, 200, 700)


# In[83]:


#build new profile_matrix using best_motifs
count_matrix = build_count_matrix(best_motifs, 0, 20)
count_matrix = include_pseudocounts(count_matrix)
profile_matrix = build_profile_matrix(count_matrix, 20)


# In[84]:


offset_seqs = []
with open("boundrandomoffset.fasta", 'r') as file:
    to_add = ""
    for line in file:
        if '>' not in line:
            line = line.replace('\n', '')
            to_add += line
            if len(line) == 1:
                offset_seqs.append(to_add)
                to_add = ""


# In[85]:


def best_kmer2(profile_matrix, seq, k):
    best_score = 0
    best_kmer_pos = 100
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        score = 0
        for j in range(len(kmer)):
            row = nuc_to_num_dict[kmer[j]]
            score += profile_matrix[row][j]
        if score > best_score:
            best_score = score
            best_kmer_pos = i
    return best_kmer_pos

for i in range(len(offset_seqs)):
    seq = offset_seqs[i]
    predicted_loc = best_kmer2(profile_matrix, seq, 20)
    with open('predictions.csv', 'a') as file:
        file.write("seq" + str(i+1) + " " + str(predicted_loc) + "\n")

