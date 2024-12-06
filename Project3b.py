#!/usr/bin/env python
# coding: utf-8

# In[43]:


from random import randint
import numpy as np


# In[29]:


bound_seqs = []
with open("bound.fasta", 'r') as file:
    to_add = ""
    for line in file:
        if '>' not in line:
            line = line.replace('\n', '')
            to_add += line
            if len(line) == 1:
                bound_seqs.append(to_add)
                to_add = ""


# In[30]:


test_seqs = []
with open("test.fasta", 'r') as file:
    to_add = ""
    for line in file:
        if '>' not in line:
            line = line.replace('\n', '')
            to_add += line
            if len(line) == 1:
                test_seqs.append(to_add)
                to_add = ""

unbound_seqs = []
with open("notbound.fasta", 'r') as file:
    to_add = ""
    for line in file:
        if '>' not in line:
            line = line.replace('\n', '')
            to_add += line
            if len(line) == 1:
                unbound_seqs.append(to_add)
                to_add = ""


# In[65]:


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
        #print(i)
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
    total_score = 1
    for motif in motifs:
        for i in range(len(motif)): #should be same as for i in range(k)
            row = nuc_to_num_dict[motif[i]]
            total_score += profile_matrix[row][i]
    return total_score

def GibbsSampler(seqs, k, t, N):
    motifs = rand_kmer_selection(seqs, k, t)
    best_motifs = motifs.copy()
    for i in range(N):
        exclude = randint(0, t-1)
        count_matrix = build_count_matrix(motifs, exclude, k)
        count_matrix = include_pseudocounts(count_matrix)
        profile_matrix = build_profile_matrix(count_matrix, k)
        motifs[exclude] = best_kmer(profile_matrix, seqs[exclude], k)
        print(Score(motifs, profile_matrix), Score(best_motifs, profile_matrix))
        if motifs == best_motifs:
            print("breaking")
            break
        if Score(motifs, profile_matrix) > Score(best_motifs, profile_matrix):
                best_motifs = motifs.copy()
    return best_motifs


# In[112]:


bound_seqs2 = []
for seq in bound_seqs:
    bound_seqs2.append(seq[80:120])
best_motifs = GibbsSampler(bound_seqs2, 20, 20, 400)


# In[113]:


best_motifs.append('x')
print(len(best_motifs))


# In[114]:


#build new profile_matrix using best_motifs
count_matrix = build_count_matrix(best_motifs, 20, 20)
count_matrix = include_pseudocounts(count_matrix)
profile_matrix = build_profile_matrix(count_matrix, 20)


# In[101]:


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
    return best_score, best_kmer_pos


# In[115]:


seq_to_score = {}
for i in range(len(test_seqs)):
    seq = test_seqs[i]
    best_score, best_loc = best_kmer2(profile_matrix, seq, 20)
    seq_to_score["seq" + str(i)] = best_score
print(len(seq_to_score))


# In[116]:


#sort dictionary by scores and then take the top 6000
sorted_seq_scores = dict(sorted(seq_to_score.items(), key=lambda x: x[1], reverse=True))
iter = 0
for seq in sorted_seq_scores.keys():
    if iter < 6000:
        with open('predictions.csv', 'a') as file:
            file.write(seq + "\n")
    iter += 1

