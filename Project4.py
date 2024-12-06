#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


input = ""
with open('input.fasta', 'r') as file:
    for line in file:
        if '>' not in line:
            line = line.replace('\n', '')
            input += line
print(len(input))


# In[126]:


def forward_backward(emission, transition_matrix, emission_matrix, states):
    #forward
    start_prob = 1
    Fwd_dict = {}
    for i in range(len(emission)):
        for state in states:
            if state not in Fwd_dict:
                Fwd_dict[state] = {}
            if i == 0:
                Fwd_dict[state][i] = 1 * start_prob * emission_matrix[state][emission[i]]
            else:
                Fwd_dict[state][i] = 0
                for state2 in states:
                    Fwd_dict[state][i] += Fwd_dict[state2][i-1] * transition_matrix[state2][state] * emission_matrix[state][emission[i]]
        #normalization to avoid divide by 0 error
        divisor = 0
        for state in states:
            divisor += Fwd_dict[state][i]
        for state in states:
            Fwd_dict[state][i] /= divisor
    
    #print(total)
    #backward
    Bkwd_dict = {}
    for i in range(len(emission) - 1, -1, -1):
        for state in states:
            if state not in Bkwd_dict:
                Bkwd_dict[state] = {}
            if i == len(emission) - 1:
                Bkwd_dict[state][i] = 1
            else:
                Bkwd_dict[state][i] = 0
                for state2 in states:
                    Bkwd_dict[state][i] += Bkwd_dict[state2][i+1] * transition_matrix[state][state2] * emission_matrix[state2][emission[i+1]]
        #Backward normalization to avoid divide by 0 error
        divisor = 0
        for state in states:
            divisor += Bkwd_dict[state][i]
        for state in states:
            Bkwd_dict[state][i] /= divisor
    return Fwd_dict, Bkwd_dict

def parameter_estimation(emission, Fwd_dict, Bkwd_dict, chars, states, tran_matrix, emm_matrix):
    #updating transition and emission probability matrix
    #transition:
    transition_matrix = np.zeros((2,2))
    pos_values = np.zeros((2,2))
    for i in range(len(emission)-1):
        divisor = 0
        for state in states:
            for state2 in states:
                value = Fwd_dict[state][i] * Bkwd_dict[state2][i+1] * tran_matrix[state][state2] * emm_matrix[state][emission[i]]#* .5 * .25
                #* tran_matrix[state][state2] * emm_matrix[state][emission[i]]
                #* .5*.25 #think i will have to change this to prev trans and emission vals
                #* tran_matrix[state][state2] * emm_matrix[state][emission[i]]
                divisor += value
                pos_values[state][state2] = value
        for state in states:
            for state2 in states:
                update = pos_values[state][state2] / divisor
                transition_matrix[state][state2] += update
                pos_values[state][state2] = 0
    for row in transition_matrix:
        norm = row[0] + row[1]
        row[0] /= norm
        row[1] /= norm
    #emission
    P_sum = 0
    notP_sum = 0
    emission_matrix = np.zeros((2, 4))
    count_matrix = np.zeros((2, 4))
    sum1 = [0,0]
    for i in range(len(emission)):
        sum1[0] += Fwd_dict[0][i]
        sum1[1] += Fwd_dict[1][i]
        count_matrix[0][emission[i]] += Fwd_dict[0][i]
        count_matrix[1][emission[i]] += Fwd_dict[1][i]
    for row in range(2):
        for item in range(4):
            emission_matrix[row][item] = count_matrix[row][item] / sum1[row]
    #print("transition: " , transition_matrix)
    #print("emission: ", emission_matrix)
    return transition_matrix, emission_matrix

def EM(emission, initial_transition, initial_emission, chars, states, num_iterations):
    transition_matrix = initial_transition
    emission_matrix = initial_emission
    for i in range(num_iterations):
        #forward-backward
        Fwd_dict, Bkwd_dict = forward_backward(emission, transition_matrix, emission_matrix, states)
        #update transition_emission
        transition_matrix, emission_matrix = parameter_estimation(emission, Fwd_dict, Bkwd_dict, chars, states, initial_transition, initial_emission)
        initial_transition = transition_matrix
        initial_emission = emission_matrix
        print("transition: ", transition_matrix)
        print("emission: ", emission_matrix)
    return transition_matrix, emission_matrix

def converter(input):
    convert = []
    letter_to_num_dict = {'x':0, 'y':1, 'z':2, 'n':3}
    for letter in input:
        convert.append(letter_to_num_dict[letter])   
    return convert

def print_transition(matrix):
    print("P", "~P")
    to_print = []
    j = 0
    for i in matrix:
        for k in matrix[i]: 
            j += 1
            to_print.append(matrix[i][k])
            if j == 2:
                j = 0
                
                print(round(to_print[0], 6), round(to_print[1], 6))#, round(to_print[3], 4))
                to_print = []

def print_emission(matrix):
    print("x", "y", "z", 'n')
    to_print = []
    j = 0
    for i in matrix:
        for k in matrix[i]: 
            j += 1
            to_print.append(matrix[i][k])
            if j == 4:
                j = 0
                print(round(to_print[0], 6), round(to_print[1], 6), round(to_print[2], 6), round(to_print[3],6))
                to_print = []

def state_estimation(fwd_dict, bkwd_dict, emission):
    pos_to_prob_dict = {}
    for i in range(len(emission)):
        value = (fwd_dict[0][i] * bkwd_dict[0][i]) / (fwd_dict[0][i] * bkwd_dict[0][i] + fwd_dict[1][i] * bkwd_dict[1][i])
        pos_to_prob_dict[i] = value
    return pos_to_prob_dict
                


# In[127]:


emission = converter(input)
num_iterations = 10
chars = [0, 1, 2, 3]
states = [0, 1]
initial_transition = [[0.9, 0.1],
                     [0.3, 0.7]]
initial_emission = [[1/3,	1/3,	1/3, 0.1],
                    [0.1,	0.1,	0.1, 0.7]]
transition_matrix, emission_matrix = EM(emission, initial_transition, initial_emission, chars, states, num_iterations)
print(transition_matrix)
print(emission_matrix)


# In[128]:


fwd_dict, bkwd_dict = forward_backward(emission, transition_matrix, emission_matrix, states)


# In[129]:


pos_to_prob_dict = state_estimation(fwd_dict, bkwd_dict, emission)
sorted_seq_scores = dict(sorted(pos_to_prob_dict.items(), key=lambda x: x[1], reverse=True))
iter = 0
for seq in sorted_seq_scores.keys():
    if iter < 50000:
        with open('predictions.csv', 'a') as file:
            file.write(str(seq+1) + "\n")
    iter += 1

