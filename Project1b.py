reference_genome = ""
with open("project1b-u_reference_genome.fasta" ,'r') as file:
    for line in file:
          if '>' not in line:
              line = line.replace('\n', '')
              reference_genome += line

reads = []
with open("project1b-u_with_error_paired_reads.fasta" ,'r') as file:
    for line in file:
          if '>' not in line:
              line = line.replace('\n', '')
              reads.append(line)

ref_map = {}
for i in range(len(reference_genome)-16):
    to_add = reference_genome[i:i+16]
    if to_add not in ref_map:
        ref_map[reference_genome[i:i+16]] = [i]
    else:
        ref_map[reference_genome[i:i+16]].append(i)

def find_substitutions(ref_kmer, read, index):
    count = 0
    mutation_index = []
    i=0
    for i in range(len(ref_kmer)):
        if ref_kmer[i] != read[i]:
            count += 1
            mutation_index.append(str(i+index+1)+ read[i])
        if count > 2:
           return [-1]
    return mutation_index


substitutions_list = []
for read in reads:
    read1 = read[0:16]
    read2 = read[16:32]
    read3 = read[32:48]
    

    #2) check for substitutions
    #if two of the reads match but not the third
    if (read1 in ref_map and read2 in ref_map and read3 not in ref_map) or (read1 in ref_map and read2 not in ref_map and read3 in ref_map):
        #look at each location that read1 maps to and run the trivial algorithm
        for location in ref_map[read1]:
            mutation_index = find_substitutions(reference_genome[location:location+48], read[:48], location)
            for index in mutation_index:
                if index != -1:
                    substitutions_list.append(index)
    if read1 not in ref_map and read2 in ref_map and read3 in ref_map:
        for location in ref_map[read2]:
            mutation_index = find_substitutions(reference_genome[location-16:location+32], read[:48], location-16)
            for index in mutation_index:
                if index != -1:
                    substitutions_list.append(index) 


sorted_list = sorted(substitutions_list)
my_dict = {}
for i in range(len(sorted_list)-1):
    if sorted_list[i] == sorted_list[i+1]:
        nuc = sorted_list[i][-1]
        loc = sorted_list[i][:-1]
        my_dict[int(loc)] = reference_genome[int(loc)-1] + " " + nuc


def insertion_finder(ref_kmer, read, index):
    count = 0
    insertion = -1
    nuc = ''
    for i in range(len(ref_kmer)-2):
        if ref_kmer[i] != read[i]:
            count += 1
            if count > 3: #should i make this one?
                return -1, ''
            if ref_kmer[i] == read[i+1]:
                insertion = i+index+1
                nuc = read[i]
                first_half = read[:i]
                second_half = read[i+1:]
                read = first_half+second_half
                
    return insertion, nuc

x=16
insertions = {}
final_insertions = {}
deletions = []
for read in reads:
    read1 = read[0:x]
    read2 = read[x:2*x]
    read3 = read[2*x:3*x]
    if read1 in ref_map and read2 not in ref_map and read3 in ref_map:
        if len(ref_map[read1]) == len(ref_map[read3]):
            for i in range(len(ref_map[read1])):
                loc_on_gen = ref_map[read1][i]
                if ref_map[read3][i]-ref_map[read1][i] == 2*x-1: #insertion
                    insertion, nuc = insertion_finder(reference_genome[loc_on_gen: loc_on_gen+48], read[:48], loc_on_gen)
                    if insertion != (-1):
                        if insertion not in insertions:
                            insertions[insertion] = nuc 
                        else:
                            final_insertions[insertion] = nuc 
                #moving all the deletion code to the section below
                
        else:
            for k in range(len(ref_map[read1])):
                pos1 = ref_map[read1][k]
                for pos3 in ref_map[read3]:
                    if pos3 - pos1 < (40) and pos3 - pos1 > 0 and pos3-pos1 != 32:
                        loc = ref_map[read1][k]
                        insertion, nuc = insertion_finder(reference_genome[loc:loc+48], read[:48], loc)
                        if insertion != (-1):
                            if insertion not in insertions:
                                insertions[insertion] = nuc 
                            else:
                                final_insertions[insertion] = nuc 
                

def deletion_finder(ref_kmer, read, index):
    count = 0
    deletion = -1
    for i in range(len(ref_kmer)-2):
        if ref_kmer[i] != read[i]:
            count += 1
            if count > 2: #should i make this one?
                return -1
            if ref_kmer[i+1] == read[i]:
                deletion = i+index+1
                first_half = read[:i]
                second_half = "_" + read[i:]
                read = first_half+second_half
    return deletion

deletions = {}
final_deletions = {}
for read in reads:
    read1 = read[0:x]
    read2 = read[x:2*x]
    read3 = read[2*x:3*x]
    if read1 in ref_map and read2 in ref_map and read3 not in ref_map:
        loc = ref_map[read1][0]
        #print(ref_map[read1], ref_map[read2])
        deletion = deletion_finder(reference_genome[loc:loc+48], read[:48], loc)
        if deletion != -1:
            if deletion not in deletions:
                deletions[deletion] = reference_genome[deletion - 1] 
            else:
                final_deletions[deletion] = reference_genome[deletion - 1]


sorted_sub_dict = dict(sorted(my_dict.items()))
for key in sorted_sub_dict:
    with open('predictions.csv', 'a') as file:
        file.write(">S" + str(key) + " " + sorted_sub_dict[key] + "\n")

sorted_insertion_dict = dict(sorted(final_insertions.items()))
for key in sorted_insertion_dict:
    with open('predictions.csv', 'a') as file:
        file.write(">I" + str(key) + " " + sorted_insertion_dict[key] + "\n")

sorted_deletion_dict = dict(sorted(final_deletions.items()))
for key in sorted_deletion_dict:
    with open('predictions.csv', 'a') as file:
        file.write(">D" + str(key) + " " + sorted_deletion_dict[key] + "\n")