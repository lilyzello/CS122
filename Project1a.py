from Bio import SeqIO
from Bio import AlignIO
import os

ref_genome_proj = SeqIO.read("project1a_no_repeats_reference_genome.fasta", "fasta")
ref_seq_proj = ref_genome_proj.seq

reads_proj = []
for record in SeqIO.parse("project1a_no_repeats_with_error_paired_reads.fasta", "fasta"):
    reads_proj.append(str(record.seq))

kept_reads = []
count = 0
for read in reads_proj:
    if len(read) == 50:
        kept_reads.append(read)

def equals(ref_kmer, read, loc):
    mutation = False
    location = -1
    ref_nuc = ""
    read_nuc = ""
    mutation2 = False
    location2 = -1
    ref_nuc2 = ""
    read_nuc2 = ""
    if read == ref_kmer:
        return mutation, location, ref_nuc, read_nuc, mutation2, location2, ref_nuc2, read_nuc2
    else:        
        count = 0
        for i in range(len(read)):
            if ref_kmer[i] != read[i]:
                count += 1
                if count > 2:
                    return False, -1, "", "", False, -1, "", ""
                mutation = True
                if location == -1:
                    location = i + loc
                    ref_nuc = ref_kmer[i]
                    read_nuc = read[i]
                else:
                    mutation2 = True
                    location2 = i + loc
                    ref_nuc2 = ref_kmer[i]
                    read_nuc2 = read[i]
            if i == len(read) - 1:
                return mutation, location, ref_nuc, read_nuc, mutation2, location2, ref_nuc2, read_nuc2

def is_match2(ref_gen, reads):
    all_mutations = dict()
    final_mut = dict()
    for read in reads:
        #if len(read) < 16:
            #break
        for i in range(len(ref_gen)-len(read)):
            mutation, location, ref_nuc, read_nuc, mutation2, location2, ref_nuc2, read_nuc2= equals(ref_gen[i:i+len(read)], read, i)
            if mutation:
                if location in all_mutations.keys():
                    final_mut[location] = ref_nuc + " " + read_nuc
                else:
                    all_mutations[location] = ref_nuc + " " + read_nuc
            if mutation2:
                if location2 in all_mutations.keys():
                    final_mut[location2] = ref_nuc2 + " " + read_nuc2
                else:
                    all_mutations[location2] = ref_nuc2 + " " + read_nuc2
    return final_mut


mut_dict_proj2 = is_match2(ref_seq_proj, kept_reads)

sorted_dict_by_keys_proj2 = dict(sorted(mut_dict_proj2.items()))
#print(len(sorted_dict_by_keys_proj2))

final = []
values = []
for key in sorted_dict_by_keys_proj2:
    final.append(str(key))
for value in sorted_dict_by_keys_proj2.values():
    values.append(value)

for i in range(len(sorted_dict_by_keys_proj2)):
    with open('predictions.csv', 'a') as file:
        file.write(">S" + final[i] + " " + values[i] + "\n")
