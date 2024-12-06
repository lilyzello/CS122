reads = []
read_nums = []
with open("project1d_reads.fasta", 'r') as file:
    for line in file:
        if '>' in line:
            line = line.replace('\n', '')
            read_nums.append(line)
        else:
            line = line.replace('\n', '')
            reads.append(line)

file_names = []
acceptable = [709, 798, 816, 933, 944, 1032, 1031, 1108, 1148, 1185, 1233, 1480, 1440, 1509, 1873, 1941, 2199, 2164, 2182, 2138, 2203, 2384, 2300, 2416, 2563, 2670, 2488, 2920, 3141, 3270, 3423, 3458, 3631, 3640, 3683, 4037, 4363, 4457, 4614, 4652, 4674, 4871]
for i in acceptable:
    file_names.append("project1d_genome_" + str(i) + ".fasta")

ref_gens = []
for gen_file in file_names:
    reference_genome = ""
    with open(gen_file, 'r') as file:
        for line in file:
          if '>' not in line:
              line = line.replace('\n', '')
              reference_genome += line
    ref_gens.append(reference_genome)

#ref_map = {}
#for j in range(3600,3900):
#    ref_gen = ref_gens[j]
#    for i in range(len(ref_gen)-16):
#        to_add = ref_gen[i:i+16]
#        if to_add not in ref_map:
#            ref_map[to_add] = [j]
#        else:
#            ref_map[to_add].append(j)
ref_map = {}
for ref_gen_num in range(len(ref_gens)):
    ref_gen = ref_gens[ref_gen_num]
    number = acceptable[ref_gen_num] #new
    for i in range(len(ref_gen)-16):
        to_add = ref_gen[i:i+16]
        if to_add not in ref_map:
            #ref_map[to_add] = [ref_gen_num]
            ref_map[to_add] = [number]
        else:
            #ref_map[to_add].append(ref_gen_num)
            ref_map[to_add].append(number)

reads_nomatch = []
reads_nomatch_nums = []
count0 = 0
count1 = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0
count6 = 0
count7 = 0
read_map = {}
#for read in reads:
for i in range(len(reads)):
    read = reads[i]
    read1 = read[0:16]
    read2 = read[16:32]
    read3 = read[32:48]
    case1 = read1 in ref_map and read2 in ref_map and read3 not in ref_map #[1,2] match
    case2 = read1 in ref_map and read2 not in ref_map and read3 in ref_map #[1,3] match
    case3 = read1 not in ref_map and read2 in ref_map and read3 in ref_map #[2,3] match
    
    if read1 in ref_map and read2 in ref_map and read3 in ref_map: #all hash
        count0 += 1
        read_map[i] = [ref_map[read1], ref_map[read2], ref_map[read3]]
    elif case1:
        count1 += 1
        read_map[i] = [ref_map[read1], ref_map[read2]]
    elif case2:
        count2 += 1
        read_map[i] = [ref_map[read1], ref_map[read3]]
    elif case3:
        count3 += 1
        read_map[i] = [ref_map[read2], ref_map[read3]]
    elif read1 in ref_map:
        count4 += 1
        read_map[i] = [ref_map[read1]]
    elif read2 in ref_map:
        count5 += 1
        read_map[i] = [ref_map[read2]]
    elif read3 in ref_map:
        count6 += 1
        read_map[i] = [ref_map[read3]]
    else: #no matches
        count7 += 1
        read_map[i] = ['N']
        reads_nomatch.append(read)
        reads_nomatch_nums.append(i)

from collections import Counter
for key in read_map:
    value = read_map[key]
    arr = []
    if len(value) == 1:
        for k in value:
            arr.append(k[0])
    if len(value) > 1:
        for x in value:
            if len(x) == 1:
                arr.append(x[0])
            else:
                for z in x:
                    arr.append(z)
    count = Counter(arr)
    value = count.most_common(1)[0][0]
    read_map[key] = value

sorted_dict = dict(sorted(read_map.items()))
for key in sorted_dict:
    with open('predictions.csv', 'a') as file:
        file.write(">read_" + str(key) + " " + "Genome_Number_" + str(sorted_dict[key]) + "\n")

test_dict = {}
for key in sorted_dict:
    if sorted_dict[key] not in test_dict:
        test_dict[sorted_dict[key]] = 0
    else:
        test_dict[sorted_dict[key]] += 1
print(dict(sorted(test_dict.items(), key=lambda x:x[1], reverse=True)))