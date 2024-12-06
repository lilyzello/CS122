reads = []
read_nums = []
with open("project1c_reads.fasta", 'r') as file:
    for line in file:
        if '>' in line:
            line = line.replace('\n', '')
            read_nums.append(line)
        else:
            line = line.replace('\n', '')
            reads.append(line)

file_names = []
acceptable = [0, 12, 31, 63, 73, 84, 85, 67, 91]
for i in acceptable:
    file_names.append("project1c_genome_" + str(i) + ".fasta")

ref_gens = []
for gen_file in file_names:
    reference_genome = ""
    with open(gen_file, 'r') as file:
        for line in file:
          if '>' not in line:
              line = line.replace('\n', '')
              reference_genome += line
    ref_gens.append(reference_genome)

ref_map = {}
for ref_gen_num in range(len(ref_gens)):
    ref_gen = ref_gens[ref_gen_num]
    number = acceptable[ref_gen_num] #new
    for i in range(len(ref_gen)-16):
        to_add = ref_gen[i:i+16]
        if to_add not in ref_map:
            #ref_map[to_add] = [ref_gen_num]
            ref_map[to_add] = number
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
    

ref_map2 = {}
for ref_gen_num in range(len(ref_gens)):
    ref_gen = ref_gens[ref_gen_num]
    for i in range(len(ref_gen)-10):
        to_add = ref_gen[i:i+10]
        if to_add not in ref_map2:
            ref_map2[to_add] = [ref_gen_num]
        else:
            ref_map2[to_add].append(ref_gen_num)

counter_none = 0
for j in range(len(reads_nomatch)):
    read = reads_nomatch[j]
    r1 = read[0:10]
    r2 = read[10:20]
    r3 = read[20:30]
    r4 = read[30:40]
    read_num = reads_nomatch_nums[j]
    if r1 in ref_map2 and r2 in ref_map2 and r3 in ref_map2 and r4 in ref_map2:
        read_map[read_num] = [ref_map2[r1], ref_map2[r2], ref_map2[r3], ref_map2[r4]]
    elif r1 in ref_map2 and r2 in ref_map2 and r3 in ref_map2:
        read_map[read_num] = [ref_map2[r1], ref_map2[r2], ref_map2[r3]]
    elif r1 in ref_map2 and r2 in ref_map2 and r4 in ref_map2:
        read_map[read_num] = [ref_map2[r1], ref_map2[r2], ref_map2[r4]]
    elif r1 in ref_map2 and r3 in ref_map2 and r4 in ref_map2:
        read_map[read_num] = [ref_map2[r1], ref_map2[r3], ref_map2[r4]]
    elif r2 in ref_map2 and r3 in ref_map2 and r4 in ref_map2:
        read_map[read_num] = [ref_map2[r2], ref_map2[r3], ref_map2[r4]]
    elif r1 in ref_map2:
        if r2 in ref_map2:
            read_map[read_num] = [ref_map2[r1], ref_map2[r2]]
        if r3 in ref_map2:
            read_map[read_num] = [ref_map2[r1], ref_map2[r3]]
        if r4 in ref_map2:
            read_map[read_num] = [ref_map2[r1], ref_map2[r4]]
        else:
            read_map[read_num] = [ref_map2[r1]]
    elif r2 in ref_map2:
        if r3 in ref_map2:
            read_map[read_num] = [ref_map2[r2], ref_map2[r3]]
        if r4 in ref_map2:
            read_map[read_num] = [ref_map2[r2], ref_map2[r4]]
        else:
            read_map[read_num] = [ref_map2[r2]]
    elif r3 in ref_map2:
        if r4 in ref_map2:
            read_map[read_num] = [ref_map2[r3], ref_map2[r4]]
        else:
            read_map[read_num] = [ref_map2[r3]]
    elif r4 in ref_map2:
        read_map[read_num] = [ref_map2[r4]]
    else: #no matches
        counter_none += 1

#print(counter_none)
        

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
    with open('predictions2.csv', 'a') as file:
        file.write(">read_" + str(key) + " " + "Genome_Number_" + str(sorted_dict[key]) + "\n")