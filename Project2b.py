# In[2]:


reads = []
with open("project2b_reads.fasta" , 'r') as file:
    for line in file:
        if '>' not in line:
            line = line.replace('\n', '')
            reads.append(line)


# In[3]:


print(len(reads))


# In[4]:


count_map = {}
read_nums = {}
reads_to_num = {}
for i in range(len(reads)):
    read = reads[i]
    if read[:30] not in reads_to_num:
        reads_to_num[read[:30]] = i
    #else: probably a repeated region
    for k in range(len(read)-20):
        kmer = read[k:k+20]
        if kmer not in count_map:
            count_map[kmer] = 1
        else:
            count_map[kmer] += 1
        if kmer not in read_nums:
            read_nums[kmer] = [i]
        else:
            read_nums[kmer].append(i)
print("first: ", len(count_map))

removed = []
spectrum = []
for key in count_map:
    if count_map[key] < 6:
        if read_nums[key][0] not in removed:
            removed.append(read_nums[key][0])
    else:
        spectrum.append(key)
print("second: ", len(spectrum))


# In[21]:


removed = []
spectrum = []
for key in count_map:
    if count_map[key] < 3:
        if read_nums[key][0] not in removed:
            removed.append(read_nums[key][0])
    else:
        spectrum.append(key)
print("second: ", len(spectrum))


# In[22]:


read_nums_to_read = {}
for i in range(len(reads)):
    read_nums_to_read[i] = reads[i]
for removed_read in removed:
    read_nums_to_read.pop(removed_read)


# In[34]:


prefix_to_suffix_map = {}
count = 0
for kmer in spectrum:
    prefix = kmer[:-1]
    suffix = kmer[1:]
    if prefix not in prefix_to_suffix_map:
        prefix_to_suffix_map[prefix] = suffix
    else:
        count += 1
print(len(prefix_to_suffix_map))


# In[33]:


adjacency = {}
missing_count = 0
for prefix in prefix_to_suffix_map:
    suffix = prefix_to_suffix_map[prefix] 
    if suffix in prefix_to_suffix_map: 
        adjacency[suffix] = prefix_to_suffix_map[suffix] 
    else:
        missing_count += 1
        adjacency[suffix] = "end"
print(missing_count)
print(len(adjacency))


# In[25]:


possible_starts = []
for key in adjacency.keys():
    if key not in adjacency.values():
        possible_starts.append(key)
print(possible_starts)


# In[26]:


start_dict = {}
for start in possible_starts:
    order = []
    test2 = True
    order.append(start)
    while test2:
        order.append(adjacency[start])
        start = adjacency[start]
        if start not in adjacency:
            test2 = False
    start_dict[order[0]] = len(order)
    #print("start: " + str(order[0]) + " end: " + str(start) + " length: " + str(len(order)))
sorted_starts = dict(sorted(start_dict.items(), key=lambda x:x[1], reverse=True)) #sort starts so that longest sequence comes first
print(sorted_starts)


# In[27]:


total = 0
for item in sorted_starts.values():
    total += item
print(total)


# In[40]:


possible_starts = list(sorted_starts.keys())
genome = ""
for start in possible_starts:
    order4 = []
    test4 = True
    genome += start
    while test4:
        #order4.append(adjacency[start])
        if adjacency[start] != 'end':
            genome += adjacency[start][-1]
        start = adjacency[start]
        if start not in adjacency:
            test4 = False
            #genome += '_'
print(len(genome))
#next step is matching reads to reference genome


# In[36]:


half_reads = {}
for i in range(len(reads)):
    read = reads[i]
    first_half = read[:20]
    second_half = read[20:40]
    if first_half not in half_reads:
        half_reads[first_half] = [i]
    else:
        half_reads[first_half].append(i)
    if second_half not in half_reads:
        half_reads[second_half] = [i]
    else:
        half_reads[second_half].append(i)


# In[37]:


final = []
for i in range(len(genome)):
    if genome[i:i+30] in reads_to_num:
        final.append(reads_to_num[genome[i:i+30]])
    elif genome[i:i+20] in half_reads:
        for j in range(len(half_reads[genome[i:i+20]])): 
            final.append(half_reads[genome[i:i+20]][j])


# In[54]:


reads = []
read_nums = []
with open("project2b_reads.fasta", 'r') as file:
    for line in file:
        if '>' in line:
            line = line.replace('\n', '')
            read_nums.append(line)
        else:
            line = line.replace('\n', '')
            reads.append(line)

ref_gen = genome

ref_map = {}
for i in range(len(ref_gen)-16):
    to_add = ref_gen[i:i+16]
    if to_add not in ref_map:
         ref_map[to_add] = [i]

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
        read_map[i] = ref_map[read1]
    elif case1:
        count1 += 1
        read_map[i] = ref_map[read1]
    elif case2:
        count2 += 1
        read_map[i] = ref_map[read1]
    elif case3:
        count3 += 1
        read_map[i] = ref_map[read2]
    elif read1 in ref_map:
        count4 += 1
        read_map[i] = ref_map[read1]
    elif read2 in ref_map:
        count5 += 1
        read_map[i] = ref_map[read2]
    elif read3 in ref_map:
        count6 += 1
        read_map[i] = ref_map[read3]
    else: #no matches
        count7 += 1
        reads_nomatch.append(read)
        reads_nomatch_nums.append(i)
    

ref_map2 = {}
for i in range(len(ref_gen)-10):
    to_add = ref_gen[i:i+10]
    if to_add not in ref_map2:
        ref_map2[to_add] = [i]
    else:
        ref_map2[to_add].append(i)

counter_none = 0
for j in range(len(reads_nomatch)):
    read = reads_nomatch[j]
    r1 = read[0:10]
    r2 = read[10:20]
    r3 = read[20:30]
    r4 = read[30:40]
    read_num = reads_nomatch_nums[j]
    if r1 in ref_map2 and r2 in ref_map2 and r3 in ref_map2 and r4 in ref_map2:
        read_map[read_num] = ref_map2[r1]
    elif r1 in ref_map2 and r2 in ref_map2 and r3 in ref_map2:
        read_map[read_num] = ref_map2[r1]
    elif r1 in ref_map2 and r2 in ref_map2 and r4 in ref_map2:
        read_map[read_num] = ref_map2[r1]
    elif r1 in ref_map2 and r3 in ref_map2 and r4 in ref_map2:
        read_map[read_num] = ref_map2[r1]
    elif r2 in ref_map2 and r3 in ref_map2 and r4 in ref_map2:
        read_map[read_num] = ref_map2[r2]
    elif r1 in ref_map2:
        if r2 in ref_map2:
            read_map[read_num] = ref_map2[r1]
        if r3 in ref_map2:
            read_map[read_num] = ref_map2[r1]
        if r4 in ref_map2:
            read_map[read_num] = ref_map2[r1]
        else:
            read_map[read_num] = ref_map2[r1]
    elif r2 in ref_map2:
        read_map[read_num] = ref_map2[r2]
    elif r3 in ref_map2:
        if r4 in ref_map2:
            read_map[read_num] = ref_map2[r3]
        else:
            read_map[read_num] = ref_map2[r3]
    elif r4 in ref_map2:
        read_map[read_num] = ref_map2[r4]
    else: #no matches
        counter_none += 1

#print(counter_none)
#print(read_map)
#print(read_map)
readmap2 = {}
for i in read_map:
    readmap2[i] = read_map[i][0]
sorted_map = dict(sorted(readmap2.items(), key=lambda x:x[1]))
#print(sorted_map)
final2 = []
for i in sorted_map:
    final2.append(i)
for i in final2:
    with open('predictions.csv', 'a') as file:
        file.write(">read_" + str(i) + "\n")

