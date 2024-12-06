reads = []
with open("project2a_spectrum.fasta" , 'r') as file:
    for line in file:
        if '>' not in line:
            line = line.replace('\n', '')
            reads.append(line)

read_map = {}
count = 0
for read in reads:
    prefix = read[:-1]
    suffix = read[1:]
    if prefix not in read_map:
        read_map[prefix] = suffix
    else:
        count += 1
        #print(prefix + ' ' + suffix) #prefix is already in read_map, could be repeated region w just small difference

prefix_to_read_number = {}
suffix_to_read_number = {}
for i in range(len(reads)):
    prefix = reads[i][:-1]
    suffix = reads[i][1:]
    prefix_to_read_number[prefix] = i
    suffix_to_read_number[suffix] = i

adjacency = {}
count_x = 0
for suffix in suffix_to_read_number:
    if suffix in prefix_to_read_number:
        suffix_number = suffix_to_read_number[suffix]
        prefix_number = prefix_to_read_number[suffix]
        adjacency[suffix_number] = prefix_number
    else:
        adjacency[suffix_number] = 4320

possible_starts = []
for key in adjacency.keys():
    if key not in adjacency.values():
        possible_starts.append(key)

for start in possible_starts:
    order = []
    test2 = True
    order.append(start)
    while test2:
        order.append(adjacency[start])
        start = adjacency[start]
        if start not in adjacency:
            test2 = False
    #print("start: " + str(order[0]) + " end: " + str(start) + " length: " + str(len(order)))

possible_starts1 = [2059, 6362, 6182, 283, 7033, 354]
for start in possible_starts1:
    order4 = []
    test4 = True
    #start = 6182
    order4.append(start)
    while test4:
        order4.append(adjacency[start])
        start = adjacency[start]
        if start not in adjacency:
            test4 = False
    for i in order4:
        with open('predictions.csv', 'a') as file:
            file.write(">read_" + str(i) + "\n")