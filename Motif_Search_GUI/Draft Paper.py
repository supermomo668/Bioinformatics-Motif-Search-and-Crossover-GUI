from itertools import combinations
import random
def random_singlept_crossover(sequences):      # cross over between 2 sequences at random point
    min_length = min_seq_length(sequences)
    random_point = random.randint(1,min_length)
    new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
    new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
    return [new_string1, new_string2]
def min_seq_length(sequences):  #
    return min([len(x) for x in sequences])

def crossover(sequences, cross_percent, iteration, max_cap):
    # cross_percent: number of population crossover per iteration #max_cap: maximum number of sequence allowed to generate
    if not max_cap: max_cap = 1000;
    new_sequences = []; generation_pool = sequences; #initialize the first pool to cross over
    for i in range(iteration):
        pairpool = list(combinations(generation_pool,2))
        pairpool_selected = [x for x in pairpool if random.random()< cross_percent] # choose a portion from the combination pool
        generation_pool=[]
        for pair in pairpool_selected:
            if len(generation_pool) >= max_cap: break ### condition to stop adding sequences
            newpair = random_singlept_crossover(pair)
            generation_pool.append(newpair)

        generation_pool= [seq for seqpair in generation_pool for seq in seqpair]
        new_sequences.append(generation_pool)
    return [seq for sublist in new_sequences for seq in sublist]
sequences = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
             'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
             'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'];

print(crossover(sequences, 0.4,3,1000))

#########################################################
import numpy as np
def Consensus(Motifs):
    count,k = CountWithPseudocounts(Motifs,1)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def CountWithPseudocounts(Motifs,psuedo):   # Count Parallel Motifs Content (get Median)
    count = {} # initializing the count dictionary
    k = len(Motifs[0]); t = len(Motifs);
    for symbol in "ACGT":
        count[symbol] = []
        if psuedo == 1:
            count[symbol] = list(np.ones(k,dtype = 'int'))
        elif psuedo == 0:
            count[symbol] = list(np.zeros(k,dtype = 'int'))
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count, k

def plainScore(Motifs):
    count, k = CountWithPseudocounts(Motifs, 0)
    median_motif = Consensus(Motifs)
    score = 0;
    # print(count)
    for i in range(k):
        for j in 'ACGT':
            if j != median_motif[i]:
                score += count[j][i]
    return score

def findmotif(single_seq,motif):
    lowest_scores = []
    for mot in motif:
        mot_scores = []; k = len(mot);
        for window in range(len(single_seq) - k + 1):
            mot_scores.append(plainScore([single_seq[window:window + k], mot]))
        lowest_scores.append(mot_scores.index(min(mot_scores)))
    return lowest_scores  # [3,2,4,..]

print(findmotif(sequences[0],['ACCGT','AGGTAC']))

