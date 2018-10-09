
sequences = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
             'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
             'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'];
#imports
import random
from itertools import combinations
import math as m
import numpy as np

########################## Small Functions
def min_seq_length(sequences):  #
    return min([len(x) for x in sequences])

def maximum_crossover_permutations(sequences,iteration): # Don't use, n>4 runtime impractical
    possibility_sum = 0; n = len(sequences);
    for i in range(iteration):
        possibility_sum += m.factorial(n)/(m.factorial(2)*m.factorial(n-2))*2
        n += possibility_sum/2 #### each new generation have additional nC2 combinations (due to redundancy) it could form
    return possibility_sum
'''
def generate_random_array(sequence_length): # k: length of random array
    return [random.randint(0, i) for i in sequence_length]
'''
######################### Computational Function
def random_singlept_crossover(sequences):      # cross over between 2 sequences at random point
    min_length = min_seq_length(sequences)
    random_point = random.randint(1,min_length)
    new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
    new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
    return [new_string1, new_string2]

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

def pointmutation(sequences,mutation_rate, skewrate): #mutation rate is number in a million # skew (A,C,G,T) rate
    if not skewrate: skewrate = (1,1,1,1)
    tot = sum(skewrate); dicerate = list(np.cumsum(skewrate)/tot-0.25); print(dicerate)
    bases = ['A','C','G','T']
    for seq_num in range(len(sequences)):
        for base in range(len(sequences[seq_num])):
            dice = random.random();
            if dice <= mutation_rate/1e6:
                dice2 = random.random()
                for i in range(len(dicerate)):
                    if dice2 > dicerate[i]:
                        seq_base = list(sequences[seq_num])
                        seq_base[base] = bases[i]
                        sequences[seq_num]=''.join(seq_base)
    return sequences

seqs = sequences[0:2]
print(min_seq_length(sequences))
print('maximum: '+str(maximum_crossover_permutations(sequences,3)))
print(newsequence)
print(pointmutation(newsequence,100000,[1,1,1,10]))
def sequence_motification_pipeline(sequences):
    crossover_choice = input('Would you like to perform crossover? (1 or 0) ')
    pointmutation_choice = input('Would you like to induce point mutation? (1 or 0)')
