import numpy as np
def RandomMotifs(Dna, k, t):
    import random
    return [text[random.randint(0,len(text)-k):][0:k] for text in Dna]

def Consensus(Motifs):
    count,k = CountWithPseudocounts(Motifs)
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

def Score(Motifs):
    count,k = CountWithPseudocounts(Motifs)
    median_motif = Consensus(Motifs)
    score = 0;
    #print(count)
    for i in range(k):
        for j in 'ACGT':
            if j != median_motif[i]:
                    score += count[j][i]
    return score

def Pr(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob *= Profile[Text[i]][i]
    return prob

def CountWithPseudocounts(Motifs):   # Count Parallel Motifs Content (get Median)
    count = {} # initializing the count dictionary
    k = len(Motifs[0]); t = len(Motifs);
    for symbol in "ACGT":
        count[symbol] = []
        count[symbol] = list(np.ones(k,dtype = 'int'))
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count, k

def ProfileWithPseudocounts(Motifs):
    count, k = CountWithPseudocounts(Motifs)
    num_seq = len(Motifs)+4
    for i in range(k):
        for j in 'ACGT':
            count[j][i] = float(count[j][i]/num_seq)
    return count

def ProfileMostProbableKmer(Text, k, Profile):
    probs = []; topprobs_seq = [];
    for i in range(len(Text)-k+1):
        probs.append(Pr(Text[i:i+k],Profile))
    #print(probs)
    for j in range(len(probs)):
        if probs[j] == max(probs):
            topprobs_seq=Text[j:j+k]
            #print(max(probs))
    #topprob_index = probs.index(topprob)
        if max(probs) == 0:
            topprobs_seq = Text[0:k]
    return topprobs_seq

def Motifs(Profile, Dna,k):
    motifs = []
    for i in range(len(Dna)):
        mot = ProfileMostProbableKmer(Dna[i],k,Profile);
        motifs.append(mot)
    return motifs

#### N is the number of reiteration
def RandomizedMotifSearch(Dna, k, N):
    t = len(Dna);
    M = RandomMotifs(Dna, k, t) # ['CCA','CCT','CTT','TTG'] # Customizable Initial kmer
    BestMotifs = M
    for i in range(N):
        Motifprofile = ProfileWithPseudocounts(M)
        #print(Profile)
        M = Motifs(Motifprofile, Dna,k)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs, Score(BestMotifs)

################ Gibb's
def GreedyMotifSearchWithPseudocounts(Dna, k):
    t = len(Dna); BestMotifs = [];
    for i in range(0, len(Dna)):
        BestMotifs.append(Dna[i][0:k])

    for i in range(len(Dna[0])-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        #print(Motifs)
        for j in range(1, len(Dna)):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
            #print(Motifs)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
        bestscore = Score(BestMotifs)
    return BestMotifs, bestscore

Dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
       'TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
       'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']