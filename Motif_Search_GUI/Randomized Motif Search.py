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

def Entropy(Motifs):
    ()
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

import math as m
def ProfileWithPseudocounts(Motifs):
    count, k = CountWithPseudocounts(Motifs)
    num_seq = len(Motifs)+4; entropy = 0;
    for i in range(k):
        for j in 'ACGT':
            count[j][i] = float((float(count[j][i])/float(num_seq)))
        f = max(count[j][i] for j in 'ACGT')
        entropy += -f*m.log(f)
    return count, entropy


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
def RandomizedMotifSearch(Dna, k_iter, N):
    k_bestmotifs = []; k_scores = []; k_entropy = [];
    for k in range(2, k_iter):
        t = len(Dna); M = RandomMotifs(Dna, k, t); BestMotifs = M
        for i in range(N):
            Motifprofile, entropy = ProfileWithPseudocounts(M)
            #print(Profile)
            M = Motifs(Motifprofile, Dna,k);
            this_score = Score(M)
            if this_score < Score(BestMotifs):
                BestMotifs = M
                BestP, Bestentropy = Motifprofile, entropy
        k_bestmotifs.append(BestMotifs)
        k_scores.append(this_score / k)
        k_entropy.append(Bestentropy / k)
    return k_bestmotifs, k_scores, k_entropy

################ Gibb's
def GreedyMotifSearchWithPseudocounts(Dna, k_iter):
    k_bestmotifs = []; k_scores = []; k_entropy = [];
    for k in range(2, k_iter):
        t = len(Dna); BestMotifs = []; loop_scores = [];
        for i in range(0, len(Dna)):
            BestMotifs.append(Dna[i][0:k])

        for i in range(len(Dna[0]) - k + 1):
            Motifs = []
            Motifs.append(Dna[0][i:i + k])
            # print(Motifs)
            for j in range(1, len(Dna)):
                P, entropy = ProfileWithPseudocounts(Motifs[0:j])
                Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
            # print(Motifs)
            this_score = Score(Motifs)
            if this_score < Score(BestMotifs):
                BestMotifs = Motifs
                loop_scores.append(this_score)
                BestP, entropy = P, entropy
        k_bestmotifs.append(BestMotifs)
        k_scores.append(loop_scores[-1]/k)
        k_entropy.append(entropy/k)
    return k_bestmotifs, k_scores, k_entropy

Dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
       'TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
       'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']


import ast
def entry_processor(Dna_entry, k_entry, N_entry): ### take .entry as input
    Dna_disp = str(Dna_entry.get());
    seq = ast.literal_eval(Dna_disp);  # seq = [n.strip() for n in seq_list];
    k = int(k_entry.get());
    if not N_entry.get(): N = 100;
    else: N = int(N_entry.get());
    return [seq, k, N]

############################################### GUI
def display_motif_Random():
    [seq, k, N] = entry_processor(Dna_entry, k_entry,N_entry);
    k_bestmotifs, k_scores, k_entropy = RandomizedMotifSearch(seq,k,N)
    T.insert(END, 'Motifs:\n '+str(k_bestmotifs)+ '\n\nScores:\n '+str(k_scores)+
             '\n\nEntropy:\n '+str(k_entropy)+'\n')
    Status_done()

def display_motif_Gibbs():
    [seq, k, N] = entry_processor(Dna_entry, k_entry,N_entry);
    k_bestmotifs, k_scores, k_entropy = GibbsSampler(seq, k, N)
    T.insert(END, 'Motifs:\n '+str(k_bestmotifs)+ '\n\nScores:\n '+str(k_scores)+
             '\n\nEntropy:\n '+str(k_entropy)+'\n')
    Status_done()

def display_motif_Greedy():
    seq, k,_ = entry_processor(Dna_entry, k_entry,N_entry);
    k_bestmotifs, k_scores, k_entropy = GreedyMotifSearchWithPseudocounts(seq, k)
    T.insert(END, 'Motifs:\n '+str(k_bestmotifs)+ '\n\nScores:\n '+str(k_scores)+
             '\n\nEntropy:\n '+str(k_entropy)+'\n')
    Status_done()

###################################################### GUI
##################### Entry
from tkinter import *
root = Tk()
root.geometry("800x700"); root.title('Holy Moti');
root.iconbitmap(r"C:\Users\mailb\PycharmProjects\Motif Search Function\apta_index_tn_13c_icon.ico")
Label(root, text= "Sequence: (Format: List)").grid(row=0); Dna_entry = Entry(root);
Dna_entry.bind("<Return>");Dna_entry.grid(row = 0, column = 1);

Label(root, text= "Target Motif Length").grid(row=1);k_entry = Entry(root);
k_entry.bind("<Return>"); k_entry.grid(row = 1, column = 1);

Label(root, text= "No. of Reiteration (optional)").grid(row=2); N_entry = Entry(root);
N_entry.bind("<Return>"); N_entry.grid(row = 2, column = 1);

##################### Buttons
frame1 = Frame(root); frame1.grid(row = 4, column = 0);
compute_button1 = Button(frame1, text = 'Randomized\nMotif Search', width = 15, height = 2, fg = 'blue',
                        command = lambda:Status_running(1));
compute_button1.pack(side = LEFT);

frame2 = Frame(root); frame2.grid(row = 4, column = 1);
compute_button2 = Button(frame2, text = 'Gibb''s Sampling\nMotif Search', width = 15, height = 2, fg = 'green',
                        command = lambda:Status_running(2));
compute_button2.pack()

frame3 = Frame(root); frame3.grid(row = 4, column = 2);
compute_button3 = Button(frame3, text = 'Greedy\nMotif Search', width = 15, height = 2, fg = 'purple',
                        command = lambda:Status_running(3));
compute_button3.pack()

def reset():
    T.delete('1.0', END); Status_done();

reset_button = Button(root, text = 'Reset Space', width = 15, height = 2, fg = 'black',
                        command = reset).grid(row = 5, column = 2);
exit_button = Button(root, text = 'Exit', width = 15, height = 2, fg = 'red',
                        command = exit).grid(row = 5, column = 1);

label_status = Label(text = 'Status: ');label_status.grid(row = 6, column = 1);

T = Text(root, height=30, width=99)
T.grid(row = 7, column = 0, columnspan=3, rowspan=2, padx=2, pady=3)

def Status_running(number):
    label_status.config(text = 'Searching for the bottom of your patience...')
    if number == 1 : display_motif_Random()
    if number == 2 : display_motif_Gibbs()
    if number == 3: display_motif_Greedy()

def Status_done():
   label_status.config(text = 'Status: Done! I''m ready.')

root.mainloop()