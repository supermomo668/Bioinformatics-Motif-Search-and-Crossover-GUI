######################################
'''
This program is created by Matthew Mo.
The program consist of a python GUI that could be use to select highly conserved or probable motifs among sequence (or
aptamer) of similar lengths. 3 alogrithmic approaches could be used to find the best motifs best represent all the sequence
given an unknown length for the k-mer. Graphical representation of motif's quality will facilitate user to narrow down
and determine best consensus sequence.
'''
import numpy as np
import random
import matplotlib.pyplot as plt
import time
from tkinter import font
from itertools import combinations
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

def Score(Motifs):
    count,k = CountWithPseudocounts(Motifs,1)
    median_motif = Consensus(Motifs)
    score = 0;
    #print(count)
    for i in range(k):
        for j in 'ACGT':
            if j != median_motif[i]:
                    score += count[j][i]
    return score
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
def Pr(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob *= Profile[Text[i]][i]
    return prob

import math as m
def ProfileWithPseudocounts(Motifs):
    count, k = CountWithPseudocounts(Motifs,1);
    num_seq = len(Motifs)+4; totentropy = 0; probcount = count;
    for i in range(k):
        for j in 'ACGT':
            probcount[j][i] = count[j][i] / num_seq
            f = probcount[j][i]
            if f != 0:
                totentropy += -f * m.log2(f)
    return probcount, totentropy
def Entropy4Logo(Motifs):
    count, k = CountWithPseudocounts(Motifs,0)
    num_seq = len(Motifs) + 4;
    entropy_logo = []
    for i in range(k):
        entropy_logo.append([])
        for j in 'ACGT':
            f = float(count[j][i]/num_seq)
            if f > 0:
                entropy_logo[i].append((j,-f*m.log(f)))
            else: entropy_logo[i].append((j,0))
    return entropy_logo
################################################################## Randomized
def Motifs(Profile, Dna,k):
    motifs = []
    for i in range(len(Dna)):
        mot = ProfileMostProbableKmer(Dna[i],k,Profile);
        motifs.append(mot)
    return motifs

def RandomMotifs(Dna, k, t):
    return [text[random.randint(0, len(text) - k):][0:k] for text in Dna]

################################################################## Gibb's
def Normalize(Probabilities):
    return {key:Probabilities[key]/sum(Probabilities.values()) for key in Probabilities}
def WeightedDie(Probabilities):
    num = random.uniform(0,1);
    cumprob = 0;
    for i in Probabilities:
        cumprob += Probabilities[i]
        if num < cumprob:
            return i
def Pr(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob *= Profile[Text[i]][i]
    return prob
def ProfileGeneratedString(Text, Profile, k):
    probabilities = {};
    for i in range(0,len(Text)-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], Profile)
    return WeightedDie(Normalize(probabilities))

################################################################### Greedy
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

################################################################## Search Main Functions
#### N is the number of reiteration
############# Randomized
def RandomizedMotifSearch(Dna, k_iter, N):
    k_bestmotifs = []; k_scores = []; k_entropy = [];
    for k in range(2, k_iter+1):
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
def GibbsSampler(Dna, k_iter, N):
    BestMotifs = []; t = len(Dna);
    k_bestmotifs = []; k_scores = []; k_entropy = [];
    for k in range(2,k_iter+1):
        RandMotifs = RandomMotifs(Dna, k, t);
        BestMotifs = RandMotifs;
        for i in range(1, N):
            tsided_dice = random.randint(0, t - 1);
            removed_seq = Dna[tsided_dice];
            remained_seqs = Dna[0:tsided_dice] + Dna[tsided_dice + 1:]
            motif_removed = RandMotifs[tsided_dice]
            motif_selected = RandMotifs[0:tsided_dice] + RandMotifs[tsided_dice + 1:]
            motifprofile, entropy = ProfileWithPseudocounts(motif_selected);
            generated_motif = ProfileGeneratedString(removed_seq, motifprofile, k)
            motif_selected.append(generated_motif);
            Motifs = motif_selected
            this_score = Score(Motifs)
            if this_score < Score(BestMotifs):
                BestMotifs = Motifs
                Bestentropy = entropy
        k_bestmotifs.append(BestMotifs);
        k_scores.append(this_score / k);
        k_entropy.append(Bestentropy / k);

    return k_bestmotifs, k_scores, k_entropy

################ Greedy
def GreedyMotifSearchWithPseudocounts(Dna, k_iter):
    k_bestmotifs = []; k_scores = []; k_entropy = [];
    for k in range(2, k_iter+1):
        t = len(Dna); BestMotifs = []; loop_scores = [];
        for i in range(0, len(Dna)):
            BestMotifs.append(Dna[i][0:k])
        for i in range(len(Dna[0]) - k + 1):
            Motifs = []
            Motifs.append(Dna[0][i:i + k])

            for j in range(1, len(Dna)):
                P, entropy = ProfileWithPseudocounts(Motifs[0:j])
                Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
                this_score = Score(Motifs)
            if this_score < Score(BestMotifs):
                BestMotifs = Motifs
                loop_scores.append(this_score)
                BestP, entropy = P, entropy
        k_bestmotifs.append(BestMotifs)
        k_scores.append(loop_scores[-1]/k)
        k_entropy.append(entropy/k)
    return k_bestmotifs, k_scores, k_entropy

Dna ='CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',\
     'TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',\
     'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
######################################################################################################### GUI
############################################################## GUI functions
import ast
import tkinter.messagebox
def upload_list():
    return seq_list.get()

def entry_processor(Dna_entry, k_entry, N_entry): ### take .entry as input
    if Dna_entry.get():
        Dna_disp = Dna_entry.get(); seq = list(ast.literal_eval(Dna_disp)); # turn a 'list string' into real list
    else:
        seq = None;
    if k_entry.get(): k = int(k_entry.get());
    else: k = None
    if not N_entry.get():
        N = 100;
    else:
        N = int(N_entry.get());
    return [seq, k, N]

def display_results(seq, k, delta_t, k_bestmotifs,k_scores, k_entropy, searchname):
    T.insert(END, '\nYou ran: '+str(searchname) +' (k = ' + str(len(seq[0])) + ')\n| k = ' + str(k) +
             '| Number of Sequence: ' + str(len(seq)) + '| Runtime: ' + delta_t + 's\n\n')
    T.insert(END, 'Motifs:\n ' + str(k_bestmotifs) + '\n\nScores:\n' + str(k_scores) + '\n\nEntropy:\n ' + str(k_entropy) + '\n')
    scores.set(k_scores); entropy.set(k_entropy);
    return scores, entropy, k_bestmotifs

def upload_choicegate(upload):
    t0 = time.time(); #print('Upload is: ' +str(upload));
    if not k_entry.get() or (not Dna_entry.get() and upload ==0) or (upload == 0 and int(k_entry.get())>len(ast.literal_eval(Dna_entry.get())[0])):
        tkinter.messagebox.showinfo('Opps','You forgot something...');
        return None
    else:
        if upload == 1:
            _, k, N = entry_processor(None, k_entry, N_entry)
            seq = seq_list.get()
        else:
            seq, k, N = entry_processor(Dna_entry, k_entry, N_entry);
        return seq, k, N,t0

def display_motif_Random(upload):
    if upload_choicegate(upload) == None: return None
    seq, k, N, t0 = upload_choicegate(upload); searchname = 'Randomized Motif Search'
    k_bestmotifs, k_scores, k_entropy = RandomizedMotifSearch(seq, k, N)
    delta_t = str(time.time() - t0);
    return T.delete('1.0',END), display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy,searchname);  ####

def display_motif_Gibbs(upload):
    if upload_choicegate(upload) == None: return None
    seq, k, N, t0 = upload_choicegate(upload); searchname = 'Gibb\'s Randomized Sampler'
    k_bestmotifs, k_scores, k_entropy = GibbsSampler(seq, k, N)
    delta_t = str(time.time() - t0);
    return T.delete('1.0',END), display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy,searchname);  ####

def display_motif_Greedy(upload):
    if upload_choicegate(upload) == None: return None
    seq,k,N,t0 = upload_choicegate(upload); searchname = 'Greedy Motif Search'
    k_bestmotifs, k_scores, k_entropy = GreedyMotifSearchWithPseudocounts(seq, k)
    delta_t = str(time.time() - t0);
    return T.delete('1.0',END), display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy,searchname);##

###################################### Crossover and mutation
def min_seq_length(sequences):  #
    return min([len(x) for x in sequences])

######################### Computational Function
def random_singlept_crossover(sequences):      # cross over between 2 sequences at random point
    min_length = min_seq_length(sequences)
    random_point = random.randint(1,min_length)
    new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
    new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
    return [new_string1, new_string2]

def defined_singlept_crossover(sequences,position_dict,lengths):   # cross over between 2 sequences at random point
    min_length = min_seq_length(sequences);
    positions1 = position_dict[sequences[0]]; positions2 = position_dict[sequence[1]]
    for pos1, pos2,length in zip(positions1, positions2,lengths):
        random_point = random.randint(1, min_length)
        if ((random_point > pos1 and (random_point - pos1) < length) and
            (random_point > pos2 and (random_point - pos2) < length)):
            random_point = random_point + length;
            if random_point > min_length:
                return sequences
            new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
            new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
            return [new_string1, new_string2]
        else:
            return sequences

def crossover(sequences, cross_percent, iteration, max_cap,defined):
    # cross_percent: number of population crossover per iteration #max_cap: maximum number of sequence allowed to generate
    if defined == 1:
        motifs = list(ast.literal_eval(motif_entry.get())); mot_lengths = [len(i) for i in motifs];
        if motifs == '':
            tkinter.messagebox.showinfo('Warm Warning', 'No, sorry. You can\'t.')
            return None
        else:
            motif_pos_dict = {}
            for seq in sequences:
                motif_pos_dict[seq] = findmotif(seq,motifs) ## list of position of the inputted motifs
    if not max_cap: max_cap = 1e6;
    new_sequences = []; generation_pool = sequences; #initialize the first pool to cross over
    for i in range(iteration):
        pairpool = list(combinations(generation_pool,2))
        pairpool_selected = [x for x in pairpool if random.random()< cross_percent] # choose a portion from the combination pool
        generation_pool=[]
        for pair in pairpool_selected:
            if len(generation_pool) >= max_cap: break ### condition to stop adding sequences
            if defined == 1:
                newpair = defined_singlept_crossover(pair,motif_pos_dict, mot_lengths)  ;### pair : ('acgt','acgt); position:[[3,2],[1,4]]
            else:
                newpair = random_singlept_crossover(pair)  ; #print(newpair) #### newpair: ['ACGGT','ACGT]
            generation_pool.append(newpair)
        generation_pool= [seq for seqpair in generation_pool for seq in seqpair]
        new_sequences.append(generation_pool)
        temp_seqs = [seq for sublist in new_sequences for seq in sublist]
        if len(temp_seqs)>max_cap:
            return temp_seqs
    return temp_seqs

def pointmutation(sequences,mutation_rate, skewrate): #mutation rate is number in a million # skew (A,C,G,T) rate
    if skewrate =='': skewrate = (1,1,1,1)
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

def findmotif(single_sequence,motif):
    for seq in sequence:
        lowest_scores = []
        for mot in motif:
            mot_scores = []; k = len(mot);
            for window in range(len(seq)-k+1):
                mot_scores.append(plainScore([seq[window:window+k],mot]))
            lowest_scores.append(mot_scores.index(min(mot_scores)))
    return lowest_scores   # [3,2,4,..]

###################################################### GUI Interface
######################################################## Cross over GUI
def cross_options():
    seq, _, _ = entry_processor(Dna_entry, k_entry, N_entry); upload = 0;
    if seq == None and seq_list.get() == None:   #### None or ''
        tkinter.messagebox.showinfo('Opps?','There is something fishy going on ...')
        return None
    elif seq_list.get()!= None:
        upload = 1;
    windowcross = Tk(); windowcross.title('Crossover Options'); windowcross.geometry('1080x700');



    option1 = Button(windowcross,text= 'Random\nsingle point',width = 15, height = 2,fg='red',
        command= lambda:gui_crossover(textbox,crosspercent_entry, gen_entry, maxcap_entry, mutation_entry,mutskew_entry,0,upload));
    option1.grid(row = 0, column = 0, padx = 5,pady = 5,rowspan=2)
    option2 = Button(windowcross,text= 'Random\nMulti point',width = 15, height = 2,fg='orange',command= lambda: None);
    option2.grid(row = 0, column = 1,padx=5, pady = 5,rowspan=2)
    option3 = Button(windowcross, text='Defined\nMotif', width=12, height=2, fg='violet',
                     command=lambda: gui_crossover(textbox,crosspercent_entry, gen_entry, maxcap_entry, mutation_entry,mutskew_entry,0,upload));
    option3.grid(row=0, column=2, padx=5, pady=5,rowspan=2)
    Label(windowcross, text="Motifs:").grid(row=0, column=3)

    motif_entry = Entry(windowcross); motif_entry.bind("<Return>"); motif_entry.grid(row=1, column=3);
    Button(windowcross, text='G-quad', command=g_gui,fg='purple').grid(row=1, column=4)
    Button(windowcross, text='Reset Space', width=17, height=2, fg='black',
                          command=lambda:reset_space(textbox)).grid(row=1, column=5)
    label = Label(windowcross,text="We got the sequence, so now these:"); label.grid(row=2,column=1)
    f = font.Font(label,label.cget("font")); f.configure(slant='italic'); label.configure(font=f);

    Label(windowcross, text = "Rate of Cross-over : (0 - 1)").grid(row=3,column = 0)
    crosspercent_entry = Entry(windowcross); crosspercent_entry.bind("<Return>"); crosspercent_entry.grid(row=3, column=1);

    Label(windowcross, text = "No. of Generations: ").grid(row=3,column = 2)
    gen_entry = Entry(windowcross); gen_entry.bind("<Return>"); gen_entry.grid(row=3, column=3);

    Label(windowcross, text="Number of sequence generated (default:1000)").grid(row=3, column=4)
    maxcap_entry = Entry(windowcross); maxcap_entry.bind("<Return>"); maxcap_entry.grid(row=3, column=5);

    Label(windowcross, text = "Rate of Mutation: (/million)").grid(row=4,column = 0)
    mutation_entry = Entry(windowcross);  mutation_entry.bind("<Return>");  mutation_entry.grid(row=4, column=1);

    Label(windowcross, text = "Skewrate: (integer) (A,C,G,T)").grid(row=4,column = 2)
    mutskew_entry = Entry(windowcross); mutskew_entry.bind("<Return>"); mutskew_entry.grid(row=4, column=3);

    textbox = Text(windowcross, height=30, width=130); textbox.grid(row=5, column=0, columnspan=6)

    ## G-quad Crossover

def g_gui():
    windowg = Tk();
    Label(windowg, text = " G-quad Sequences: ").pack()
    gseq_entry = Entry(windowg, width=80); gseq_entry.bind("<Return>"); gseq_entry.pack()
    seq = ['ACGTAGTAGATATA', 'CACAGATAATATA']
    checknames = []  # Positions for G-quadruplexes
    Label(windowg,text = 'Then, select \'em Gs').pack()
    for x in range(1, len(seq[0])):  # Generate Variables
        exec("gcheck_bool%s = IntVar()" % (x))
        checknames.append("gcheck_bool{0}".format(x))
    checkbuttonframe = Frame(windowg); checkbuttonframe.pack(side=TOP)
    for i in range(0, len(seq[0]) - 1):
        exec("Checkbutton(checkbuttonframe, text = 'P%d'%(i+1), variable=locals()[checknames[i]]).pack(side=LEFT)")

    gbool_array = ""
    button_frame = Frame(windowg); button_frame.pack(side=LEFT)
    Button(button_frame,text = 'Confirm Option',command=lambda:getgcheckbools(checknames)).pack(side = TOP,pady=5)
    Button(button_frame,text = 'Generate crossover sequence',command = gcrossover).pack(side = TOP, pady=5)
    textboxg = Text(windowg, height=30, width=80); textboxg.insert(END, 'Input Example: \n' +str(seq))
    textboxg.pack(side = BOTTOM, anchor=S)

    windowg.mainloop()
def gcrossover():
    sequences = gseq_entry.get()
    motif_blocks = sequence_blockseparate(sequences)

    return textboxg.insert(END, sequences)
def getgcheckbools(checknames):
    gcheckbool_array = []
    for i in range(len(checknames)):
        gcheckbool_array.append(exec('gcheck_bool%s.get()' % (i + 1)))
    gbool_array.set(gcheckbool_array)

def markbreaks(gcheck_bools):  ## return partition points by given array
    breakpoints = [0]
    for i in range(len(gcheck_bools) - 1):
        if gcheck_bools[i] != gcheck_bools[i + 1]:
            breakpoints.append(i + 1)
    return breakpoints  # return start-point index of those point of differences, added first position

def sequence_blockseparate(sequences, gcheck_bools):  ### Input: List of sequences of t-length, boolean array.
    partitionpoints = markbreaks(gcheck_bools);
    sequences_motif_blocks = [];
    for i in range(len(sequences)):
        thisseq = sequences[i];
        thisseq_motifs = []
        for i in range(len(partitionpoints)):
            if gcheck_bools[partitionpoints[i]] == 0:
                thisseq_motifs.append(thisseq[partitionpoints[i]:partitionpoints[i + 1]])
        sequences_motif_blocks.append(thisseq_motifs)
    return sequences_motif_blocks

def gui_crossover(textbox,crosspercent_entry, gen_entry, maxcap_entry, mutation_entry, mutskew_entry,defined,upload):
    if upload == 0:
        seq = ast.literal_eval(Dna_entry.get());
    elif upload ==1:
        seq = seq_list.get()
    if gquad_bool:
        if cross_percent.get() != '':
            cross_percent = float(crosspercent_entry.get())

    else:
        cross_percent = float(crosspercent_entry.get());
        iteration = int(gen_entry.get());
        if maxcap_entry.get() != None:
            max_cap = 1000;
        else:
            max_cap = int(maxcap_entry.get())
        new_sequence = crossover(seq, cross_percent, iteration, max_cap, defined);
        if mutation_entry.get() != '':
            mutation_rate = float(mutation_entry.get());  # skewrate = mutskew_entry.get()
            mutation_skew = mutskew_entry.get()
            if mutskew_entry.get() != '':
                new_sequence = pointmutation(new_sequence, mutation_rate, list(mutation_skew))
        textbox.insert(END, 'produced: ' + str(len(new_sequence)) + ' sequences.\n');
        for x in new_sequence:
            textbox.insert(END, x + '\n')

##################### Entry
from tkinter import *
root = Tk()
root.geometry("950x700"); root.title('Holy Moti');
root.iconbitmap(r"C:\Users\Mo\Dropbox\notes\Bioinformatics\Motif Search GUI\MM_icon.ico")
#root.iconbitmap(r"C:\malib")
Label(root, text= "Sequence: (Format: Python 'List')").grid(row=0); Dna_entry = Entry(root);
Dna_entry.bind("<Return>");Dna_entry.grid(row = 0, column = 1);

Label(root, text= "Motif Length Up to: (Format: Integer < Kmax)").grid(row=1);k_entry = Entry(root);
k_entry.bind("<Return>"); k_entry.grid(row = 1, column = 1);

Label(root, text= "No. of Reiteration (Default: 100): (Format: Integer)").grid(row=2); N_entry = Entry(root);
N_entry.bind("<Return>"); N_entry.grid(row = 2, column = 1);

others_label = Label(root, text = 'Other Functions'); others_label.grid(row = 2,column=3);  ##  Make underlined text
f = font.Font(others_label,others_label.cget("font"));
f.configure(underline=True, slant = 'italic');others_label.configure(font=f);

scores = Variable() ; entropy = Variable(); BestMotifs = Variable();

from tkinter.filedialog import askopenfilename
################################### Import Functions
import csv
def fixup_thelist(mylist):
    min = len(mylist[0]); lengthlist = [min]; cleaned_list = []
    for x in (x for x in mylist if 'N' not in x):
        cleaned_list.append(x[:x.find(' ')]);
        lengthlist.append(len(x))
        if lengthlist[-1]<=lengthlist[-2]:
            min = lengthlist[-1]
    return [i[0:min-1] for i in cleaned_list]
def csv2list(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = [];
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return fixup_thelist(my_seq)

def file_upload():           ############# Import csv or fastq
    Tk().withdraw();filename = askopenfilename(); T.insert(END, '\n\nImporting...');
    T.insert(END, '\nImported file:\n '+ filename +'\n')
    if filename.endswith('.fastq'): T.insert(END, '\nThis function is not ready T_T\'\n')
    if filename.endswith('.csv'):
        my_seq = csv2list(filename)
        T.insert(END, 'The file you uploaded contains: ' + str(len(my_seq))+' processed sequences.\n')
    else:
        T.insert(END,'Help!! Something has gone horribly wrong :< ')
        return None
    return seq_list.set(my_seq)
############################### Buttons
seq_list = Variable(); status_var = StringVar(); status_var = 'Search Awaits';
upload_button_fastq = Button(root, text = 'Upload Fastq(.fastq!)', command=file_upload).grid(row=0, column = 2)
upload_button_csv = Button(root, text = 'Upload CSV(.csv!)', fg='navy blue',command=file_upload).grid(row=1, column = 2)

###################### Motifs Button
frame1 = Frame(root); frame1.grid(row = 4, column = 0);
compute_button1 = Button(frame1, text = 'Randomized\nMotif Search', width = 15, height = 2, fg = 'blue',
                        command = lambda:Status_search(1)); compute_button1.pack(side = LEFT);
frame2 = Frame(root); frame2.grid(row = 4, column = 1);
compute_button2 = Button(frame2, text = 'Gibb\'s Sampling\nMotif Search', width = 15, height = 2, fg = 'green',
                        command = lambda:Status_search(2)); compute_button2.pack()
frame3 = Frame(root); frame3.grid(row = 4, column = 2);
compute_button3 = Button(frame3, text = 'Greedy\nMotif Search', width = 15, height = 2, fg = 'purple',
                        command = lambda:Status_search(3)); compute_button3.pack()
frame4 = Frame(root); frame4.grid(row = 4, column = 3);
compute_button4 = Button(frame4, text = 'Crossover and Mutation', width = 20, height = 2, fg = 'orange',
                        command = cross_options); compute_button4.pack()

def reset():
    T.delete('1.0', END);
    return label_status.config(text='cleared') ,seq_list.set('')
def reset_space(textbox):
    return textbox.delete('1.0',END);
frame5 = Frame(root, width = 20, height = 4).grid(row = 5, column = 2);
reset_button = Button(frame5, text = 'Reset Space &  Import', width = 17, height = 1, fg = 'black', pady = 3,
                        command = reset).grid(row = 5, column = 2,sticky='S');
exit_button = Button(root, text = 'Exit', width = 10, height = 1, fg = 'red',
                        command = exit).grid(row = 5, column = 1,sticky='S');
graph_button = Button(root, text = 'Show Graph', width = 10, height = 1, fg = 'black', pady = 15, anchor = 's',
                        command = lambda:open_graph(scores,entropy)).grid(row = 5, column = 0,sticky='S');

label_status = Label(text = 'Status: '+ status_var);label_status.grid(row = 6, column = 1);
#### Textbox
T = Text(root, height=30, width=117); T.insert(END, 'Input Example: \n' +str(Dna));## T.insert(END, str(Dna));
T.grid(row = 7, column = 0, columnspan=4, rowspan=2, padx=2, pady=3)

################################################################ Basic GUI functions
def Status_search(number):
    Status_running(); upload_status = 0;
    if seq_list.get() != '': upload_status =1; #print('sequence is: '+ str(seq_list.get()))
    print('Upload is: '+str(upload_status))
    if number == 1:
        if display_motif_Random(upload_status) == None: return None
        _, (scores, entropy, motifs) = display_motif_Random(upload_status); BestMotifs.set(motifs)
    if number == 2:
        if display_motif_Gibbs(upload_status) == None: return None
        _, (scores, entropy, motifs)= display_motif_Gibbs(upload_status); BestMotifs.set(motifs)
    if number == 3:
        if display_motif_Greedy(upload_status) == None: return None
        _, (scores, entropy, motifs) = display_motif_Greedy(upload_status); BestMotifs.set(motifs)
    return label_status.config(text='Status: Done! I''m ready.')
def Status_running():
    status_var = 'Searching for the bottom of your patience...'
    return status_var

################################################################## Graph Plotting Window
###### creating variable x,y output variable
k_motif = Variable()

def exit_graph(window):
    window.destroy();
    T.insert(END, '\nGraph Closed \n\n');
def display_sequence(seq_T,entry):
    k_input = entry.get(); int_k = int(k_input);
    if not k_input or int(k_input)> int(k_entry.get()) or int(k_input)< 2:
        tkinter.messagebox.showinfo('Wait?!','Sorry, you cannot do that :(')
    else:
        motifs_k = list(BestMotifs.get()); motif = motifs_k[int_k-2];
        [seq_T.insert(END, 'The Motifs with length k = ' +k_input+' are: \n'+ str(motif)+'\n'+
                      'Consensus Motif: '+ str(Consensus(motif))+'\n\n')];
        return k_motif.set(motifs_k[int_k-2])

def open_graph(scores, entropy):
    scores = list(scores.get()); entropy = list(entropy.get())
    T.insert(END, '\n\nGraph Requested')
    if not scores or not entropy:
        label_status.config(text='Status: Except it is not working. See why.');
        T.insert(END, '\nError: No data found.');
    else:
        window = Tk(); window.title("Graphs"); window.geometry('820x740')
        canvas = Canvas(window, width= 800, height=510); canvas.pack(expand =1)
# data
        num = len(scores)
        x = np.linspace(2, num+1, num);
        y1 = scores; y2 = entropy;
# figure
        fig = plt.figure(figsize=(8, 5))
        plt.subplot(2,1,1); plt.plot(x,y1,'b-*');
        plt.title('Motif Search Results'); plt.ylabel('Motif Scores (/base)')
        plt.xticks(np.arange(2,int(k_entry.get())+1,step=1))

        plt.subplot(2,1,2); plt.plot(x,y2,'r-*');
        plt.xlabel('K of k-mers');plt.ylabel('Motif Entropy (/base)')
        plt.xticks(np.arange(2, int(k_entry.get()) + 1, step=1))
# Canvas
        fig_photo = draw_figure(canvas, fig, loc=(5, 3))
## Textbox
        window_frame2 = Frame(window, width=650, height=120);
        window_frame2.pack(side=BOTTOM, expand=1);
        seq_T = Text(window_frame2, height=8, width=98);
        seq_T.pack(expand=1);
        seq_T.insert(END, 'The motifs resulted are as follow...\n\n')
        window_frame3 = Frame(window,width=300, height = 20); window_frame3.pack(side=BOTTOM)
# Entry
        window_frame1 = Frame(window, width =800, height = 20); window_frame1.pack()
        label_seq = Label(window_frame1, text=u"Enter the \'K\' of K-mers to be examined: (<K max)");label_seq.pack(side=LEFT)
        k_input = Entry(window_frame1);k_input.bind("<Return>"); k_input.pack(side=LEFT);
        enter_button = Button(window_frame1, text = 'Enter',
                              command=lambda:display_sequence(seq_T,k_input)); enter_button.pack();
        cluster_button = Button(window_frame3, text = "Display Motif Dendrogram",command=lambda:open_dendrogram(k_motif,int(k_input.get())));
        cluster_button.pack(side=RIGHT)
        window.mainloop()
###################################### Dendrogram with Levenstein Distance of motifs
from scipy.cluster import hierarchy
import pandas as pd
def open_dendrogram(motif, k_s):
    if k_s<3:
        tkinter.messagebox.showinfo('Actually...','Actually, nobody does that.')
        return None
    else:
        window2 = Tk();
        window2.title("Cluster");
        window2.geometry('820x540')
        canvas1 = Canvas(window2, width=800, height=720);
        canvas1.pack(expand=1)

        motif = list(k_motif.get());  print(motif);
        dist_matrix = levenstein_distances(motif, (1, 1, 1))
        fig = plt.figure(figsize=(8.5, 4.3), dpi=100, facecolor='w', edgecolor='k')
        Z = hierarchy.linkage(dist_matrix);
        plt.xlabel('Selected Motifs');
        plt.ylabel('Levinstein Edit Distances (No. of Edits)');
        plt.title('Motif Candidate Cluster')
        hierarchy.dendrogram(Z, leaf_rotation=10, leaf_font_size=7, labels=motif)
        fig_x, fig_y = 1, 1
        fig_photo = draw_figure(canvas1, fig, loc=(fig_x, fig_y))
        ### motif logo
        # print(Entropy4Logo(motif))
        draw_logo2(Entropy4Logo(motif))
        window2.mainloop()

def iterative_levenshtein(s, t, costs=(2, 2, 2)):
    rows = len(s) + 1
    cols = len(t) + 1
    deletes, inserts, substitutes = costs
    dist = [[0 for x in range(cols)] for x in range(rows)]
    ### Setup 1st row and 1st column
    for row in range(1, rows):
        dist[row][0] = row * deletes
    for col in range(1, cols):
        dist[0][col] = col * inserts

    for col in range(1, cols):
        for row in range(1, rows):
            if s[row - 1] == t[col - 1]:
                cost = 0
            else:
                cost = substitutes
            dist[row][col] = min(dist[row - 1][col] + deletes,
                                 dist[row][col - 1] + inserts,
                                 dist[row - 1][col - 1] + cost)  # substitution
    return dist[row][col]
def levenstein_distances(consensus_motifs, cost): ### matrix of motif
    motif_dict = {}
    for i in consensus_motifs:
        try: motif_dict[i] += 1
        except KeyError: motif_dict[i] = 1
    motif_list = list(motif_dict.keys())
    #### Matrix for levenstein_distances
    cols = len(motif_list); rows = cols;
    distance_matrix = np.zeros((rows,cols))
    for i in range(rows):
        distance_matrix[i][i] = 0;
        for ii in range(i+1,cols):
            distance_matrix[i][ii] = \
                iterative_levenshtein(motif_list[i],motif_list[ii],costs =cost)
            ################# flip matrix
            distance_matrix[ii][i] = distance_matrix[i][ii]
    return distance_matrix

######################################################## Matlab plot library
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg

def draw_figure(canvas, figure, loc=(0, 0)):

    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = PhotoImage(master=canvas, width=figure_w, height=figure_h)
# Position: convert from top-left anchor to center anchor
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)
# Unfortunately, there's no accessor for the pointer to the native renderer
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)

    return photo
#################################################### Drawing logo Class from GIST

import seaborn
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import numpy as np

COLOR_SCHEME = {'G': 'orange',
                'A': 'red',
                'C': 'blue',
                'T': 'darkgreen'}

BASES = list(COLOR_SCHEME.keys())

class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def draw_logo2(all_scores, fontfamily='Arial', size=80):
    if fontfamily == 'xkcd':
        plt.xkcd()
    else:
        mpl.rcParams['font.family'] = fontfamily
    window = Tk();     window.title("Cluster");   window.geometry('1120x300')
    canvas2 = Canvas(window, width=1100, height=420);   canvas2.pack(expand=1, side=LEFT)

    fig, ax = plt.subplots(figsize=(len(all_scores), 2.5))

    font = FontProperties()
    font.set_size(size)
    font.set_weight('bold')

    # font.set_family(fontfamily)

    ax.set_xticks(range(1, len(all_scores) + 1))
    ax.set_yticks(range(0, 3))
    ax.set_xticklabels(range(1, len(all_scores) + 1), rotation=90)
    ax.set_yticklabels(np.arange(0, 3, 1))
    seaborn.despine(ax=ax, trim=True)

    trans_offset = transforms.offset_copy(ax.transData,
                                          fig=fig,
                                          x=1,
                                          y=0,
                                          units='dots')

    for index, scores in enumerate(all_scores):
        yshift = 0
        for base, score in scores:
            txt = ax.text(index + 1,
                          0,
                          base,
                          transform=trans_offset,
                          fontsize=80,
                          color=COLOR_SCHEME[base],
                          ha='center',
                          fontproperties=font,

                          )
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height * score
            trans_offset = transforms.offset_copy(txt._transform,
                                                  fig=fig,
                                                  y=yshift,
                                                  units='points')
        trans_offset = transforms.offset_copy(ax.transData,
                                              fig=fig,
                                              x=1,
                                              y=0,
                                              units='points')
    fig_photo2 = draw_figure(canvas2, fig)
    window.mainloop()
root.mainloop()
