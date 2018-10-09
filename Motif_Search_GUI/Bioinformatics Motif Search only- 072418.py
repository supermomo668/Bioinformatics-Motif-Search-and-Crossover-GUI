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

def Pr(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob *= Profile[Text[i]][i]
    return prob

import math as m
def ProfileWithPseudocounts(Motifs):
    count, k = CountWithPseudocounts(Motifs,1);
    num_seq = len(Motifs)+4; entropy = 0;
    for i in range(k):
        for j in 'ACGT':
            count[j][i] = float(count[j][i]/num_seq)
        f = max(count[j][i] for j in 'ACGT')
        entropy += -f*m.log(f)
    return count, entropy
def Entropy4Logo(Motifs):
    count, k = CountWithPseudocounts(Motifs,0)
    num_seq = len(Motifs) + 4;
    entropy_logo = []
    for i in range(k):
        entropy_logo.append([])
        for j in 'ACGT':
            f = float(count[j][i]/num_seq)
            print(f)
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

Dna =\
    ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
     'TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
     'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
######################################################################################################### GUI
############################################################## GUI functions
import ast
import tkinter.messagebox
def upload_list():
    return seq_list.get()

def entry_processor(Dna_entry, k_entry, N_entry): ### take .entry as input
    if Dna_entry != None: Dna_disp = Dna_entry.get();seq = ast.literal_eval(Dna_disp); # turn a 'list string' into real list
    else: seq = None;
    k = int(k_entry.get());
    if not N_entry.get(): N = 100;
    else: N = int(N_entry.get());
    return [seq, k, N]

def display_results(seq, k, delta_t, k_bestmotifs,k_scores, k_entropy, searchname):
    T.insert(END, '\nYou ran: '+str(searchname) +' (k = ' + str(len(seq[0])) + ')\n| k = ' + str(k) +
             '| Number of Sequence: ' + str(len(seq)) + '| Runtime: ' + delta_t + 's\n\n')
    T.insert(END, 'Motifs:\n ' + str(k_bestmotifs) + '\n\nScores:\n' + str(k_scores) + '\n\nEntropy:\n ' + str(k_entropy) + '\n')
    scores.set(k_scores); entropy.set(k_entropy);
    return scores, entropy, k_bestmotifs

def upload_choicegate(upload):
    t0 = time.time(); print('Upload is: ' +str(upload));
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
    return display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy,searchname);  ####

def display_motif_Gibbs(upload):
    if upload_choicegate(upload) == None: return None
    seq, k, N, t0 = upload_choicegate(upload); searchname = 'Gibb\'s Randomized Sampler'
    k_bestmotifs, k_scores, k_entropy = GibbsSampler(seq, k, N)
    delta_t = str(time.time() - t0);
    return display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy,searchname);  ####

def display_motif_Greedy(upload):
    if upload_choicegate(upload) == None: return None
    seq,k,N,t0 = upload_choicegate(upload); searchname = 'Greedy Motif Search'
    k_bestmotifs, k_scores, k_entropy = GreedyMotifSearchWithPseudocounts(seq, k)
    delta_t = str(time.time() - t0);
    return display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy,searchname);##

###################################################### GUI Interface
##################### Entry
from tkinter import *
root = Tk()
root.geometry("800x700"); root.title('Holy Moti');
root.iconbitmap(r"C:\Users\Mo\Dropbox\notes\Bioinformatics\Motif Search GUI\MM_icon.ico")
Label(root, text= "Sequence: (Format: Python 'List')").grid(row=0); Dna_entry = Entry(root);
Dna_entry.bind("<Return>");Dna_entry.grid(row = 0, column = 1);

Label(root, text= "Motif Length Up to: (Format: Integer < Kmax)").grid(row=1);k_entry = Entry(root);
k_entry.bind("<Return>"); k_entry.grid(row = 1, column = 1);

Label(root, text= "No. of Reiteration (Default: 100): (Format: Integer)").grid(row=2); N_entry = Entry(root);
N_entry.bind("<Return>"); N_entry.grid(row = 2, column = 1);

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
    else: T.insert(END,'Help, something has gone horribly wrong ')
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

def reset():
    T.delete('1.0', END);
    return label_status.config(text='cleared') ,seq_list.set('')

reset_button = Button(root, text = 'Reset Space &  Import', width = 17, height = 2, fg = 'black',
                        command = reset).grid(row = 5, column = 2);
exit_button = Button(root, text = 'Exit', width = 15, height = 2, fg = 'red',
                        command = exit).grid(row = 5, column = 1);
graph_button = Button(root, text = 'Show Graph', width = 10, height = 1, fg = 'black',
                        command = lambda:open_graph(scores,entropy)).grid(row = 5, column = 0);

label_status = Label(text = 'Status: '+ status_var);label_status.grid(row = 6, column = 1);

T = Text(root, height=30, width=99); T.insert(END, 'Input Example: \n' +str(Dna));## T.insert(END, str(Dna));
T.grid(row = 7, column = 0, columnspan=3, rowspan=2, padx=2, pady=3)

################################################################ Basic GUI functions
def Status_search(number):
    Status_running(); upload_status = 0;
    if seq_list.get() !='': upload_status =1; #print('sequence is: '+ str(seq_list.get()))
    if number == 1:
        if display_motif_Random(upload_status) == None: return None
        scores, entropy, motifs = display_motif_Random(upload_status); BestMotifs.set(motifs)
    if number == 2:
        if display_motif_Gibbs(upload_status) == None: return None
        scores, entropy, motifs = display_motif_Gibbs(upload_status); BestMotifs.set(motifs)
    if number == 3:
        if display_motif_Greedy(upload_status) == None: return None
        scores, entropy, motifs = display_motif_Greedy(upload_status); BestMotifs.set(motifs)
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
    window2 = Tk(); window2.title("Cluster"); window2.geometry('820x540')
    canvas1 = Canvas(window2, width=800, height=720); canvas1.pack(expand=1)

    motif = list(k_motif.get()); print(motif);
    dist_matrix = levenstein_distances(motif, (1, 1, 1))
    fig = plt.figure(figsize=(8.5, 4.3), dpi = 100, facecolor='w', edgecolor='k')
    Z = hierarchy.linkage(dist_matrix);
    plt.xlabel('Selected Motifs'); plt.ylabel('Levinstein Edit Distances (No. of Edits)'); plt.title('Motif Candidate Cluster')
    hierarchy.dendrogram(Z, leaf_rotation=10, leaf_font_size=7,labels=motif)
    fig_x, fig_y = 1, 1
    fig_photo = draw_figure(canvas1, fig, loc=(fig_x, fig_y))
    ### motif logo
    print(Entropy4Logo(motif))
    draw_logo(Entropy4Logo(motif))
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
plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
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
        affine = affine.identity().scale(self._sx, self._sy) + affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def draw_logo(all_scores):
    window = Tk(); window.title("Cluster"); window.geometry('1120x300')
    canvas2 = Canvas(window, width=1100, height=420); canvas2.pack(expand=1,side=LEFT)
    fig = plt.figure(figsize=(3,4.5),dpi = 100, facecolor='w',edgecolor='k')
    fig.set_size_inches(len(all_scores), 2.5)
    ax = fig.add_subplot(111)
    ax.set_xticks(range(len(all_scores)))
    plt.xticks(range(len(all_scores)))

    xshift = 0
    trans_offset = transforms.offset_copy(ax.transAxes,
                                          fig=fig,
                                          x=0,
                                          y=0,
                                          units='points')
    for scores in all_scores:
        yshift = 0
        for base, score in scores:
            txt = ax.text(0,
                          0,
                          base,
                          transform=trans_offset,
                          fontsize=40,
                          color=COLOR_SCHEME[base],
                          weight='bold',
                          ha='center',
                          family='sans-serif'
                          )
            txt.set_clip_on(False)
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height * score
            trans_offset = transforms.offset_copy(txt._transform, fig=fig, y=yshift, units='points')
        xshift += window_ext.width
        trans_offset = transforms.offset_copy(ax.transAxes, fig=fig, x=xshift, units='points')

    ax.set_yticks(range(0, 3))

    seaborn.despine(ax=ax, offset=30, trim=True)
    ax.set_xticklabels(range(1, len(all_scores) + 1), rotation=90)
    ax.set_yticklabels(np.arange(0, 3, 1))
    fig_photo2 = draw_figure(canvas2, fig)
    window.mainloop()

root.mainloop()
