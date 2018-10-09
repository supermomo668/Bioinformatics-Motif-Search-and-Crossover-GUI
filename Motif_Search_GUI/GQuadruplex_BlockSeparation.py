seq = ['ACGTAGTAGATATA','CACAGATACGTATG','GACTACTAGCTAGA'];
gcheck_bools = [0,0,1,1,0,0,1,1,1,0,0,0,1,1]

print('sequence length: '+str(len(seq[1])))
print('check bool length: '+str(len(gcheck_bools)))
def markbreaks(gcheck_bools):   ## return partition points by given array
    breakpoints = [0]
    for i in range(len(gcheck_bools)-1):
        if gcheck_bools[i] != gcheck_bools[i+1]:
            breakpoints.append(i+1)
    return breakpoints                      #return start-point index of those point of differences, added first position

def sequence_blockseparate(sequences, gcheck_bools):  ### Input: List of sequences of t-length, boolean array.
    partitionpoints = markbreaks(gcheck_bools);
    sequences_motif_blocks = [];
    for i in range(len(sequences)):
        thisseq = sequences[i]; thisseq_motifs = []
        for i in range(len(partitionpoints)):
            if gcheck_bools[partitionpoints[i]] == 0:
                thisseq_motifs.append(thisseq[partitionpoints[i]:partitionpoints[i+1]])
        sequences_motif_blocks.append(thisseq_motifs)
    return sequences_motif_blocks                                          ### Output list of motifs and their position in sequence


print(markbreaks(gcheck_bools))
print(sequence_blockseparate(seq, gcheck_bools))

def motifs_join_4mutation(sequences):
    return [''.join(i) for i in sequences]

print(motifs_join_4mutation(sequence_blockseparate(seq, gcheck_bools)))