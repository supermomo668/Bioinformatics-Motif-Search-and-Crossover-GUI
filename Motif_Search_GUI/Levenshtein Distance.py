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
################################ Function Test
motif = ['ACGTTTT','AGCTTTA','ACGAATACA','ACAGAGCCAATA','ACCGAGTTTAATA','ACGTTTA']


################################
import numpy as np
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

dist_matrix = levenstein_distances(motif,(2,2,2))
print(dist_matrix)

#dist_matrix_unroll = [item for sublist in dist_matrix for item in sublist]

from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

Z = hierarchy.linkage(dist_matrix);
plt.xlabel('Selected Motifs'); plt.ylabel('Levinstein Edit Distances (No. of Edits)'); plt.title('Motif Candidate Cluster')
hierarchy.dendrogram(Z, leaf_rotation=10, leaf_font_size=7,labels=motif);
plt.show()
