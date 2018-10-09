from Bio import motifs
import seaborn
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib import transforms
import matplotlib.patheffects
import matplotlib.pyplot as plt

import numpy as np

COLOR_SCHEME = {'G': 'orange',
                'A': 'red',
                'C': 'blue',
                'T': 'darkgreen'}
bases = list(COLOR_SCHEME.keys())

def approximate_error(motif):
    """Calculate approximate error"""
    pwm = motif.pwm
    bases = list(pwm.keys())
    n = sum(motif.counts[bases[0]])
    approx_error = (len(bases)-1)/(2 * np.log(2) * n)
    return approx_error


def exact_error(motif):
    """Calculate exact error, using multinomial(na,nc,ng,nt)"""
    ## Super Slow. O(n^3)
    pwm = motif.pwm
    bases = list(pwm.keys())
    na = sum(motif.counts['A'])
    n = na
    nc = 0
    ng = 0
    nt = 0
    done = False
    exact_error = 0
    while not done:
        print (na,nc,ng,nt)
        exact_error += sum([-p*np.log2(p) for p in [na/n, nc/n, ng/n, nt/n]])
        if nt<=0:
            ## iterate inner loop
            if ng > 0:
                ## g => t
                ng = ng - 1
                nt = nt + 1
            elif nc > 0:
                ## c -> g
                nc = nc - 1;
                ng = ng + 1;
            else:
                ## a->c
                na = na - 1
                nc = nc + 1
        else:
            if ng > 0:
                ## g => t
                ng = ng - 1
                nt = nt + 1
            elif nc>0:
                ## c => g; all t -> g
                nc = nc - 1
                ng = nt + 1
                nt = 0
            elif na>0:
                ## a => c; all g,t -> c
                nc = nt + 1
                na = na - 1
                nt = 0
            else:
                done = True
    return exact_correction

def calc_info_matrix(motif, correction_type='approx'):
    """Calculate information matrix with small sample correction"""
    pwm = motif.pwm
    bases = list(pwm.keys())
    if correction_type=='approx':
        error = approximate_error(motif)
    else:
        error = exact_error(motif)
    info_matrix = [2-error+sum([pwm[b][l]*np.nan_to_num(np.log2(pwm[b][l])) for b in bases]) for l in range(0, len(motif))]
    return info_matrix

def calc_relative_information(motif, correction_type='approx'):
    """Calculate relative information matrix"""
    pwm = motif.pwm
    bases = list(pwm.keys())
    if correction_type=='approx':
        info_matrix = calc_info_matrix(motif)
    else:
        info_matrix = calc_info_matrix(motif, 'exact')
    relative_info = {base: [prob*info for prob,info in zip(pwm[base], info_matrix)]  for base in bases}
    return relative_info