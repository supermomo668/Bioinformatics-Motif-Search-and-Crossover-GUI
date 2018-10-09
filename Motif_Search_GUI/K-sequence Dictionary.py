import csv

def fixup_thelist(mylist):
    #min = len(mylist[0]); lengthlist = [min];
    cleaned_list = []
    for x in (x for x in mylist if 'N' not in x):
        cleaned_list.append(x[:x.find(' ')]);
        #lengthlist.append(len(x))
        #if lengthlist[-1]<=lengthlist[-2]:
        #    min = lengthlist[-1]
    return cleaned_list #[i[0:min-1] for i in cleaned_list]
def csv2list(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = [];
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return fixup_thelist(my_seq)

def k_sequencedictionary(seqs,k):
    totaldict = {}
    for i in seqs:
        for ii in range(len(seqs)-k+1):
            try:
                totaldict[i[ii:ii+k]] += 1
            except KeyError:
                totaldict[i[ii:ii+k]] = 1
    return totaldict

k = 5
seq_list = csv2list('testseq_MM.csv')
output = k_sequencedictionary(seq_list, 5)
print(output)




