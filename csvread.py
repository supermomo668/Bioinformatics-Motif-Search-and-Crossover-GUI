import csv
def csv2list(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = [];
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return my_seq

mylist = csv2list('myseqs.csv')
print(mylist[0])

def even_thelist(mylist):
    min = len(mylist[0]); lengthlist = [min]; cleaned_list = []
    for x in (x for x in mylist if 'N' not in mylist):
        cleaned_list.append(x);
        lengthlist.append(len(x))
        if lengthlist[-1]<=lengthlist[-2]:
            min = lengthlist[-1]
        print(min)
        print(lengthlist)
        print (x)
    newlist = []
    return [newlist.append(i[0:min-1]) for i in cleaned_list]