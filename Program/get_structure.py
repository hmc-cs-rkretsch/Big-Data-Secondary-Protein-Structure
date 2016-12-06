'''
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

get the secondary structure data!
Uses the dssp database which is a database created
by Wolfgang Kabsch and Chris Sander to standardize
the secondary strcuture assignments from the Protein
data bank. This is all put into a folder and then the
files can be read.
'''


#first run the ./get_structure file and give the list of sequence
#names as input
#now the files are all located in /DSSP/ folder.
#now need to decipher the groups

'''
We will seperate the secondary structures into three broad categories:
A: helices, includes alpha-helices and 310-helices
B: strands, includes beta-strands and isolated beta-bridges
C: coils, includes turns, bends, pi-helices, and others/not identified
'''

import scoring_matrices as score

def get_structure(seqids):
    structures = []
    for seqid in seqids:
        file = "DSSP/"+seqid+".dssp"
        i=0
        structure = []
        for line in open(file):
            i+=1
            if i>28:
                structure += [ score.sec_struct_dic[line[16]] ]
        structures += [structure]
    return structures
