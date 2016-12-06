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
A: (1) helices, includes alpha-helices and 310-helices
B: (-1) strands, includes beta-strands and isolated beta-bridges
C: (0) coils, includes turns, bends, pi-helices, and others/not identified
'''

import scoring_matrices as score
import pickle
import os.path

def check_struct_files(filename):
    '''checks if get_structure.py was run correctly
    returns all the structure files that were not added.
    we can then try and manually find them (perhaps name
    differs) or not'''
    new_file = filename + "_seqs_data.pik"
    with open(new_file, 'rb') as f:
        sequences,seqids,names,descriptions =pickle.load(f)
    missing = []
    for seqid in seqids:
        file = "DSSP/"+seqid+".dssp"
        if not os.path.isfile(file):
            missing += [seqid]
    return missing

def get_structure(filename):
    '''get the structure from the file, represents as 1, -1, or 0'''
    new_file = filename + "_seqs_data.pik"
    with open(new_file, 'rb') as f:
        sequences,seqids,names,descriptions =pickle.load(f)
    structures = []
    sequences = []
    for seqid in seqids:
        file = "DSSP/"+seqid+".dssp"
        i=0
        structure = []
        sequence = []
        for line in open(file):
            i+=1
            if i>28:
                structure += [ score.sec_struct_dic[line[16]] ]
                sequence += [ line[13] ]
        structures += [structure]
        sequences += [sequence]
    struct_file = filename + "_struct_data.pik"
    with open(struct_file, 'wb') as f:
        pickle.dump([structures,sequences],f,-1)
    return structures


def check_structure(filename):
    '''checks that the strucutre data and the sequence data
    match in length. Returns any that do not.'''
    struct_file = filename + "_struct_data.pik"
    with open(struct_file, 'rb') as f:
        structures, sequences = pickle.load(f)
    incorrect = []
    for seqid in range(len(sequences)):
        if len(sequences[seqid])!=len(structures[seqid]):
               incorrect += [seqid]
               incorrect += [len(sequences[seqid])-len(structures[seqid])]
    return incorrect
