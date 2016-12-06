'''
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

interpret the fasta files
'''


import scoring_matrices as score
import pickle

def convert_seq(seq):
    '''takes in a string sequence and converts to number representation
    of amino acids according to amino_dic'''
    seq_list =[]
    for amino_acid in seq:
        if amino_acid in score.amino_dic:
            seq_list += [score.amino_dic[amino_acid]]
    return seq_list

def read_fasta(filename):
    '''reads in the fasta files and creates the desired data structures'''
    names = []
    seqids = []
    descriptions = []
    sequences = []
    for line in open(filename) :
        if line[0]=='>':
            #this is the fasta indication of a new sequence
            #these may have sequence id, name, and description
            line=line[1:].strip()
            #for now just this need to edit when line gets complex!!
            #TODO
            descriptions += ['']
            names += [line]
            seqids += [line]
            sequences += [[]]
        elif line!='' and line!='\n' :
            # if sequence data add to my sequence
            sequences[-1]+=convert_seq(line.strip())
    new_file = filename + "_seqs_data.pik"
    with open(new_file, 'wb') as f:
        pickle.dump([sequences,seqids,names,descriptions],f,-1)
    return sequences,seqids,names,descriptions

    
