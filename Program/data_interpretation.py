'''
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

interpret the fasta files
'''


import get_scores
import pickle

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
            if (line[:4].lower() not in seqids):
                if len(line)<9:
                    descriptions += ['']
                    names += [line]
                    seqids += [line]
                    sequences += [[]]
                else:
                    descriptions += ['']
                    names += [line]
                    seqids += [line[:4].lower()]
                    sequences += [[]] 
        elif line!='' and line!='\n' :
            # if sequence data add to my sequence
            sequences[-1]+=get_scores.convert_seq(line.strip())
    new_file = filename + "_seqs_data.pik"
    with open(new_file, 'wb') as f:
        pickle.dump([sequences,seqids,names,descriptions],f,-1)
    return sequences,seqids,names,descriptions

    
