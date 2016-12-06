'''
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

Making the data files as numpy matrices
now we have how to represent each amino acid!
'''

import numpy as np
import pickle
import scoring_matrices as score

def make_matrices(filename):
    '''using the structure, sequence and the scores made in get_scores
    to create the matrices'''
    
    score_file = filename + "_score_data.pik"
    with open(score_file, 'rb') as f:
        current_loc, one_off, two_off, three_off = pickle.load(f)
    
    struct_file = filename + "_struct_data.pik"
    with open(struct_file, 'rb') as f:
        structures, sequences = pickle.load(f)

    data = []
    j=0
    for sequence in sequences:
        structure = structures[j]
        j+=1
        seq_length = len(sequence)
        for i in range(seq_length):
            if sequence[i] in score.amino_acid_list:
                data_to_add = [structure[i],current_loc[sequence[i]]]
                if i>0:
                    if sequence[i-1] in score.amino_acid_list:
                        data_to_add += [one_off[sequence[i-1]]]
                    else:
                        data_to_add += [one_off[sequence[i]]]
                else:
                    data_to_add += [one_off[sequence[i]]]

                if i>1:
                    if sequence[i-2] in score.amino_acid_list:
                        data_to_add += [two_off[sequence[i-2]]]
                    else:
                        data_to_add += [two_off[sequence[i]]]
                else:
                    data_to_add += [two_off[sequence[i]]]

                if i>2:
                    if sequence[i-3] in score.amino_acid_list:
                        data_to_add += [three_off[sequence[i-3]]]
                    else:
                        data_to_add += [three_off[sequence[i]]]
                else:
                    data_to_add += [three_off[sequence[i]]]

                if i<seq_length-1:
                    if sequence[i+1] in score.amino_acid_list:
                        data_to_add += [one_off[sequence[i+1]]]
                    else:
                        data_to_add += [one_off[sequence[i]]]
                else:
                    data_to_add += [one_off[sequence[i]]]

                if i<seq_length-2:
                    if sequence[i+2] in score.amino_acid_list:
                        data_to_add += [two_off[sequence[i+2]]]
                    else:
                        data_to_add += [two_off[sequence[i]]]
                else:
                    data_to_add += [two_off[sequence[i]]]

                if i<seq_length-3:
                    if sequence[i+3] in score.amino_acid_list:
                        data_to_add += [three_off[sequence[i+3]]]
                    else:
                        data_to_add += [three_off[sequence[i]]]
                else:
                    data_to_add += [three_off[sequence[i]]]
                                    

                data+=[data_to_add]
    data = np.array(data)
    data_file = filename + "_data.pik"
    with open(data_file, 'wb') as f:
        pickle.dump(data,f,-1)
    return data
