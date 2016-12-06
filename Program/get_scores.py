'''
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

Method to get the scores of each amino acid for the current position,
and one, two or three amino acid away
'''
import scoring_matrices as score
import pickle

def get_scored(filename):
    '''takes in our sequence and structure and gives each base
    its score, negative indicating beta behavior, positive alpha
    and near zero coil'''

    current_loc = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
    one_off ={'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
    two_off ={'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
    three_off = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}

    struct_file = filename + "_struct_data.pik"
    with open(struct_file, 'rb') as f:
        structures, sequences = pickle.load(f)

    j = 0
    num_amino_acids = 0
    for sequence in sequences:
        structure = structures[j]
        j+=1
        seq_length = len(sequence)
        for i in range(seq_length):
            num_amino_acids += 1
            score_position = structure[i]
            #looking at how the amino acids around
            #this position effect the score at this
            #position
            if sequence[i] in score.amino_acid_list:
                current_loc[sequence[i]] += score_position
            if i>0:
                if sequence[i-1] in score.amino_acid_list:
                    one_off[sequence[i-1]] += score_position
            if i>1:
                if sequence[i-2] in score.amino_acid_list:
                    two_off[sequence[i-2]] += score_position
            if i>2:
                if sequence[i-3] in score.amino_acid_list:
                    three_off[sequence[i-3]] += score_position
            if i<seq_length-1:
                if sequence[i+1] in score.amino_acid_list:
                    one_off[sequence[i+1]] += score_position
            if i<seq_length-2:
                if sequence[i+2] in score.amino_acid_list:
                    two_off[sequence[i+2]] += score_position
            if i<seq_length-3:
                if sequence[i+3] in score.amino_acid_list:
                    three_off[sequence[i+3]] += score_position   
            
        
    #lets normalize the data!
    for amino_acid in score.amino_acid_list:
        current_loc[amino_acid] = current_loc[amino_acid]/(1.0*num_amino_acids)
        one_off[amino_acid] = one_off[amino_acid]/(num_amino_acids*2.0)
        two_off[amino_acid] = two_off[amino_acid]/(num_amino_acids*2.0)
        three_off[amino_acid] = three_off[amino_acid]/(num_amino_acids*2.0)

    score_file = filename + "_score_data.pik"
    with open(score_file, 'wb') as f:
        pickle.dump([current_loc, one_off, two_off, three_off],f,-1)
        
    return current_loc, one_off, two_off, three_off


def convert_seq(seq):
    '''takes in a string sequence and converts to number representation
    of amino acids according to amino_dic'''
    seq_list =[]
    for amino_acid in seq:
        if amino_acid in score.amino_dic:
            seq_list += [score.amino_dic[amino_acid]]
    return seq_list
