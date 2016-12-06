'''
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

Method to get the scores of each amino acid for the current position,
and one, two or three amino acid away
'''
import scoring_matrices as score

def get_scored(seqs,structures):
    '''takes in our sequence and structure and gives each base
    its score, negative indicating beta behavior, positive alpha
    and near zero coil'''

    current_loc = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
    one_off ={'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
    two_off ={'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
    three_off = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}

    #lets normalize the data!
    num_proteins = len(seqs)
    for amino_acid in score.amino_acid_list:
        current_loc[amino_acid] = current_loc[amino_acid]/(1.0*num_proteins)
        one_off[amino_acid] = one_off[amino_acid]/(num_proteins*2.0)
        two_off[amino_acid] = two_off[amino_acid]/(num_proteins*2.0)
        three_off[amino_acid] = three_off[amino_acid]/(num_proteins*2.0)
