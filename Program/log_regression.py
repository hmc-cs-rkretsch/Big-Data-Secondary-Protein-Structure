'''
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

finally perform logistic regression on our data!
'''


import numpy as np
import pickle

def log_reg(filename, penalty='l2', reg=1.0):
    '''grabs the data file and performs a 3 group logistic
    regression!'''    
    
    data_file = filename + "_data.pik"
    with open(data_file, 'rb') as f:
        data = pickle.load(f)
    
    np.random.shuffle(data)
    
    data_test, data_train = data[:len(data)/4,:], data[len(data)/4:,:]
    
    x_train = data_test[:,1:]
    y_train = data_test[:,:1]
    x_test = data_test[:,1:]
    y_test = data_test[:,:1]
    
    