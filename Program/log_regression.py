'''
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

finally perform logistic regression on our data!
'''


import numpy as np
import pickle
import matplotlib.pyplot as plt
from sklearn import linear_model

def log_reg(filename, penalty='l2', reg=1.0):
    '''grabs the data file and performs a 3 group logistic
    regression!'''    
    
    data_file = filename + "_data.pik"
    with open(data_file, 'rb') as f:
        data = pickle.load(f)
    
    np.random.shuffle(data)
    
    data_test, data_train = data[:len(data)/4,:], data[len(data)/4:,:]
    
    x_train = data_train[:,1:]
    y_train = data_train[:,:1]
    x_test = data_test[:,1:]
    y_test = data_test[:,:1]
    
    #trying out logistic regression...  
    
    accuracy=[] 
    for regl in np.linspace(1e-4,5,20):
        logreg = linear_model.LogisticRegression(C=regl)

        # we create an instance of Neighbours Classifier and fit the data.
        logreg.fit(x_train, y_train)

        result = logreg.predict(x_test)
        accuracy += [(result==y_test.flatten()).sum()/len(result)]

    plt.plot(np.linspace(1e-4,5,20), accuracy)
    
    #best reg is like  3.4 but any where in this range was fine
    
    accuracy=[] 
    for regl in np.linspace(1e-4,5,20):
        logreg = linear_model.LogisticRegression(C=regl,penalty='l1')

        # we create an instance of Neighbours Classifier and fit the data.
        logreg.fit(x_train, y_train)

        result = logreg.predict(x_test)
        accuracy += [(result==y_test.flatten()).sum()/len(result)]

    plt.plot(np.linspace(1e-4,5,20), accuracy)
   
   #L2 performs better than L1!