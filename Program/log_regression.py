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
#==============================================================================
#   
#     accuracy=[] 
#     for regl in np.linspace(1e-4,5,20):
#         logreg = linear_model.LogisticRegression(C=regl)
# 
#         # we create an instance of Neighbours Classifier and fit the data.
#         logreg.fit(x_train, y_train)
# 
#         result = logreg.predict(x_test)
#         accuracy += [(result==y_test.flatten()).sum()/len(result)]
# 
#     plt.plot(np.linspace(1e-4,5,20), accuracy, label='l2 regression')
#     
#     #best reg is like  3.4 but any where in this range was fine
#     
#     accuracy=[] 
#     for regl in np.linspace(1e-4,5,20):
#         logreg = linear_model.LogisticRegression(C=regl,penalty='l1')
# 
#         # we create an instance of Neighbours Classifier and fit the data.
#         logreg.fit(x_train, y_train)
# 
#         result = logreg.predict(x_test)
#         accuracy += [(result==y_test.flatten()).sum()/len(result)]
# 
#     plt.plot(np.linspace(1e-4,5,20), accuracy, label='l1 reggression', linestyle='--')
#     plt.title("Accuracy of prediction")
#     plt.xlabel("regression coefficient")
#     plt.ylabel("Accuracy")
#     plt.legend(loc=4)
#     
#     #L2 performs better than L1!
#     #get graph for this...
#   
#==============================================================================
 
    logreg = linear_model.LogisticRegression(C=3.4)
    logreg.fit(x_train, y_train)
    result = logreg.predict(x_test) 
    accuracy = [(result==y_test.flatten()).sum()/len(result)]
    prob = logreg.predict_proba(x_test)
    coef = logreg.coef_
    print(coef)
    print(accuracy)
    
#==============================================================================
#     coef_0 = [coef[0][3],coef[0][2],coef[0][1],coef[0][0],coef[0][4],coef[0][5],coef[0][6]]
#     coef_1 = [coef[1][3],coef[1][2],coef[1][1],coef[1][0],coef[1][4],coef[1][5],coef[1][6]]
#     coef_2 = [coef[2][3],coef[2][2],coef[2][1],coef[2][0],coef[2][4],coef[2][5],coef[2][6]]
#     
#     position = [-3,-2,-1,0,1,2,3]
#     
#     plt.plot(position,coef_0,label="beta strand (-1)")
#     plt.plot(position,coef_1,label="coil (0)")
#     plt.plot(position,coef_2,label='alpha helix (1)')
#     plt.title("Regression coefficients")
#     plt.xlabel("Distance from amino acid being predicted")
#     plt.legend(bbox_to_anchor=(1.05, 1), loc=2,)
#==============================================================================
    
    print(prob)
    print(result)
   
   
    data_file = filename + "_data.pik"
    with open(data_file, 'rb') as f:
        data_full = pickle.load(f)
       
    struct_file = filename + "_struct_data.pik"
    with open(struct_file, 'rb') as f:
        structures, sequences = pickle.load(f)
   
    x_full = data_full[:,1:]
    y_full = data_full[:,:1]
    prob_full = logreg.predict_proba(x_full)
    result_full = logreg.predict(x_full)
    
    
    new_file = filename + "_seqs_data.pik"
    with open(new_file, 'rb') as f:
        sequences_2,seqids,names,descriptions = pickle.load(f)
    
    i = 0
    j = -1
    largest_errors  = [] #ORDERED
    smallest_errors = []
    accurcacies = []
    probs = []
    for seqid in seqids:
        j+=1
        length_sequence = len(sequences[j])
        result_temp = result_full[i:length_sequence+i]
        y_temp = y_full[i:length_sequence+i]
        accuracy_temp = (result_temp==y_temp.flatten()).sum()/len(result_temp)
        accurcacies += [accuracy_temp]
        probs += [prob_full[i:length_sequence+i]]
        if (len(smallest_errors)<10 or accuracy_temp<smallest_errors[0][1]) and seqid != '1unj':
            smallest_errors += [[seqid,accuracy_temp,j]]
            smallest_errors = sorted(smallest_errors,key=lambda l:l[1], reverse=True)
            if len(smallest_errors)>10:
                smallest_errors=smallest_errors[:10]
        if (len(largest_errors)<10 or accuracy_temp>largest_errors[0][1]) and seqid != '1unj':
            largest_errors += [[seqid,accuracy_temp,j]]
            largest_errors = sorted(largest_errors,key=lambda l:l[1])
            if len(largest_errors)>10:
                largest_errors=largest_errors[:10]
        i+=length_sequence
    
#==============================================================================
#     plt.hist(np.array(accurcacies)[np.isfinite(accurcacies)], np.linspace(0.0,0.8,20), alpha = 0.5)
#     plt.title("Range of accuracies for individual proteins")
#     plt.xlabel("Accuracy")
#     plt.ylabel("Frequency")
# 
#     
#==============================================================================
    print(largest_errors)
    print(smallest_errors)
    

    print('small')    
    for res in smallest_errors:
        num = res[2]
        print(seqids[num])
        print(''.join(sequences[num]))
        graph_protein_prob(probs[num],sequences[num],res[0],res[1])

    print('large')    
    for res in largest_errors:
        num = res[2]
        print(seqids[num])
        print(''.join(sequences[num]))
        graph_protein_prob(probs[num],sequences[num],res[0],res[1])



def graph_protein_prob(prob,seq,seqid,accuracy):
    '''graphs the prob of the three structures against
    the sequence'''
    x = range(len(seq))
    prob = np.array(prob)
    beta = prob[:,:1].flatten()
    coil = prob[:,1:2].flatten()
    alpha = prob[:,2:].flatten()
    plt.plot(x,beta,label='beta',linestyle='--')
    plt.plot(x,coil,label='coil')
    plt.plot(x,alpha,label='alpha',linestyle=':')
    plt.title('Secondary structure prediction '+seqid)
    plt.xlabel('Amino acid position')
    plt.ylabel('Probability')
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    ax=plt.gca()
    ax.set_xlim([0,len(seq)])
    accuracy=round(accuracy,4)
    plt.text(len(seq)+10,0.5,"accuracy = "+str(accuracy))
    plt.savefig('../Data/'+seqid+'.png',bbox_extra_artists=(lgd,),dpi=600,bbox_inches='tight')
    plt.close()
        
    


    
            
        

    
    
   