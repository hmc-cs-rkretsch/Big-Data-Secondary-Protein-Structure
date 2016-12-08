# -*- coding: utf-8 -*-
"""
Author: Rachael Kretsch
Big Data Final Project
Secondary Protein Structure

analyse the literature data!!
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model

#==============================================================================
#  
# struct_file = '../Data/color_fasta.txt' + "_struct_data.pik"
# with open(struct_file, 'rb') as f:
#     structures, sequences = pickle.load(f)   
#      
# new_file = '../Data/color_fasta.txt' + "_seqs_data.pik"
# with open(new_file, 'rb') as f:
#     sequences_2,seqids,names,descriptions = pickle.load(f)
#     
# good_seqids=['1bfp','3ekh','3ned','4k3g','3cfc',
#         '1xkh','2wht','4w6b','4xvp','3dqh',
#         '1bgp','4q7t','4qgw','5h88','4l1s',
#         '5h89','3s0f','4q9w','3rwt','5hzo']
# 
# s2D_accuracies = {}
# s2D_accuracy = 0
# 
# j=-1
# 
# s2D_predictions = []
# for seqid in seqids:
#     j+=1
#     if seqid in good_seqids:
#         struct = structures[j]
#         i=0
#         prediction = []
#         for line in open('../Data/S2D/'+seqid + '_s2d_out.txt'):
#             if line[0]!='>' and line[0]!='#':
#                 if i<11:
#                     pred = line[40]
#                 elif i<101:
#                     pred= line[41]
#                 else:
#                     pred = line[42]
#                 if pred=='H':
#                     prediction+=[1]
#                 elif pred=='E':
#                     prediction+=[-1]
#                 else:
#                     prediction+=[0]
#             i+=1
#         print(seqid)
#         x = range(len(prediction))
#         beta = []
#         alpha = []
#         coil = []
#         for amino in prediction:
#             if amino == -1:
#                 beta += [1]
#                 alpha += [0]
#                 coil += [0]
#             elif amino == 1:
#                 beta += [0]
#                 alpha += [1]
#                 coil += [0]
#             else:
#                 beta += [0]
#                 alpha += [0]
#                 coil += [1]
#         plt.scatter(x,beta,label='beta',marker = 'o',color='blue')
#         plt.scatter(x,coil,label='coil',marker='x', color='green')
#         plt.scatter(x,alpha,label='alpha',color='red')
#         plt.title('Secondary structure prediction s2D '+seqid)
#         plt.xlabel('Amino acid position')
#         plt.ylabel('Probability')
#         lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
#         ax=plt.gca()
#         fig = plt.gcf()
#         fig.set_size_inches
#         ax.set_xlim([0,len(prediction)])
#         ax.set_ylim([0.9,1.1])
#         plt.savefig('../Data/S2D/'+seqid+'_actual.png',bbox_extra_artists=(lgd,),dpi=600,bbox_inches='tight')
#         plt.close()               
#         s2D_predictions+=[[beta,coil,alpha]]
#         struct=struct[:len(prediction)]
#         acc=(np.array(prediction)==np.array(struct)).sum()/len(prediction)
#         s2D_accuracy+=acc
#         s2D_accuracies[seqid]=acc
# 
# s2D_accuracy=s2D_accuracy/len(good_seqids)
# print("accuracy s2D: "+str(s2D_accuracy))
# 
# 
# SOPM_accuracies = {}
# SOPM_accuracy = 0
# 
# j=-1
# 
# SOPM_predictions=[]
# for seqid in seqids:
#     j+=1
#     if seqid in good_seqids:
#         struct = structures[j]
#         prediction = []
#         for line in open('../Data/SOPM/'+seqid + '.sopm.txt'):
#             if line[0]in ['H','C','E'] and len(prediction)<len(struct):
#                 pred = line[0]
#                 if pred=='H':
#                     prediction+=[1]
#                 elif pred=='E':
#                     prediction+=[-1]
#                 else:
#                     prediction+=[0]
#         print(seqid)
#         
#         x = range(len(prediction))
#         beta = []
#         alpha = []
#         coil = []
#         for amino in prediction:
#             if amino == -1:
#                 beta += [1]
#                 alpha += [0]
#                 coil += [0]
#             elif amino == 1:
#                 beta += [0]
#                 alpha += [1]
#                 coil += [0]
#             else:
#                 beta += [0]
#                 alpha += [0]
#                 coil += [1]
#         plt.scatter(x,beta,label='beta',marker = 'o',color='blue')
#         plt.scatter(x,coil,label='coil',marker='x', color='green')
#         plt.scatter(x,alpha,label='alpha',color='red')
#         plt.title('Secondary structure prediction SOPM '+seqid)
#         plt.xlabel('Amino acid position')
#         plt.ylabel('Probability')
#         lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
#         ax=plt.gca()
#         fig = plt.gcf()
#         fig.set_size_inches
#         ax.set_xlim([0,len(prediction)])
#         ax.set_ylim([0.9,1.1])
#         plt.savefig('../Data/SOPM/'+seqid+'_actual.png',bbox_extra_artists=(lgd,),dpi=600,bbox_inches='tight')
#         plt.close()  
#            
#         SOPM_predictions+=[[beta,coil,alpha]]
#         struct=struct[:len(prediction)]
#         acc=(np.array(prediction)==np.array(struct)).sum()/len(prediction)
#         SOPM_accuracy+=acc
#         SOPM_accuracies[seqid]=acc
# 
# SOPM_accuracy=SOPM_accuracy/len(good_seqids)
# print("accuracy SOPM: "+str(SOPM_accuracy))
# GOR4_accuracies = {}
# GOR4_accuracy = 0
# 
# j=-1
#  
# GOR4_predictions = []
# for seqid in seqids:
#     j+=1
#     if seqid in good_seqids:
#         struct = structures[j]
#         prediction = []
#         for line in open('../Data/GORIV/'+seqid + '.gor4.mpsa.txt'):
#             if line[0]in ['H','C','E']  and len(prediction)<len(struct):
#                 pred = line[0]
#                 if pred=='H':
#                     prediction+=[1]
#                 elif pred=='E':
#                     prediction+=[-1]
#                 else:
#                     prediction+=[0]
#         print(seqid)
#         
#         x = range(len(prediction))
#         beta = []
#         alpha = []
#         coil = []
#         for amino in prediction:
#             if amino == -1:
#                 beta += [1]
#                 alpha += [0]
#                 coil += [0]
#             elif amino == 1:
#                 beta += [0]
#                 alpha += [1]
#                 coil += [0]
#             else:
#                 beta += [0]
#                 alpha += [0]
#                 coil += [1]
#         plt.scatter(x,beta,label='beta',marker = 'o',color='blue')
#         plt.scatter(x,coil,label='coil',marker='x', color='green')
#         plt.scatter(x,alpha,label='alpha',color='red')
#         plt.title('Secondary structure prediction GOR4 '+seqid)
#         plt.xlabel('Amino acid position')
#         plt.ylabel('Probability')
#         lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
#         ax=plt.gca()
#         fig = plt.gcf()
#         fig.set_size_inches
#         ax.set_xlim([0,len(prediction)])
#         ax.set_ylim([0.9,1.1])
#         plt.savefig('../Data/GORIV/'+seqid+'_actual.png',bbox_extra_artists=(lgd,),dpi=600,bbox_inches='tight')
#         plt.close()        
#         
#         GOR4_predictions+=[[beta,coil,alpha]]
#         struct=struct[:len(prediction)]
#         acc=(np.array(prediction)==np.array(struct)).sum()/len(prediction)
#         GOR4_accuracy+=acc
#         GOR4_accuracies[seqid]=acc
#  
# GOR4_accuracy=GOR4_accuracy/len(good_seqids)
# print("accuracy GORIV: "+str(GOR4_accuracy))
# 
# 
# data_file = '../Data/color_fasta.txt' + "_data.pik"
# with open(data_file, 'rb') as f:
#     data = pickle.load(f)
#     
# np.random.shuffle(data)
#     
# data_test, data_train = data[:len(data)/4,:], data[len(data)/4:,:]
#     
# x_train = data_train[:,1:]
# y_train = data_train[:,:1]
# x_test = data_test[:,1:]
# y_test = data_test[:,:1]
# 
# 
# logreg = linear_model.LogisticRegression(C=3.4)
# logreg.fit(x_train,y_train)
# data_file = '../Data/color_fasta.txt' + "_data.pik"
# with open(data_file, 'rb') as f:
#      data_full = pickle.load(f)
# 
# x_full = data_full[:,1:]
# y_full = data_full[:,:1]
# result_full = logreg.predict(x_full)
# 
# 
#==============================================================================


j=-1
i=0
k=-1

for seqid in seqids:
    j+=1
    if seqid in good_seqids:
        k+=1
        x = range(len(structures[j]))
        beta_real =[]
        alpha_real = []
        coil_real = []
        for amino in structures[j]:
            if amino == -1:
                beta_real += [1]
                alpha_real += [0]
                coil_real += [0]
            elif amino == 1:
                beta_real += [0]
                alpha_real += [1]
                coil_real += [0]
            else:
                beta_real += [0]
                alpha_real += [0]
                coil_real += [1]
        plt.scatter(x,beta_real,label='beta_real',marker = 'o',color='blue')
        plt.scatter(x,coil_real,label='coil_real',marker='x', color='green')
        plt.scatter(x,alpha_real,label='alpha_real',color='red')
        
        log_ = result_full[i:i+len(sequences[j])]
        beta_log =[]
        alpha_log = []
        coil_log = []
        for amino in log_:
            if amino == -1:
                beta_log += [0.9]
                alpha_log += [0]
                coil_log += [0]
            elif amino == 1:
                beta_log += [0]
                alpha_log += [0.9]
                coil_log += [0]
            else:
                beta_log += [0]
                alpha_log += [0]
                coil_log += [0.9]
        plt.scatter(x,beta_log,label='beta_log',marker = 'o',color='blue')
        plt.scatter(x,coil_log,label='coil_log',marker='x', color='green')
        plt.scatter(x,alpha_log,label='alpha_log',color='red')
        
        
        GOR4 = np.array(GOR4_predictions[k])*0.8
        x = range(len(GOR4[0]))
        plt.scatter(x,GOR4[0],label='beta_GOR4',marker = 'o',color='blue')
        plt.scatter(x,GOR4[1],label='coil_GOR4',marker='x', color='green')
        plt.scatter(x,GOR4[2],label='alpha_GOR4',color='red')

        SOPM = np.array(SOPM_predictions[k])*0.7
        x = range(len(SOPM[0]))
        plt.scatter(x,SOPM[0],label='beta_SOPM',marker = 'o',color='blue')
        plt.scatter(x,SOPM[1],label='coil_SOPM',marker='x', color='green')
        plt.scatter(x,SOPM[2],label='alpha_SOPM',color='red')
        
                
        s2D = np.array(s2D_predictions[k])*0.6
        x = range(len(s2D[0]))
        plt.scatter(x,s2D[0],label='beta_s2D',marker = 'o',color='blue')
        plt.scatter(x,s2D[1],label='coil_s2D',marker='x', color='green')
        plt.scatter(x,s2D[2],label='alpha_s2D',color='red')
        
        plt.title('Secondary structure prediction '+seqid)
        plt.xlabel('Amino acid position')
        plt.ylabel('Probability')
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
        ax=plt.gca()
        fig = plt.gcf()
        fig.set_size_inches
        ax.set_xlim([0,len(x)])
        ax.set_ylim([0.5,1.1])
        
        
        plt.savefig('../Data/'+seqid+'.png',bbox_extra_artists=(lgd,),dpi=600,bbox_inches='tight')
        plt.close() 
    i+=len(sequences[j])


