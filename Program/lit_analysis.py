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


struct_file = '../Data/color_fasta.txt' + "_struct_data.pik"
with open(struct_file, 'rb') as f:
    structures, sequences = pickle.load(f)   
     
new_file = '../Data/color_fasta.txt' + "_seqs_data.pik"
with open(new_file, 'rb') as f:
    sequences_2,seqids,names,descriptions = pickle.load(f)
    
good_seqids=['1bfp','3ekh','3ned','4k3g','3cfc',
        '1xkh','2wht','4w6b','4xvp','3dqh',
        '1bgp','4q7t','4qgw','5h88','4l1s',
        '5h89','3s0f','4q9w','3rwt','5hzo']

s2D_accuracies = {}
s2D_accuracy = 0

j=-1


for seqid in seqids:
    j+=1
    if seqid in good_seqids:
        struct = structures[j]
        i=0
        prediction = []
        for line in open('../Data/S2D/'+seqid + '_s2d_out.txt'):
            if line[0]!='>' and line[0]!='#':
                if i<11:
                    pred = line[40]
                elif i<101:
                    pred= line[41]
                else:
                    pred = line[42]
                if pred=='H':
                    prediction+=[1]
                elif pred=='E':
                    prediction+=[-1]
                else:
                    prediction+=[0]
            i+=1
        print(seqid)
       
        x = range(len(prediction))
        beta = []
        alpha = []
        coil = []
        for amino in prediction:
            if struct[amino] == -1:
                beta += [1]
                alpha += [0]
                coil += [0]
            elif struct[amino] == 1:
                beta += [0]
                alpha += [1]
                coil += [0]
            else:
                beta += [0]
                alpha += [0]
                coil += [1]
        plt.scatter(x,beta,label='beta',marker = 'o',color='blue')
        plt.scatter(x,coil,label='coil',marker='x', color='green')
        plt.scatter(x,alpha,label='alpha',color='red')
        plt.title('Secondary structure prediction GOR4 '+seqid)
        plt.xlabel('Amino acid position')
        plt.ylabel('Probability')
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
        ax=plt.gca()
        fig = plt.gcf()
        fig.set_size_inches
        ax.set_xlim([0,len(prediction)])
        ax.set_ylim([0.9,1.1])
        plt.savefig('../Data/S2D/'+seqid+'_actual.png',bbox_extra_artists=(lgd,),dpi=600,bbox_inches='tight')
        plt.close()               
       
        acc=(prediction==struct.flatten()).sum()/len(prediction)
        s2D_accuracy+=acc
        s2D_accuracies[seqid]=acc

s2D_accuracy=s2D_accuracy/len(seqids)


SOPM_accuracies = {}
SOPM_accuracy = 0

j=-1


for seqid in seqids:
    j+=1
    if seqid in good_seqids:
        struct = structures[j]
        prediction = []
        for line in open('../Data/SOPM/'+seqid + '.sopm.txt'):
            if line[0]in ['H','C','E'] and len(prediction)<len(struct):
                pred = line[0]
                if pred=='H':
                    prediction+=[1]
                elif pred=='E':
                    prediction+=[-1]
                else:
                    prediction+=[0]
        print(seqid)
        
        x = range(len(prediction))
        beta = []
        alpha = []
        coil = []
        for amino in prediction:
            if struct[amino] == -1:
                beta += [1]
                alpha += [0]
                coil += [0]
            elif struct[amino] == 1:
                beta += [0]
                alpha += [1]
                coil += [0]
            else:
                beta += [0]
                alpha += [0]
                coil += [1]
        plt.scatter(x,beta,label='beta',marker = 'o',color='blue')
        plt.scatter(x,coil,label='coil',marker='x', color='green')
        plt.scatter(x,alpha,label='alpha',color='red')
        plt.title('Secondary structure prediction GOR4 '+seqid)
        plt.xlabel('Amino acid position')
        plt.ylabel('Probability')
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
        ax=plt.gca()
        fig = plt.gcf()
        fig.set_size_inches
        ax.set_xlim([0,len(prediction)])
        ax.set_ylim([0.9,1.1])
        plt.savefig('../Data/SOPM/'+seqid+'_actual.png',bbox_extra_artists=(lgd,),dpi=600,bbox_inches='tight')
        plt.close()            
        
        struct=struct[:len(prediction)]
        acc=(np.array(prediction)==np.array(struct)).sum()/len(prediction)
        SOPM_accuracy+=acc
        SOPM_accuracies[seqid]=acc

SOPM_accuracy=SOPM_accuracy/len(good_seqids)
print("accuracy SOPM: "+str(SOPM_accuracy))

GOR4_accuracies = {}
GOR4_accuracy = 0

j=-1

for seqid in seqids:
    j+=1
    if seqid in good_seqids:
        struct = structures[j]
        prediction = []
        for line in open('../Data/GORIV/'+seqid + '.gor4.mpsa.txt'):
            if line[0]in ['H','C','E']  and len(prediction)<len(struct):
                pred = line[0]
                if pred=='H':
                    prediction+=[1]
                elif pred=='E':
                    prediction+=[-1]
                else:
                    prediction+=[0]
        print(seqid)
        
        x = range(len(prediction))
        beta = []
        alpha = []
        coil = []
        for amino in prediction:
            if struct[amino] == -1:
                beta += [1]
                alpha += [0]
                coil += [0]
            elif struct[amino] == 1:
                beta += [0]
                alpha += [1]
                coil += [0]
            else:
                beta += [0]
                alpha += [0]
                coil += [1]
        plt.scatter(x,beta,label='beta',marker = 'o',color='blue')
        plt.scatter(x,coil,label='coil',marker='x', color='green')
        plt.scatter(x,alpha,label='alpha',color='red')
        plt.title('Secondary structure prediction GOR4 '+seqid)
        plt.xlabel('Amino acid position')
        plt.ylabel('Probability')
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
        ax=plt.gca()
        fig = plt.gcf()
        fig.set_size_inches
        ax.set_xlim([0,len(prediction)])
        ax.set_ylim([0.9,1.1])
        plt.savefig('../Data/GORIV/'+seqid+'_actual.png',bbox_extra_artists=(lgd,),dpi=600,bbox_inches='tight')
        plt.close()        
        
        struct=struct[:len(prediction)]
        acc=(np.array(prediction)==np.array(struct)).sum()/len(prediction)
        GOR4_accuracy+=acc
        GOR4_accuracies[seqid]=acc

GOR4_accuracy=GOR4_accuracy/len(good_seqids)
print("accuracy GORIV: "+str(GOR4_accuracy))

