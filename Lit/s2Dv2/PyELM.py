'''
    Copyright (C) 2013 Piero Fariselli

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    PHMM lib Copyright (C) 2007  Piero Fariselli
    contacts: Piero Fariselli
              e-mail: piero@biocomp.unibo.it, piero.fariselli@unibo.it
              Dept. of Computer Science and Engineering
              University of Bologna
              via Mura Anteo Zamboni 7
              40127 Bologna
              Italy
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it
    under the condition of preserving this text
'''

import cPickle
import numpy as np
import sys,os
#from sklearn import linear_model

DEBUG=False

def logistic(x):
    ''' obj.actFun(x)
         => return 1/(1+exp(-x))
    '''
    return 0.5*(1.0+np.tanh(0.5*x))

# fragment the sequence with sliding windows
def buildWindSeq(X,w):
    '''
      buildWindSeq(X,w)
      adds zeros caps to the X ends of length w/2 and returns
      a sequence of the same length but w*len(X[0]) wide  
      X is the input vector, w is the sliding window size.
    '''
    Hw=w/2 # half width of the sliding window, recall that w should be even as the network predicts the value for the central position.
    vdim=len(X[0]) # will be the number of neurons per position in the sequence 
    D=np.vstack([np.zeros((Hw,vdim)),X,np.zeros((Hw,vdim))]) # extend the input vector with empty zeros corresponding to Hw positions before and after the beginning and the end of the sequence.
    D=np.hstack([D[i:len(D)-w+1+i] for i in xrange(w)])  # flatten the input vector. This yields an array of shape (sequence length , total number of inputs), where total number of inputs is window size * number of neurones per position. 
                                                        # in this way each position in the sequence has is own seqNetwork which will predict the property of that position
    return D

def tanh(x, out_scale=0.02,reverse=False):
    '''
    can be used to scale the coefficients in -1,1 and than expand them back (reverse=True)
    '''
    if not reverse : return np.tanh(out_scale*x)
    x=np.array(x)
    x[x>=1.]=0.9999999
    x[x<=-1.]=-0.9999999
    return (1./out_scale)*np.arctanh(x)
def identity(x, out_scale=1, reverse=False):
    return x

def smoothing_function(index_in_win,winSize,gaussian_exp=2):
    hW=int(winSize/2)
    den= hW**gaussian_exp / 1.386294 # it would be -X**2 /ln(0.25) so that at the end of the sliding window the weight will be 0.25
    return np.exp( -((index_in_win-hW)**gaussian_exp)/den)
def get_guided_weights_correction(win, numID, function=smoothing_function):
    '''
    for Seq networks.
    following E. Faraggi, B. Xue, Y. Zhou, Proteins 2009, 74, 847.
     it applies the guided learning correction to the random weights in the first layer of the ELM. the idea 
     is that residues closer to the central one will be more important in determining the property to be predicted
    j=1,...,J is the residue position in the sliding window (J=win)
    h=1,...,H is the number of hidden neurons0 (H=numH)
    '''
    return np.array( numID*[function(i,win) for i in xrange(win)] ).reshape(numID*win,1)

# reads input and output vectors. It expects a file with an input and the corresponding output.  
# you can have multiple input/output per file.
# it reads A file with a list of the above mentioned files.
class READER(object):
    def __init__(self,fname,indim,outdim):
        '''
        __init__(self,fname,indim,outdim)
         - fname the file name to open
         - indim input dimension of the input vector
         - output dimension of the output vector
        
          the file is expected the have input_k output_k 
          in the same k-th line  
        
        self.fp file pointer self.fp=open(fname)
        self.I current input vector
        self.O current terget vector
        '''
        self.nI=indim
        self.nO=outdim
        self.fp=open(fname)
        self.eof=False
        self.I=[]
        self.O=[]

    def close(self):
        ''' close() close the file '''
        self.fp.close()
        self.fp=None

    def next(self):
        ''' next() -> reads the next input/output'''
        line=self.fp.readline()
        if line =="": # EOF
            self.I=[]
            self.O=[]
            self.eof=True
            self.close()
        else:
            v=map(float, line.split())
            self.I=v[:self.nI]
            self.O=v[self.nI:self.nI+self.nO]
    def getIO(self):
        return self.I,self.O 
#-----------------------

# similar to the above. It assumes you have a lot of files. Each file has an input and the corresponding output but only ONE input
# this is handier if you change often input and output
class SeqREADER(READER):
    def __init__(self,fname,indim,outdim):
        '''
        __init__(self,fname,indim,outdim)
         - fname the file name to open
         - indim input dimension of the input vector
         - output dimension of the output vector
        
          the file is expected to contain a list of file
          and return one sequence at a time each reading step

        self.fp file pointer self.fp=open(fname). This contain a list of file names
        self.I current input vector
        self.O current terget vector
        '''
        self.nI=indim
        self.nO=outdim
        self.fname=fname
        self.fp=open(fname)
        self.eof=False
        self.I=[]
        self.O=[]

    def next(self):
        ''' next() -> reads the next input/output'''
        line=self.fp.readline()
        if line =="": # EOF
            self.I=[]
            self.O=[]
            self.eof=True
            self.close()
        else:
            v=line.split()
            name=v[0]  # name of SeqFile
            self.I=[]
            self.O=[]
            for fline in open(name):  # read lines in opened file
                v=map(float, fline.split())
                self.I.append(v[:self.nI])
                self.O.append(v[self.nI:self.nI+self.nO])
                if len(self.O[-1])!= self.nO :
                    sys.stderr.write('*Warn in SeqREADER() read %d output neurones but expecting %d, file is %s\n' % (len(self.O[-1]),self.nO,self.fname))
                    sys.stderr.flush()


#-----------------------


# special reader for Nto1
# it assumes to have a file with a list of files. Then each file is returned in the Nto1 format.
class Nto1READER(READER):
    def __init__(self,fname,indim,outdim):
        '''
        __init__(self,fname,indim,outdim)
         - fname the file name to open
         - indim input dimension of the input vector
         - output dimension of the output vector
        
          the file is expected the have a format as:
          filename value1 [value2 [value 3 .. ]]
          Where 
           - filename is the name of a file containing the sequence
           - value(s) is(are) real number(s) or a class information
      
        self.fp file pointer self.fp=open(fname). This contain a list of file names
        self.I current input vector
        self.O current terget vector
        '''
        self.nI=indim
        self.nO=outdim
        self.fp=open(fname)
        self.eof=False
        self.I=[]
        self.O=[]

    def next(self):
        ''' next() -> reads the next input/output'''
        line=self.fp.readline()
        if line =="": # EOF
            self.I=[]
            self.O=[]
            self.eof=True
            self.close()
        else:
            v=line.split()
            name=v[0]
            self.I=[]
            for fline in open(name):
                self.I.append(map(float, fline.split()))
            self.O=[float(x) for x in v[1:]]

class Nto1SeqREADER(READER):
    def __init__(self,fname,indim,outdim,NFcoeff=10,length_limit=100,preferred_slice_length=90,min_slices_overlap=3,exclude_Imag_0=True, avoid_triple_overlaps=True,out_ActFun=identity ):
        '''
        __init__(self,fname,indim,outdim)
         - fname the file name to open
         - indim input dimension of the input vector
         - output dimension of the output vector
        
          the file is expected to contain a list of files
          and return one sequence at a time each reading step
        NFcoeff is the number of complex Fourier coefficients to be used in the prediction.
        
        self.fp file pointer self.fp=open(fname). This contain a list of file names
        self.I current input vector
        self.O current terget vector
        '''
        self.nI=indim
        self.nO=outdim
        self.fname=fname
        self.fp=open(fname,'r')
        self.NFcoeff=NFcoeff # number of complex Fourier coefficient calculated from the profile that will be predicted 
        self.length_limit=length_limit # above this limit the sequence is split in overlapping subsets.
        self.min_slices_overlap=min_slices_overlap # mimum overlap between segments of the same sequence that got sliced because longer than length_limit
        self.preferred_slice_length=preferred_slice_length # length of the segments in which sequences longer than length_limit get sliced
        self.exclude_Imag_0=exclude_Imag_0
        self.out_ActFun=out_ActFun
        if self.preferred_slice_length==None : self.preferred_slice_length=self.length_limit
        self.avoid_triple_overlaps=avoid_triple_overlaps
        self.eof=False
        self.seqfile=None # can be used to debug
        self.stop_warnings=False
        self.profile=None
        self.fragment_storage=[]
        self.I=[]
        self.O=[]

    def next(self):
        ''' next() -> reads the next input/output'''
        if self.fragment_storage!=[] :
            self.I,self.O =self.fragment_storage[0]
            del self.fragment_storage[0]
            return
        self.stop_warnings=False
        line=self.fp.readline()
        if line =="": # EOF
            self.I=[]
            self.O=[]
            self.eof=True
            self.close()
        else:
            v=line.split()
            name=v[0]  # name of SeqFile
            self.seqfile=name
            self.profile=[] # profile of the Seq, it will be converted into Fourier coefficients.
            self.I=[]
            self.O=[]
            for fline in open(name,'r'):  # read lines in opened file, each line reperesents a position in the Seq and the output is at the end.
                v=map(float, fline.split())
                self.I.append(v[:self.nI])
                self.profile.append(v[self.nI:self.nI+self.nO])
                if not self.stop_warnings :
                    if len(self.profile[-1])!= self.nO :
                        sys.stderr.write('***WARNING*** in Nto1SeqREADER() read %d output neurones but expecting %d, file is %s [aborting further warnings for this file\n' % (len(self.profile[-1]),self.nO,self.fname))
                        sys.stderr.flush()
                        self.stop_warnings=True
                    elif len(v)!=self.nI + self.nO :
                        sys.stderr.write('***WARNING*** in Nto1SeqREADER() read %d possible input/outputs in line, but expecting %d input and %d outputs, file is %s [aborting further warnings for this file\n' % (len(v),self.nI,self.nO,self.fname))
                        sys.stderr.flush()
                        self.stop_warnings=True
            if self.length_limit!=None and len(self.profile)>self.length_limit :
                # split the sequence, save the first bit in self.I, self.O the other bits in fragment_storage so that in future next() call they are returned
                slices,_ = segmentize(self.I , self.preferred_slice_length, min_slices_overlap=self.min_slices_overlap, second_sliceable=self.profile, avoid_triple_overlaps=self.avoid_triple_overlaps)
                self.I,self.O = slices[0][0], self.out_ActFun(profile_to_Fourier_coeffs(slices[0][1],self.NFcoeff ) ) # save the first slice at this will be returned by getIO(), the rest is saved to be returned in the future calls of next()
                for sl in slices[1:] :
                    self.fragment_storage.append( ( sl[0], self.out_ActFun( profile_to_Fourier_coeffs(sl[1],self.NFcoeff, exclude_Imag_0=self.exclude_Imag_0) ) ) ) # append tuple that gets returned in future calls of next()
            else :
                #coeff=np.fft.rfft(np.array(profile),axis=0)[:self.NFcoeff:] # this are complex numbers, we need to flatten them to Real + Imag
                #self.O= flatten_FourierCoeffs(coeff, Ncoeff_if_reverse=False , exclude_Imag_0=self.exclude_Imag_0)
                self.O= self.out_ActFun( profile_to_Fourier_coeffs(self.profile,self.NFcoeff, exclude_Imag_0=self.exclude_Imag_0) )
                
            # memo from numpy man, coeff2=numpy.real(coeff)+1j*numpy.imag(coeff) ==> coeff and coeff2 are equal and this allows to use coefficients as real number in net 
            # recall also to give the length of the origianl length of the profile to np.fft.irfft( coeff, length)
            # coeff=np.fft.rfft(profile,axis=0)
            # profile=np.fft.irfft(coeff,len(profile),axis=0)




def segmentize(sliceable, Seglen, min_slices_overlap=0, second_sliceable=None, avoid_triple_overlaps=True) :
    '''
    avoid_triple_overlaps is necessary if you want to remerge the fragments
    get a sliceable of length slen and separate it in the minimum possible number of
    segments of length Seglen that includes all the positions with a minimum overlap of min_slices_overlap.
    Distributing the overlapping residues on every segment in the most possible equal way.
    if second_sliceable is not None it does the same for the second_sliceable but in this case it
      returns slices,istart, where slices is a list of tuples like (segment_from_sliceable , corresponding_segment_from_second_sliceable)
    returns slices, istart where istart is a list with the indices of the positions at which sliceable was sliced (first index of every segment)
    '''
    slen= len(sliceable)
    if slen< Seglen : return [ sliceable ]
    ns=(slen-min_slices_overlap)/(Seglen -min_slices_overlap)+1 # number of generated segments
    if avoid_triple_overlaps and Seglen*(ns-1)>slen : # avoid triple overlaps
        ns= slen/Seglen +1
        min_slices_overlap=0 
    
    ov= (Seglen -min_slices_overlap) - (slen-min_slices_overlap)%(Seglen -min_slices_overlap); # number of extra positions that will overlap
    nPov=ov/(ns-1) # Number of positions per overlap region (ns -1 overlap regions), plus some regions will have an extra position due to the remainders below
    r = ov%(ns-1)  # remainder, residues that will be distributed on the first overlapping regions, one for each till the r-region
    if slen==Seglen : ns,nPov,ov,r=1,0,0,0
    
    j=0
    inc=1
    if r==0 : inc=0 # in this case nPov will distribute them evenly
    istart=[0]
    if second_sliceable!=None :
        if len(second_sliceable)!=slen :
            sys.stderr.write("**ERROR** in segmentize, second_sliceable given but its length is different from the one of the main sliceable (%d %d)\n" % (len(second_sliceable),slen))
        slices =[ ( sliceable[0:Seglen], second_sliceable[0:Seglen] ) ]# save the first segment of length Seglen, it goes from position 0 to position Seglen-1 included
    else :    
        slices =[ sliceable[0:Seglen] ]# save the first segment of length Seglen, it goes from position 0 to position Seglen-1 included
    
    for i in range(1,ns) :
        istart+=[ istart[i-1]+Seglen-min_slices_overlap-nPov-inc ]
        j+=inc
        if j>= r : inc=0
        iend= istart[i]+Seglen-1
        istart[i]+Seglen-1
        if i==ns-1 and iend!=slen-1 :
            sys.stderr.write("**ERROR** in segmentize, last segments ends at %d but seqence has %d residues [ns=%d seglen=%d ov=%d nPov=%d r=%d min_slices_overlap=%d]\n" % (iend,slen,ns,Seglen,ov,nPov,r,min_slices_overlap))
        if second_sliceable!=None :  slices+= [ (sliceable[istart[i]:iend+1],second_sliceable[istart[i]:iend+1]) ]
        else : slices+= [ sliceable[istart[i]:iend+1] ]
        
    return slices, istart


def MergeSegments( slen, istart, slices, slice_seqtest=None, return_overlaps=None) :
    '''
    it does the opposite thing than the function segmentize
    slices must be numeric..
    it merges the elements in slices averaging the overlapping parts of consecutive segments
    slen will be the resulting length after the merging
    slice_seqtest can be given for debug and must be an object whose overapping parts are the same
    returns the merged profile
    if return_overlaps!=None returns also a list of the same length with 1 in the overlapping parts and 0 elesewhere
    '''
    ns=len(slices)
    Seglen=len(slices[0])
    seqspectrum=[]
    if return_overlaps!=None : ov_test=[]
    ev1=-1
    ev2=0
    od1=-1
    od2=0
    l1=0
    l2=0
    k=0

    for j in xrange(slen) :
        if j==ev2+Seglen : # end of segment even (recall at first j is zero)
            ev1=-1
            l1=0
        if j==od2+Seglen : #end of segment odd
            od1=-1
            l2=0 
        #print istart,k
        if k<ns and (j==istart[k] or k==0)  : #begin of segment
            k+=1
            if (k-1)%2==0 and (k-1)!=ns : # save initial position (along the sequence) in ev2 and segment id in ev1
                ev1=k-1 
                ev2=j 
            elif(k-1)<ns : # save initial position in od2 and segment id in od1 (necessary to account for the regions in which consecutive segments overlap)
                od1=k-1
                od2=j
        
        if ev1>=0 : #if we are looking at an even segment
            if(od1>=0) : #if it is an overlapping region
                seqspectrum += [ (slices[od1][l2] + slices[ev1][l1])/2. ] # average the two segments results
                if return_overlaps!=None : ov_test+=[1]

                if slice_seqtest!=None and slice_seqtest[od1][l2] != slice_seqtest[ev1][l1] :
                    sys.stderr.write("ERROR in MergeSegments overlapping regions do not match (input[%d][%d]%s!=input[%d][%d]%s) [od2=%d ev2=%d k=%d, istart[k]=%d slen=%d ns=%d j=%d istarts=" % (od1,l2,str(input[od1][l2]),ev1,l1,str(input[ev1][l1]),od2,ev2,k,istart[k],slen,ns,j))
                    for i in xrange(ns) :
                        sys.stderr.write("%d " % (istart[i]))
                    sys.stderr.write("]\n\n")
            
            else : # if it is not an overlapping region
            
                seqspectrum += [ slices[ev1][l1] ]
                if return_overlaps!=None : ov_test+=[0]
        
        elif od1>=0 : # if we are looking at an odd segment (note that this condition is checked after the even, so the check for overlaps is sufficient above).
            seqspectrum += [ slices[od1][l2] ]
            if return_overlaps!=None : ov_test+=[0]

        else :
            sys.stderr.write("ERROR problems in matching segments in MergeSegments [k=%d, istart[k]=%d slen=%d ns=%d j=%d istarts=" % (k,istart[k],slen,ns,j)) 
            for i in xrange(ns) : sys.stderr.write("%d " % (istart[i])) 
            sys.stderr.write("\n\n") 
            return 1
        
        if ev1>=0 : l1+=1
        if od1>=0 : l2+=1
    
    if return_overlaps!=None :
        return np.array(seqspectrum),ov_test
    return np.array(seqspectrum)

def Fourier_coeffs_to_profile(flattended_coeffs, NFcoeff, profile_length, exclude_Imag_0=True):
    complex_coeffs=flatten_FourierCoeffs(flattended_coeffs, Ncoeff_if_reverse=NFcoeff, exclude_Imag_0=exclude_Imag_0)
    return np.fft.irfft( complex_coeffs, profile_length, axis=0)
def profile_to_Fourier_coeffs(profile, NFcoeff, exclude_Imag_0=True):
    '''
    from a profile it returns a flattened array of its first NFcoeff Fourier coefficients.
    exclude_Imag_0=True exclude the imaginary part of the first coefficient, which is always zero for real profiles.
    '''
    coeff=np.fft.rfft(np.array(profile),axis=0)[:NFcoeff:] # this are complex numbers, we need to flatten them to Real + Imag
    return flatten_FourierCoeffs(coeff, Ncoeff_if_reverse=False, exclude_Imag_0=exclude_Imag_0)

def flatten_FourierCoeffs(coeffs, Ncoeff_if_reverse=False, exclude_Imag_0=True):
    '''
    converts ndarray of complex number to list of real number and vice versa (if Ncoeff_if_reverse instead of False is the number of complex Fourier coefficients per output profile)
     if you use it for reverse than profile=np.fft.irfft(coeff,len(profile),axis=0)
    exclude_Imag_0=True exclude the imaginary part of the first coefficient, which is always zero for real profiles.
    '''
    if False==Ncoeff_if_reverse and type(Ncoeff_if_reverse) is bool :
        return np.hstack( (np.real(coeffs.T),np.imag(coeffs.T)[:,int(exclude_Imag_0):]) ).flatten()
    else :
        n_profiles=len(coeffs)/(2*Ncoeff_if_reverse-int(exclude_Imag_0))
        #print 'n_profiles',n_profiles,coeffs.shape,Ncoeff_if_reverse,2*Ncoeff_if_reverse-int(exclude_Imag_0)
        all_real=coeffs.reshape(  (n_profiles, 2*Ncoeff_if_reverse-int(exclude_Imag_0) ) )
        if exclude_Imag_0 : return (all_real[:,:Ncoeff_if_reverse] + 1j* np.hstack( (np.zeros((n_profiles,1)),all_real[:,Ncoeff_if_reverse:] ))).T
        else :  return (all_real[:,:Ncoeff_if_reverse] + 1j* all_real[:,Ncoeff_if_reverse:]).T
#----------------------
# NEURAL NETWORK
#----------------------

# single layer feedforward network for extreme learning machine.
# 2 optioal parameters. The range in which you want to scale the initial random weights (the one that won't change within the Extreme Learning Machine).
#  if you have a large input layer (really large!) it is sometimes convenient to put a smaller wscale, 1 will do for most of the cases. 
#  IMPORTANT if you have e.g. 800 input neurons with a range in 0,1 each then something like 0.1 or even 0.01 could be more appropriate
# nbias, is a parameter that can give a shift to the output neurones. Since the hyperplane passes through the origin you can want it to pass elsewhere. This is the equivalent of
# the bias neuron in standard backpropagation networks. If you add it, it is authomatically added as a bias.
# In the context of Extreme Learning Machine sometimes it is not necessary but it help in some cases
# IMPORTANT in the typical situation with ELM you have to put a large number of numHnodes since the only free parameters are the elemnt of a matrix of size numHnodes X numOnodes
#  this is THE parameter you want to tune to your problem.
class SLFN(object):
    def __init__(self,numInodes,numHnodes,numOnodes,wscale=None,hbias=True):
        ''' 
             __init__(numInodes,numHnodes,numOnodes,wscale=None,hbias=True)
            self.numI=numInodes
            self.numH=numHnodes
            self.numO=numOnodes
            self.hbias=hbias
            self.numHB = self.numH + 1 # if hbias==True
            self.numHB = self.numH # if hbias!=True
            self.W      = array(self.numI x self.numH)
            self.Beta	= array(self.numHB x self.numO)
 
        '''
        self.actFun=logistic
        self.numI=numInodes
        self.numH=numHnodes
        self.numO=numOnodes
        self.hbias=hbias
        if hbias == True:
            self.numHB= self.numH + 1
        else:
            self.numHB= self.numH
        self.W=None
        self.Beta=None
        if type(wscale) == float: 
            self.inintWeights(wscale) 
   
    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the SLFN class
        '''
        self.reader=READER(fname,self.numI,self.numO)
    
    def inintWeights(self,wscale=1.0,wtype='symmetric'):
        ''' inintWeights(self,wscale=1.0,wtype='symmetric')
            sets W randomly in the range wrange
        '''
        self.initW(wscale,wtype)
        self.initBeta(wscale,wtype)
         
    def initBeta(self,wscale=1.0,wtype='symmetric'):
        ''' initBeta(self,wscale=1.0,wtype='symmetric')
            sets Beta weights randomly in the range wrange 
        '''
        if wtype=='symmetric':   
            self.Beta=np.array((np.random.random((self.numHB,self.numO))-0.5)*wscale)
        else:
            self.Beta=np.array((np.random.random((self.numHB,self.numO)))*wscale)
         
    def initW(self,wscale=1.0,wtype='symmetric'):
        ''' initW(self,wscale=1.0,wtype='symmetric')
            sets both W and Beta weights randomly in the range wrange 
        '''
        if wtype=='symmetric':   
            self.W=np.array((np.random.random((self.numI,self.numH))-0.5)*wscale) 
        else:
            self.W=np.array((np.random.random((self.numI,self.numH)))*wscale)

    def gotoHidden(self,X):
        ''' gotoHidden(self,X)
            retunrs the hidden from input 
        '''
        X=np.array(X)
        if self.hbias:
            N,dI=X.shape
            # assert dI == self.numI 
            Ones=np.ones((N,1))
            return np.hstack( (self.actFun(np.dot(X,self.W)), Ones) ) # if we consider a seq network the dot returns an array with shape (length of sequence, number of hidden neurones +1) where +1 is the bias neurone.
        else:
            return self.actFun(np.dot(X,self.W))

    def gotoOut(self,H):
        ''' gotoOut(self,H) 
            retuns the output from hidden 
        '''
        return np.dot(H,self.Beta)

    def predict(self,X):
        ''' self.predict(X)
            propagates the input X and returns the corresponding output 
        '''
        return np.array(self.gotoOut( self.gotoHidden(X) ))

    def propagate(self,X):
        ''' self.propagate(X)
            propagates the input X and returns the corresponding output 
        '''
        return self.gotoOut( self.gotoHidden(X) )

    def sqErr(self,X,T):
        ''' sqErr(self,X,T)
            computes the square errors between T and net.propagete(X) 
        '''
        P=self.propagate(X)
        return np.sum((np.array(P)-np.array(T))**2/self.numO)/len(P)
        
    def absErr(self,X,T):
        ''' absErr(self,X,T)
            computes mean abosulte error between T and net.propagete(X)
        '''
        P=self.propagate(X)
        return np.sum(abs(np.array(P)-np.array(T))/self.numO)/len(P)
        
 
    def save(self, fname):
        ''' obj.save(fname) 
            save obj on the file fname
        '''
        self.name=fname
        f=open(fname,"w")
        cPickle.dump(self,f)
        f.close()
#-------------------------------


#----------------------

# class with SLFN that assumes there is difference between the size of the INPUT of the problem and the input of the network.
# This ones already does the wrap for the various sequences. win is the Size of the sliding windows.
# the size of the neural network will be  self.numI=self.win*self.numID
# and each amino acid is represented by numID neurones.
# it ASSUMES that the size win is odd since there must be an amino acid in the middle for a sliding window
class SeqSLFN(SLFN):
    ''' extends SLFN assuming to read one-sequence at a time ''' 

    def __init__(self,numID,win,numH,numO,wscale=None,hbias=True,use_guided_learning=False):
        ''' 
             __init__(numI,win,numH,numO,wscale=None,hbias=True)
            self.numID=numID # dimension of the input vectors different from the network input
            self.win=win # window to be used to build the real sequence of input
            self.numI=self.win*self.numID
            self.numH=numH #-> total # of hidden neurons
            self.numO=numO # # of output nodes
            self.hbias=hbias # if bias on the output is set
            self.numHB = self.numH + 1 # if hbias==True
            self.numHB = self.numH     # if hbias!=True
            self.Beta   = array(self.numHB x self.numO)
 
        '''
        self.actFun=logistic
        self.numID=numID # number of neurons that represent a given position in the sequence, different from the network input size numI
        if win%2==0 :
            sys.stderr.write("***ERROR*** in SeqSLFN() the sliding window size must be odd (%d given) ADDING ONE!!\n", win)
            win+=1
        self.win=win # window to be used to build the real sequence of input
        self.numI=self.win*self.numID
        super (SeqSLFN,self).__init__(self.numI,numH,numO,wscale,hbias)
        if use_guided_learning : self.W *= get_guided_weights_correction(self.win, self.numID )

    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the SLFN class
        '''
        self.reader=SeqREADER(fname,self.numID,self.numO)

    def gotoHidden(self,X):
        ''' gotoHidden(self,X)
            retunrs the hidden from input 
            The input is expected to be a sequence, so that also the output
            must be a single sequence.
        '''
        X=buildWindSeq(np.array(X),self.win)
        return super (SeqSLFN,self).gotoHidden(X)


#----------------------
class SeqAvSLFN(SLFN):
    ''' 
        WARNINIG !!! Only partially written NEVER tested
        WARNINIG !!! Only partially written NEVER tested
        SeqAvSLFN(SeqSLFN) extend SeqSLFN addin average values of all the sequence positions 
        before and after each position
        A hidden layer of 3*H is used
    '''
    def __init__(self,numInodes,numOneThirdH,numOnodes,wscale=None,hbias=True):
        ''' 
             __init__(numInodes,numOneThirdH,numOnodes,wscale=None,hbias=True)
            self.numI=numInodes
            self.numH=3*numOneThirdH
            self.numO=numOnodes
            self.hbias=hbias
            self.numHB = self.numH + 1 # if hbias==True
            self.numHB = self.numH # if hbias!=True
            self.W      = array(self.numI x self.numH)
            self.Beta   = array(self.numHB x self.numO)
 
        '''
        super(SeqAvSLFN,self).__init__(numInodes,numOneThirdH,numOnodes,wscale,hbias)       
        self.numH=3*numOneThirdH
        if self.hbias:
            self.numHB = self.numH + 1
        else: 
            self.numHB = self.numH 

    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the SLFN class
        '''
        self.reader=SeqREADER(fname,self.numI,self.numO)

    def gotoHidden(self,X):
        ''' gotoHidden(self,X)
            retunrs the hidden from input 
        '''
        H=self.actFun(np.dot(X,self.W))
        N,dH=H.shape
        V=np.zeros((N+1,dH))
        V[0]=H[0]
        nH=[]
        for i in range(1,N):
            V[i]=V[i-1]+H[i]
        nH.append(np.hstack((V[-1],H[0],(V[N-1]-V[0])/(N-1))))
        for i in range(1,N-1):
            nH.append(np.hstack((V[i-1]/i,H[i],(V[N-1]-V[i])/(N-1-i))))
        nH.append(np.hstack((V[N-2]/(N-1),H[N-1],(V[N]))))
        if self.hbias:
            # assert dI == self.numI 
            Ones=np.ones((N,1))
            return np.hstack((np.array(H),Ones))
        else:
            return np.array(H)



#-------------------------------

class ForwardNet(SLFN):
    ''' ForwardNet extends SeqSLFN '''

    def __init__(self,numI,windFn,windCn,numHFn,numHCn,numO,wscale=None,hbias=True):
        ''' 
             __init__(numI,windFn,windCn,numHFn,numHCn,numO,wscale=None,hbias=True)
            self.numI=numI # dimension fo the input vectors
            self.windFn=windFn # window to be used to build Forward input
            self.numIFn=numI*windFn +self.numHFn#-> # of input nodes of the recurrent net
            self.numHFn=numHFn #-> # of hidden nodes of the recurrent Forward part
            self.numICn=numI*windCn #-> # of input nodes of the FF non-recurrent part 
            self.numHCn=numHCn #-> # of hidden nodes  of the FF non-recurrent part
            self.numH=numHFn+numHCn #-> total # of hidden neurons
            self.numO=numOnodes # # of output nodes
            self.hbias=hbias # if bias on the output is set
            self.numHB = self.numH + 1 # if hbias==True
            self.numHB = self.numH     # if hbias!=True
            self.WFn =   # array( (self.numIF+self.numHF)x self.numHF) # w I-to-H Forward part
            self.Cnet =   # SLFN(numIC,numHC,numO=None,wscale=None,hbias=False)
            self.Beta   = array(self.numHB x self.numO)
 
        '''
        self.actFun=logistic
        self.numI=numI # dimension fo the input vectors
        self.windFn=windFn # window to be used to build Forward input
        self.numHFn=numHFn #-> # of hidden nodes of the recurrent Forward part
        if self.numHFn == 0 or self.windFn == 0:
            self.hasToWheelF=False
            self.numHFn=self.windFn=0
        else:
            self.hasToWheelF=True   
        self.numIFn=numI*self.windFn+self.numHFn #-> # of input nodes of the recurrent I_j part
        self.windCn=windCn # window to be used to build Forward input
        self.numICn=numI*windCn #-> # of input nodes of the FF non-recurrent part 
        self.numHCn=numHCn #-> # of hidden nodes  of the FF non-recurrent part
        self.numH=self.numHFn+self.numHCn #-> total # of hidden neurons
        self.numO=numO # # of output nodes
        self.hbias=hbias # if bias on the output is set
        if hbias == True:
            self.numHB= self.numH + 1
        else:
            self.numHB= self.numH
        self.WFn=None
        self.Beta=None
        self.Cnet=SLFN(self.numICn,self.numHCn,None,None,False)
        if type(wscale) == float:
            self.inintWeights(wscale)
        
    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the SLFN class
        '''
        self.reader=SeqREADER(fname,self.numI,self.numO)

    def initWF(self,wscale=1.0,wtype='symmetric'):
        ''' initWF(self,wscale=1.0,wtype='symmetric')
        '''
        if not self.hasToWheelF:
            return
        if wtype=='symmetric':
            self.WFn=np.array((np.random.random((self.numIFn,self.numHFn))-0.5)*wscale)
        else:
            self.WFn=np.array((np.random.random((self.numIFn,self.numHFn)))*wscale)
    

    def initW(self,wscale=1.0,wtype='symmetric'):
        ''' initW(self,wscale=1.0,wtype='symmetric')
            sets both W and Beta weights randomly in the range wrange 
        '''
        self.initWF(wscale,wtype)
        self.Cnet.initW(wscale)

    def wheelF(self,X):
        '''
            wheelF(self,X) 
            assumes that X is a sequence of vectors of dimension self.numI
            returns a sequence of the same length of vectors of dimesion self.numHF
        '''
        if not self.hasToWheelF:
            return None
        xL=len(X)
        F=np.zeros((xL,self.numHFn))
        start=np.zeros(self.numHFn)
        SEQ=buildWindSeq(X,self.windFn)
        F[0]=self.actFun(np.dot(np.hstack((start,SEQ[0])),self.WFn))
        for i in range(1,xL):
            F[i]=self.actFun(np.dot(np.hstack((F[i-1],SEQ[i])),self.WFn))
        SEQ=None
        return F
        
    def gotoHidden(self,X):
        ''' gotoHidden(self,X)
            retunrs the hidden from input 
            The input is expected to be a sequence, so that also the output
            must be a single sequence.
        '''
        X=np.array(X)
        F=self.wheelF(X)
        C=self.Cnet.gotoHidden(buildWindSeq(X,self.windCn))
        N,dI=X.shape
        if not self.hasToWheelF and self.hbias:
            Ones=np.ones((N,1))
            return np.hstack((C,Ones))
        elif not self.hasToWheelF:
            return C
          
        if self.hbias:
            # assert dI == self.numI 
            Ones=np.ones((N,1))
            return np.hstack((F,C,Ones))
        else:
            return np.hstack((F,C))

#-------------------------------


class BackwardNet(SLFN):
    ''' BackwardNet extends SeqSLFN '''

    def __init__(self,numI,windBn,windCn,numHBn,numHCn,numO,wscale=None,hbias=True):
        ''' 
             __init__(numI,winBn,windCn,numHBn,numHCn,numO,wscale=None,hbias=True)
            self.numI=numI # dimension fo the input vectors
            self.windBn=windBn # window to be used to build Backward input
            self.numIBn=numI*windBn +self.numHFn#-> # of input nodes of the recurrent net
            self.numHBn=numHBn #-> # of hidden nodes of the recurrent Backward part
            self.numICn=numI*windCn #-> # of input nodes of the FF non-recurrent part 
            self.numHCn=numHCn #-> # of hidden nodes  of the FF non-recurrent part
            self.numH=numHFn+numHCn #-> total # of hidden neurons
            self.numO=numOnodes # # of output nodes
            self.hbias=hbias # if bias on the output is set
            self.numHB = self.numH + 1 # if hbias==True
            self.numHB = self.numH     # if hbias!=True
            self.WBn =   # array( (self.numIBn+self.numHBn)x self.numHBn) # w I-to-H Backward part
            self.Cnet =   # SLFN(numICn,numHCn,numO=None,wscale=None,hbias=False)
            self.Beta   = array(self.numHB x self.numO)
 
        '''
        self.actFun=logistic
        self.numI=numI # dimension fo the input vectors
        self.windBn=windBn # window to be used to build Backward input
        self.numHBn=numHBn #-> # of hidden nodes of the recurrent Backward part
        if self.windBn==0 or self.numHBn==0:
            self.hasToWheelB=False
            self.windBn=self.numHBn=0
        else:
            self.hasToWheelB=True 
        self.numIBn=numI*self.windBn+self.numHBn #-> # of input nodes of the recurrent I_j part
        self.windCn=windCn # window to be used to build Backward input
        self.numICn=numI*windCn #-> # of input nodes of the FF non-recurrent part 
        self.numHCn=numHCn #-> # of hidden nodes  of the FF non-recurrent part
        self.numH=self.numHBn+self.numHCn #-> total # of hidden neurons
        self.numO=numO # # of output nodes
        self.hbias=hbias # if bias on the output is set
        if hbias == True:
            self.numHB= self.numH + 1
        else:
            self.numHB= self.numH
        self.WBn=None
        self.Beta=None
        self.Cnet=SLFN(self.numICn,self.numHCn,None,None,False)
        if type(wscale) == float:
            self.inintWeights(wscale)
        
    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the SLFN class
        '''
        self.reader=SeqREADER(fname,self.numI,self.numO)

    def initWB(self,wscale=1.0,wtype='symmetric'):
        ''' initWB(self,wscale=1.0,wtype='symmetric')
        '''
        if not self.hasToWheelB:
            return
        if wtype=='symmetric':
            self.WBn=np.array((np.random.random((self.numIBn,self.numHBn))-0.5)*wscale)
        else:
            self.WBn=np.array((np.random.random((self.numIBn,self.numHBn)))*wscale)
        

    def initW(self,wscale=1.0,wtype='symmetric'):
        ''' initW(self,wscale=1.0,wtype='symmetric')
        '''
        self.initWB(wscale,wtype)
        self.Cnet.initW(wscale)

    def wheelB(self,X):
        '''
            wheelB(self,X) 
            assumes that X is a sequence of vectors of dimension self.numI
            returns a sequence of the same length of vectors of dimesion self.numHBn
        '''
        if not self.hasToWheelB:
            return None
        xL=len(X)
        B=np.zeros((xL,self.numHBn))
        start=np.zeros(self.numHBn)
        SEQ=buildWindSeq(X,self.windBn)
        B[-1]=self.actFun(np.dot(np.hstack((start,SEQ[-1])),self.WBn))
        for i in range(xL-1,0,-1):
            B[i-1]=self.actFun(np.dot(np.hstack((B[i],SEQ[i-1])),self.WBn))
        SEQ=None
        return B
        
    def gotoHidden(self,X):
        ''' gotoHidden(self,X)
            retunrs the hidden from input 
            The input is expected to be a sequence, so that also the output
            must be a single sequence.
        '''
        X=np.array(X)
        C=self.Cnet.gotoHidden(buildWindSeq(X,self.windCn))
        B=self.wheelB(X)
        N,dI=X.shape
        if not self.hasToWheelB and self.hbias:
            Ones=np.ones((N,1))
            return np.hstack((C,Ones))
        elif not self.hasToWheelB:
            return C
        if self.hbias:
            # assert dI == self.numI 
            Ones=np.ones((N,1))
            return np.hstack((C,B,Ones))
        else:
            return np.hstack((C,B))


#-------------------------------

class BRNet(ForwardNet,BackwardNet):
    ''' BRNet extends ForwardNet,BackwardNet '''

    def __init__(self,numI,windFn,windCn,windBn,numHFn,numHCn,numHBn,numO,wscale=None,hbias=True):
        ''' 
             __init__(numI,windFn,windCn,numHFn,numHCn,numO,wscale=None,hbias=True)
            self.numI=numI # dimension fo the input vectors
            self.windFn=windFn # window to be used to build Forward input
            self.numIFn=numI*windFn +self.numHFn#-> # of input nodes of the recurrent net
            self.numHFn=numHFn #-> # of hidden nodes of the recurrent Forward part
            self.numICn=numI*windCn #-> # of input nodes of the FF non-recurrent part 
            self.numHCn=numHCn #-> # of hidden nodes  of the FF non-recurrent part
            self.windBn=windBn # window to be used to build Backward input
            self.numIBn=numI*windBn +self.numHFn#-> # of input nodes of the recurrent net
            self.numHBn=numHBn #-> # of hidden nodes of the recurrent Backward part
            self.numH=numHFn+numHCn #-> total # of hidden neurons
            self.numO=numOnodes # # of output nodes
            self.hbias=hbias # if bias on the output is set
            self.numHB = self.numH + 1 # if hbias==True
            self.numHB = self.numH     # if hbias!=True
            self.WFn =   # array( (self.numIF+self.numHF)x self.numHF) # w I-to-H Forward part
            self.Cnet =   # SLFN(numIC,numHC,numO=None,wscale=None,hbias=False)
            self.Beta   = array(self.numHB x self.numO)
 
        '''
        super(BRNet,self).__init__(numI,windFn,windCn,numHFn,numHCn,numO,wscale=None,hbias=False)
        self.windBn=windBn # window to be used to build Backward input
        self.numHBn=numHBn #-> # of hidden nodes of the recurrent Backward part
        if self.windBn==0 or self.numHBn==0:
            self.hasToWheelB=False
            self.windBn=self.numHBn=0
        else:
            self.hasToWheelB=True
        self.numIBn=numI*self.windBn+self.numHBn #-> # of input nodes of the recurrent I_j part
        self.numH=self.numHFn+self.numHCn+self.numHBn #-> total # of hidden neurons
        self.hbias=hbias # if bias on the output is set
        if hbias == True:
            self.numHB= self.numH + 1
        else:
            self.numHB= self.numH
        self.WCn=None
        self.WBn=None
        self.Beta=None
        self.Cnet=SLFN(self.numICn,self.numHCn,None,None,False)
        if type(wscale) == float:
            self.inintWeights(wscale)

    def initW(self,wscale=1.0,wtype='symmetric'):
        ''' initW(self,wscale=1.0,wtype='symmetric')
        '''
        self.initWF(wscale,wtype)
        self.initWB(wscale,wtype)
        self.Cnet.initW(wscale)

    def gotoHidden(self,X):
        ''' gotoHidden(self,X)
            retunrs the hidden from input 
            The input is expected to be a sequence, so that also the output
            must be a single sequence.
        '''
        X=np.array(X)
        F=self.wheelF(X)
        C=self.Cnet.gotoHidden(buildWindSeq(X,self.windCn))
        B=self.wheelB(X)
        if self.hbias:
            N,dI=X.shape
            # assert dI == self.numI 
            Ones=np.ones((N,1))
            if self.hasToWheelF and self.hasToWheelB:
                return np.hstack((F,C,B,Ones))
            elif self.hasToWheelF:
                return np.hstack((F,C,Ones))
            elif self.hasToWheelB:
                return np.hstack((C,B,Ones))
            else:
                return np.hstack((C,Ones))
        else:
            if self.hasToWheelF and self.hasToWheelB:
                return np.hstack((F,C,B))
            elif self.hasToWheelF:
                return np.hstack((F,C))
            elif self.hasToWheelB:
                return np.hstack((C,B))
            else:
                return C





def flinear(self,x):
        return self.hscale*(self.actFun(np.dot(x,self.W)).sum(axis=0))

def flog(self,x):
        return self.hscale*np.log(self.actFun(np.dot(x,self.W)).sum(axis=0)+1.0)

def fsqrt(self,x):
        return self.hscale*np.sqrt(self.actFun(np.dot(x,self.W)).sum(axis=0))

def faverage(self,x):
        return self.actFun(np.dot(x,self.W)).mean(axis=0)

#-------------------------------
# takes input of SLFN, similar to seq_SLFN but without the assumption of having an odd win (here window)
"""
ora quando inizializzi una Nto1 o WNto1 puo usare come prima:
PyELM.WNto1(windows,dimInput,numHw,numO,wscale=1.0,hbias=True)

Oppure specificare che vuoi una funzione di trasferimento diverso con:

PyELM.WNto1(windows,dimInput,numHw,numO,wscale=1.0,hbias=True,hfun="average")
PyELM.WNto1(windows,dimInput,numHw,numO,wscale=1.0,hbias=True,hfun="linear",hscale=1.0)
PyELM.WNto1(windows,dimInput,numHw,numO,wscale=1.0,hbias=True,hfun="sqrt",hscale=1.0)
PyELM.WNto1(windows,dimInput,numHw,numO,wscale=1.0,hbias=True,hfun="log",hscale=1.0)
Dove hai:
      if hfun=="linear":
              hfun(x)= hscale*x
      elif hfun=="log":
              hfun(x)= hscale*log(x+1)
      elif hfun=="sqrt":
              hfun(x)= hscale*sqrt(x)
      elsif hfun=="average": # default "average"
              hfun(x)= x.mean(axis=0)

Cosi' puoi provare 4 modi  diversi
"""
class Nto1(SLFN):
    def __init__(self,window,dimInput,numHnodes,numOnodes,wscale=None,hbias=True,hfun="average",hscale=1.0):
        '''
            __init__(self,window,dimInput,numHnodes,numOnodes,wscale=None,hbias=True,hfun="average",hscale=1.0)
            if hfun=="linear":
               hfun(x)= hscale*x
            elif hfun=="log":
               hfun(x)= hscale*log(x+1)
            elif hfun=="sqrt":
               hfun(x)= hscale*sqrt(x)
            else: # default "average"
               hfun(x)= x.mean(axis=0)
 
        '''
        self.window=window
        self.hscale=hscale
        if hfun=="linear":
            self.Hfun= flinear
        elif hfun=="log":
            self.Hfun= flog
        elif hfun=="sqrt":
            self.Hfun= fsqrt
        elif hfun=="average":
            self.Hfun= faverage
        else:
            raise Exception("ERROR: the hfun = "+str(hfun)+" is undefined\n")
        self.dimInput=dimInput # dimension fo the read input vector
        super(Nto1,self).__init__(dimInput*window,numHnodes,numOnodes,wscale,hbias)
    
    
    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the SLFN class
        '''
        self.reader=Nto1READER(fname,self.numI,self.numO)

    def gotoHidden(self,A):
        ''' '''
        # we assume A is 3D: 
        #  - first index equal the number of targets, 
        #  - second index the sequence length
        #  - third index the input dimensions self.dimInput
        w=self.window
        N=len(A)
        H=[]
        for X in A:
            lX=len(X)
            if lX<w :
                raise Exception("sequence shorter than window! (%d %d)\n" % (lX,w))
                #sys.exit(lX)
            X=np.array(X)
            #print 'self.W.shape,X.shape,N',self.W.shape,X.shape,N
            D=np.hstack([X[i:lX-w+1+i] for i in range(w) ])
            #print 'D.shape',D.shape
            H.append(self.Hfun(self,D))
        if self.hbias==True:
            Ones=np.ones((N,1))
            return  np.hstack((H,Ones))
        else:
            return np.array(H)

    def deconvolve(self,seq):
        ''' deconvolve(seq)
            assumes seq is a single sequence that is going to be deconvolved
            i.e. for each position a prediction is made as in the case of SLFN 
        ''' 
        n=SLFN(self.numI,self.numH,self.numO,wscale=None,hbias=True)
        n.Beta=self.Beta
        n.W=self.W
        X=np.array(seq)
        lX=len(X)
        w=self.window
        X=np.hstack([X[i:lX-w+1+i] for i in range(w) ])
        return n.predict(X)

#--------------------------



class Nto1Seq(Nto1): # THERE SHOULD BE A MINIMUM LENGTH MAYBE
    def __init__(self,window,dimInput,numHnodes,numOnodes,wscale=None,hbias=True,hfun="average",hscale=1.0, NFcoeff=10,length_limit=100,preferred_slice_length=90,min_slices_overlap=3,exclude_Imag_0=True, maximize_training_size=True, out_ActFun=identity):
        '''
            __init__(self,windows,dimInput,numHnodes,numOnodes,wscale=None,hbias=True,hfun="average",hscale=1.0)
             
        '''
        self.NFcoeff=NFcoeff # number of complex Fourier coefficient calculated from the profile that will be predicted 
        self.length_limit=length_limit # above this limit the sequence is split in overlapping subsets.
        self.min_slices_overlap=min_slices_overlap # mimum overlap between segments of the same sequence that got sliced because longer than length_limit
        self.preferred_slice_length=preferred_slice_length # length of the segments in which sequences longer than length_limit get sliced
        self.exclude_Imag_0=exclude_Imag_0
        self.nProfiles=numOnodes # the number of different outputs
        self.numID=dimInput # so that this is compatible with a seq network
        self.out_ActFun=out_ActFun # scaling for the Fourier coefficients
        self.maximize_training_size=maximize_training_size
            
        numOnodes= self.nProfiles*( 2*NFcoeff-int(self.exclude_Imag_0))
        super(Nto1Seq,self).__init__(window,dimInput,numHnodes,numOnodes,wscale,hbias,hfun,hscale)
        

    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the Nto1Seq class
        '''
        if self.maximize_training_size :
            min_slices_overlap=self.preferred_slice_length-1
        else :
            min_slices_overlap=self.min_slices_overlap
        self.reader=Nto1SeqREADER(fname,self.dimInput,self.nProfiles, NFcoeff=self.NFcoeff, length_limit=self.length_limit, preferred_slice_length=self.preferred_slice_length, min_slices_overlap=min_slices_overlap,exclude_Imag_0=self.exclude_Imag_0,avoid_triple_overlaps=not self.maximize_training_size, out_ActFun=self.out_ActFun)
    
    def predict(self,X):
        ''' self.predict(X)
            propagates the input X and returns the corresponding output 
        '''
        if len(X) > self.length_limit :
            slices, istart= segmentize( X, self.preferred_slice_length, min_slices_overlap=self.min_slices_overlap, second_sliceable=None,avoid_triple_overlaps=True)
            out_slices=[]
            for sl in slices :
                flattened_coefficients= self.out_ActFun( np.array(self.gotoOut( self.gotoHidden(np.array([sl])) )) , reverse=True )
                out_slices.append( Fourier_coeffs_to_profile(flattened_coefficients.T, self.NFcoeff, len(sl), exclude_Imag_0=self.exclude_Imag_0) )
            return MergeSegments(len(X), istart, out_slices, slice_seqtest=None, return_overlaps=None)
        else :
            flattened_coefficients= self.out_ActFun( np.array(self.gotoOut( self.gotoHidden(np.array([X])) )) , reverse=True )
            return Fourier_coeffs_to_profile(flattened_coefficients.T, self.NFcoeff, len(X), exclude_Imag_0=self.exclude_Imag_0)
    def propagate(self,X):
        '''
        like predict but returns the Fourier coefficients, not the profile
         return out_slices,istart,len(X)  where out_slices contains the coefficient for each slice. If the sequence is not sliced istart is None
        '''
        if len(X) > self.length_limit :
            slices, istart= segmentize( X, self.preferred_slice_length, min_slices_overlap=self.min_slices_overlap, second_sliceable=None,avoid_triple_overlaps=True)
            out_slices=[]
            for sl in slices :
                flattened_coefficients= self.out_ActFun( np.array(self.gotoOut( self.gotoHidden(np.array([sl])) )) ,reverse=True)
                # rather than converting them to complex number we returned the real coefficients, but reshaped according to which profile they belong to
                n_profiles=len(flattened_coefficients.T)/(2*self.NFcoeff-int(self.exclude_Imag_0))
                all_real=(flattened_coefficients.T).reshape(  (n_profiles, 2*self.NFcoeff-int(self.exclude_Imag_0) ) )
                out_slices.append( all_real )
            return out_slices,istart,len(X)
        else :
            flattened_coefficients= self.out_ActFun( np.array(self.gotoOut( self.gotoHidden(np.array([X])) )) ,reverse=True)
            # rather than converting them to complex number we returned the real coefficients, but reshaped according to which profile they belong to
            n_profiles=len(flattened_coefficients.T)/(2*self.NFcoeff-int(self.exclude_Imag_0))
            all_real=(flattened_coefficients.T).reshape(  (n_profiles, 2*self.NFcoeff-int(self.exclude_Imag_0) ) )
            #complex_coeffs=flatten_FourierCoeffs(flattened_coefficients.T, Ncoeff_if_reverse=self.NFcoeff, exclude_Imag_0=self.exclude_Imag_0)
            return [all_real],None, len(X)

class Nto1R(Nto1):
    def __init__(self,window,dimInput,numHnodes,numOnodes,Rcut,wscale=None,hbias=True,hfun="average",hscale=1.0):
        '''
            Nto1CutR as Nto1 but takes only the k-th rightmost elements
            if k >=len(seq) it takes the whole sequence
            __init__(self,window,dimInput,numHnodes,numOnodes,Rcut,wscale=None,hbias=True,hfun="average",hscale=1.0)
 
        '''
        self.rcut=Rcut # dimension fo the read input vector
        super(Nto1R,self).__init__(window,dimInput,numHnodes,numOnodes,wscale,hbias,hfun,hscale)

    def gotoHidden(self,A):
        ''' '''
        # we assume A is 3D: 
        #  - first index equal the number of targets, 
        #  - second index the sequence length
        #  - third index the input dimensions self.dimInput
        w=self.window
        N=len(A)
        H=[]
        for Xs in A:
            lX=len(Xs)
            if self.rcut + w < lX:
                X=Xs[-self.rcut:]
                lX=self.rcut
            else:
                X=Xs
            if lX<w :
                raise Exception("sequence shorter than window! (%d %d)\n" % (lX,w))
                #sys.exit(lX)
            X=np.array(X)
            D=np.hstack([X[i:lX-w+1+i] for i in range(w) ])
            H.append(self.Hfun(self,D))
        if self.hbias==True:
            Ones=np.ones((N,1))
            return  np.hstack((H,Ones))
        else:
            return np.array(H)




class Nto1L(Nto1):
    def __init__(self,window,dimInput,numHnodes,numOnodes,Lcut,wscale=None,hbias=True,hfun="average",hscale=1.0):
        '''
            Nto1CutR as Nto1 but takes only the k-th leftmost elements
            if k >=len(seq) it takes the whole sequence
            __init__(self,window,dimInput,numHnodes,numOnodes,Rcut,wscale=None,hbias=True,hfun="average",hscale=1.0)
 
        '''
        self.lcut=Lcut # dimension fo the read input vector
        super(Nto1L,self).__init__(window,dimInput,numHnodes,numOnodes,wscale,hbias,hfun,hscale)

    def gotoHidden(self,A):
        ''' '''
        # we assume A is 3D: 
        #  - first index equal the number of targets, 
        #  - second index the sequence length
        #  - third index the input dimensions self.dimInput
        w=self.window
        N=len(A)
        H=[]
        for Xs in A:
            lX=len(Xs)
            if self.lcut + w < lX:
                X=Xs[:self.lcut]
                lX=self.lcut
            else:
                X=Xs
            if lX<w :
                raise Exception("sequence shorter than window! (%d %d)\n" % (lX,w))
                #sys.exit(lX)
            X=np.array(X)
            D=np.hstack([X[i:lX-w+1+i] for i in range(w) ])
            H.append(self.Hfun(self,D))
        if self.hbias==True:
            Ones=np.ones((N,1))
            return  np.hstack((H,Ones))
        else:
            return np.array(H)

           
#--------------------------------------

# instead of having a single window it has a series of window.
# it assumes that you can slide on the sequence with a window of 1, one of 5 one of 10 (for example!)
# you have to give a list of windows, and consequently a list of nodes numHnodes.
# This are not SEPARATED window. EACH window fills the same HIDDEN layer before the final propagation to the output.
# This Hidden layer will be given by the sum(numHnodes) (which is a list). In fact each window fills a different region of 
# the hidden layer, then the final propagation is done all at once!
class WNto1(SLFN):
    def __init__(self,windows,dimInput,numHnodes,numOnodes,wscale=None,hbias=True,hfun="average",hscale=1.0):
        '''
            __init__(self,windows,dimInput,numHnodes,numOnodes,wscale=None,hbias=True,hfun="average",hscale=1.0)
 
        '''
        self.windows=windows
        self.Nwindows=len(windows)
        self.dimInput=dimInput # dimension fo the read input vector
        super(WNto1,self).__init__(None,sum(numHnodes),numOnodes,wscale=False,hbias=hbias)
        self.initBeta(wscale)
        self.nets=[None]*self.Nwindows
        for i in range(self.Nwindows):
            self.nets[i]=Nto1(windows[i],dimInput,numHnodes[i],numOnodes,wscale=None,hbias=False,hfun=hfun,hscale=hscale)
            self.nets[i].initW(wscale)

    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the SLFN class
        '''
        self.reader=Nto1READER(fname,self.numI,self.numO)

    def gotoHidden(self,A):
        ''' '''
        # we assume A is 3D: 
        #  - first index equal the number of targets, 
        #  - second index the sequence length
        #  - third index the input dimensions self.dimInput
        H=np.hstack((self.nets[i].gotoHidden(A) for i in range(self.Nwindows)))
        N=len(A)
        if self.hbias==True:
            Ones=np.ones((N,1))
            return  np.hstack((H,Ones))
        else:
            return np.array(H)


#--------------------------------------

class WNto1SLR(SLFN):
    def __init__(self,windows,dimInput,numHnodes,numOnodes,Lcut=0,Rcut=0,wscale=None,hbias=True,hfun="average",hscale=1.0):
        '''
           SLR= Sequence, LeftSequence, RightSequence
           This is a N-to-1 network that uses also the same winodows and hidden for the letfmost (Lcut)
           or the rightmost (Rcut) sequence part     
            __init__(self,windows,dimInput,numHnodes,numOnodes,Lcut=0,Rcut=0,wscale=None,hbias=True,hfun="average",hscale=1.0)
 
        '''
        self.lcut=Lcut
        self.rcut=Rcut
        self.windows=windows
        self.numWindows=len(windows)
        self.dimInput=dimInput # dimension fo the read input vector
        count=1
        self.SnumNets=self.numWindows
        self.LnumNets=self.RnumNets=0
        if self.lcut > 0:
            count+=1
            self.LnumNets=self.numWindows 
        if self.rcut > 0:
            count+=1
            self.RnumNets=self.numWindows
        totHnodes=sum(numHnodes)*count
        self.numNets=self.numWindows*count 
        super(WNto1SLR,self).__init__(None,totHnodes,numOnodes,wscale=False,hbias=hbias)
        self.initBeta(wscale)
        self.nets=[None]*self.numNets
        iNet=0
        for i in range(self.numWindows):
            self.nets[iNet]=Nto1(windows[i],dimInput,numHnodes[i],numOnodes,wscale=None,hbias=False,hfun=hfun,hscale=hscale)
            self.nets[iNet].initW(wscale)
            iNet+=1
        if self.lcut > 0:
            for i in range(self.numWindows):
                self.nets[iNet]=Nto1L(windows[i],dimInput,numHnodes[i],numOnodes,self.lcut,wscale=None,hbias=False,hfun=hfun,hscale=hscale)
                self.nets[iNet].initW(wscale)
                iNet+=1
        if self.rcut > 0:
            for i in range(self.numWindows):
                self.nets[iNet]=Nto1R(windows[i],dimInput,numHnodes[i],numOnodes,self.rcut,wscale=None,hbias=False,hfun=hfun,hscale=hscale)
                self.nets[iNet].initW(wscale)
                iNet+=1

    def setReader(self,fname):
        ''' setReader(self,fname)
            set a reader for the SLFN class
        '''
        self.reader=Nto1READER(fname,self.numI,self.numO)


    def gotoHidden(self,A):
        ''' '''
        # we assume A is 3D: 
        #  - first index equal the number of targets, 
        #  - second index the sequence length
        #  - third index the input dimensions self.dimInput
        H=np.hstack((self.nets[i].gotoHidden(A) for i in range(self.numNets)))
        N=len(A)
        if self.hbias==True:
            Ones=np.ones((N,1))
            return  np.hstack((H,Ones))
        else:
            return np.array(H)



#--------------------------------------
def loadModel(fname,verbose=True):        
    ''' loadModel(fname) 
        => retuns the previously saved model
    '''
    gzip=False
    if fname[-3:]=='.gz' :
        os.system('gunzip '+fname)
        fname=fname[:-3]
        gzip=True
    f=open(fname,"r")
    obj = cPickle.load(f)
    f.close()
    if gzip : os.system('gzip '+fname)
    if verbose:
        sys.stdout.write("loaded ELM network from %s\n" % (fname))
        sys.stdout.write("   %s input neurones, %s output, %s hidden\n" % (str(obj.numI),str(obj.numO),str(obj.numH)))
        if 'win' in dir(obj) :
            sys.stdout.write("   %s actual input neurones, %s sliding window\n" % (str(obj.numID),str(obj.win)))
        sys.stdout.flush()
    obj.name=fname
    return obj
#======================================
        


# This is the learning class. It does all the training.

class ELM:
    def __init__(self,net):
        ''' __init__(self,net) 
            initialize an Extreme Learning Machine based on a network net
        '''
        self.net=net
        self.P=None  # used in the online learning
    # computes beta (the matrix of the weights) by loading everything, then inverting a matrix and getting it.
    # norm L2 (typical of ELM)
    def lrDirect(self,vectors, targets):
        ''' lrDirect(self,vectors, targets)
            learns the self.net.Beta using the direct method
        ''' 
        targets = np.array(targets)
#        if targets.shape[0] == 1:
#            targets = targets.transpose()
        H=self.net.gotoHidden(vectors)
        Hplus = np.linalg.pinv(H)
        self.net.Beta = np.dot(Hplus,targets)    
        if DEBUG:
            print "Direct ELM learning "
            print " H shape=",H.shape,"H+ shape=",Hplus.shape

    """
    def lrLasso(self,vectors, targets,alpha=0.01):
        ''' lrLasso(self,vectors, targets,alpha=0.01)
            learns the self.net.Beta using Lasso methods
        '''
        targets = np.array(targets)
        H=self.net.gotoHidden(vectors)
        #clf = linear_model.Lasso(alpha, max_iter=1000, tol=0.0001, fit_intercept=False) 
        clf = linear_model.Lasso(alpha, max_iter=1000, tol=0.01, fit_intercept=False) 
        clf.fit(H,targets)
        self.net.Beta =np.array(np.matrix(clf.coef_).T)
        if DEBUG:
            print "Lasso ELM learning with alpha=",alpha
            print " H shape=",H.shape
    
    # combo of norm L1 and L2
    def lrElasticNet(self,vectors, targets,alpha=0.01,l1_ratio=0.7):
        ''' lrElasticNet(self,vectors, targets,alpha=0.01,ro=)
            learns the self.net.Beta using Lasso methods
        '''
        targets = np.array(targets)
        H=self.net.gotoHidden(vectors)
        clf = linear_model.ElasticNet(alpha=alpha, l1_ratio=l1_ratio,fit_intercept=False,max_iter=2000, tol=0.01)
        clf.fit(H,targets)
        self.net.Beta =np.array(np.matrix(clf.coef_).T)
        if DEBUG:
            print "ElasticNet ELM learning with alpha=",alpha
            print " H shape=",H.shape
    """
    
    # improved version of lrDirect. it first multiplies by the transpose before inverting.
    def lrOrtProj(self, vectors, targets, lam=0.0):
        ''' lrOrtProj(self, vectors, targets, lam=None)
            learns the self.net.Beta using the Orthogonal Projection
            if lam > 0.0  H^T (I*lam + HH^T)^{-1} T is used  
        '''
        targets = np.array(targets)
#        if targets.shape[0] == 1:
#            targets = targets.transpose()
        H=self.net.gotoHidden(vectors)
        if lam > 0.0:
            Hplus =np.linalg.pinv((np.diag([lam]*H.shape[1])+np.dot(H.T,H)))
            self.net.Beta=np.dot(np.dot(Hplus,H.T),targets)
        else:
            Hplus = np.linalg.pinv(np.dot(H.T,H)) 
            #self.net.Beta=np.dot(np.dot(Hplus,H.T),targets)
            self.net.Beta=np.dot(Hplus,np.dot(H.T,targets))
        if DEBUG:
            print "Orthogonal Projection ELM learning "
            print " H shape=",H.shape,"H+ shape=",Hplus.shape

    # this reads from disc and compute each matrix element at a time. 
    # saves a lot of ram and it is very prone to be parallelized (not done yet).
    def lrHD(self,lam=0.0):
        '''
        lrHD(fp,reader, lam=0.0)
           fp=file handler (as in fp=open(filename))
           reader= a READER object
           if lam > 0.0 H^T (I*lam + HH^T)^{-1} T is used  
 
           Defining HtH=H^tH and HtT=H^t T.
           while fileReder returns a pairs of input/output vector pairs, x_k, t_k
                HtH += numpy.outer(g(x_k),g(x_k))
                HtT += numpy.outer(g(x_k),t_k)
           then build Beta=np.dot(np.linalg.pinv(HtH),)
        '''
        h=self.net.numHB
        reader=self.net.reader
        if lam > 0.0:
            HtH=np.diag([lam]*h)
        else:
            HtH=np.zeros((h,h))
        HtT=np.zeros((h,self.net.numO))
        reader.next()
# WHAT HAPPENS TO FIRST I/O
        #k=0
        #print 'first next'
        #sys.stdout.flush()
        while not reader.eof :
            x,t=reader.I,reader.O            
            #k+=1
            #print 'k=%d %s %s %d %d %s' % (k, str(type(x)),str(type(t)) ,len(x),len(t), str(reader.seqfile)) 
            #sys.stdout.flush()

            g=self.net.gotoHidden(np.array([x]))
            #if type(t) is list: print reader.seqfile,t
            #print 'HtH.shape,t.shape,len(x),g.shape',HtH.shape,t.shape,len(x),g.shape,reader.seqfile
            HtH+=np.outer(g,g)
            HtT+=np.outer(g, np.array(t))
            reader.next()
        #print 'done while'
        #sys.stdout.flush()
        Hplus=np.linalg.pinv(HtH) 
        self.net.Beta=np.dot(Hplus,HtT)
        if DEBUG:
            print "External Memory HD ELM learning "
            print " HtH shape=",HtH.shape,"HtH+ shape=",Hplus.shape
    # similar to the above but for sequence learning... 
    def lrSetHD(self,lam=0.0):
        '''
        lrSetHD(fp,reader, lam=0.0)
           differently from lrHD, here a set of vectors
           is returned after each call.
           fp=file handler (as in fp=open(filename))
           reader= a SeqREADER object
           if lam > 0.0 H^T (I*lam + HH^T)^{-1} T is used  (for the transpose and inverse stuff...)
 
           Defining HtH=H^tH and HtT=H^t T.
           while fileReder returns a pairs of input/output vector pairs, x_k, t_k
                HtH += numpy.outer(g(x_k),g(x_k))
                HtT += numpy.outer(g(x_k),t_k)
           then build Beta=np.dot(np.linalg.pinv(HtH),)
        '''
        h=self.net.numHB # load number of hidden neurones (plus 1 for bias if any)
        reader=self.net.reader # get the Reader object from the network. If it is a Seq network then this is a SeqREADER (and it must be)
        if not isinstance(reader,SeqREADER):
            print "NOT isinstance of SeqREADER"
            sys.exit()
        if lam > 0.0: # initialize matrix (this is at most number of Hidden* number of Hidden in size)
            HtH=np.diag([lam]*h)
        else:
            HtH=np.zeros((h,h))
        HtT=np.zeros((h , self.net.numO))
        reader.next()  # this loads the next (in this case the first) sequence
        while not reader.eof :
            for g,t in zip(self.net.gotoHidden(reader.I),reader.O): # a loop on the sliding window to fill the matrix
                HtH+=np.outer(g,g)
                HtT+=np.outer(g, np.array(t))
            reader.next() # this loads the next sequence
        Hplus=np.linalg.pinv(HtH)  # invert the matrix (use pseudo inverse as there is no guarantee the matrix is actually invertible..)
        self.net.Beta=np.dot(Hplus,HtT) # determine the Beta!!
        if DEBUG:
            print "External Memory HD ELM learning "
            print " HtH shape=",HtH.shape,"HtH+ shape=",Hplus.shape
    
    def lrSetRAM(self,input_output_list, lam=0.0):
        '''
        for sequence learning...
        lrSetRAM(input_output_list, lam=0.0)
           differently from lrHD, here a set of vectors
           is returned after each call.
           but differently from lrSetHD 
           this expects the input and output of all sequences to train on to be loaded in a list of tuples (or list o list)
           Therefore it doesn't read the file, which is very advantageous in the case one whant to repeat multiple trainings
           with different parameters (such as number of hidden neurones wscale or sliding window size)
           Clearly if there is a large number of trainins sequence it could be still convenient to use lrSetHD as it doesn't fill the RAM
           
           if lam > 0.0 H^T (I*lam + HH^T)^{-1} T is used  (for the transpose and inverse stuff...)
 
           Defining HtH=H^tH and HtT=H^t T.
           while fileReder returns a pairs of input/output vector pairs, x_k, t_k
                HtH += numpy.outer(g(x_k),g(x_k))
                HtT += numpy.outer(g(x_k),t_k)
           then build Beta=np.dot(np.linalg.pinv(HtH),)
        '''
        h=self.net.numHB # load number of hidden neurones (plus 1 for bias if any)
        # initialize matrices
        if lam > 0.0: # initialize matrix (this is at most number of Hidden* number of Hidden in size)
            HtH=np.diag([lam]*h)
        else:
            HtH=np.zeros((h,h))
        HtT=np.zeros((h , self.net.numO))
        
        for seq_in,seq_out in input_output_list :
            for g,t in zip(self.net.gotoHidden(seq_in),seq_out): # a loop on the sliding window to fill the matrix
                HtH+=np.outer(g,g)
                HtT+=np.outer(g, np.array(t))
        Hplus=np.linalg.pinv(HtH)  # invert the matrix (use pseudo inverse as there is no guarantee the matrix is actually invertible..)
        self.net.Beta=np.dot(Hplus,HtT) # determine the Beta!!
        if DEBUG:
            print "RAM External Memory HD ELM learning "
            print " HtH shape=",HtH.shape,"HtH+ shape=",Hplus.shape

 
    # don't use, it is slower and approximate (very advertised in the literature though).
    def lrOnline(self,vectors, targets, first=True):
        ''' lrOnline(self,vectors, targets, first=True)
            online learning
        '''
        targets = np.array(targets)
        H=self.net.gotoHidden(vectors)
        if first==True:
            self.P=np.linalg.pinv(np.dot(H.T,H))
            self.net.Beta=np.dot(np.dot(self.P,H.T),targets)
        else:
            P=self.P
            # Inv=H^T (I+ H*P*H^T)^{-1} 
            HtInvH=np.dot(np.dot(H.T,np.linalg.pinv(np.identity(H.shape[0])+np.dot(H,np.dot(P,H.T)))),H)
            # P= P  -P*HtInvH*P = 
            self.P-= np.dot(P,np.dot(HtInvH,P))
            # beta= beta + P* H^T (T - H * brta)
            self.net.Beta+=np.dot(self.P,np.dot(H.T,(targets-np.dot(H,self.net.Beta))))


    def lrIter(self,vectors, targets, first=True):
        ''' lrIter(self,vectors, targets, first=True)
            online learning one step at a time
        '''
        targets = np.array(targets)
#        if target.shape[0] == 1:
#            target = target.transpose()
        H=self.net.gotoHidden(vectors)
        if first==True:
            self.P=np.linalg.pinv(np.dot(H.T,H))
            self.net.Beta=np.dot(np.dot(self.P,H.T),targets)
        else:
            P=self.P
            D=1+np.dot(np.dot(H,P),H.T)
            self.P+=-np.dot(P,np.dot(np.outer(H,H),P))/D
            self.net.Beta+=np.dot(self.P,np.dot(H.T,(targets-np.dot(H,self.net.Beta))))
            #self.net.Beta+=np.dot(P,np.dot(H,(targets-np.dot(H,self.net.Beta))))





    def save(self, fname):
        ''' obj.save(fname) 
            save obj.net on the file fname
        '''
        f=open(fname,"w")
        cPickle.dump(self.net,f)
        f.close()



#======================================

if __name__=='__main__':
    DEBUG=True
    N=3
    I=2
    H=2
    O=1
    T=np.array([[-3],[3.0],[3.0]]) # or
    X=np.array([[0.0,0.0],[1.0,0.0],[0.0,1.0]])
    net=SLFN(I,H,O)
    net.inintWeights(0.01)       
    elmModel = ELM(net)
    elmModel.lrOrtProj(X, T)
    print "Taining"
    P=net.actFun(net.propagate(X))
    for i in range(len(X)):
        print X[i],P[i]
    print "Testing"
    Y=np.array([[0.0,0.0],[1.0,0.0],[0.0,1.0],[1.0,1.0]])
    P=net.actFun(net.propagate(Y))
    for i in range(len(Y)):
        print Y[i],P[i]

