'''
    Copyright (C) 2014 Pietro Sormanni

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

    contacts: Pietro Sormanni
              e-mail: ps589@cam.ac.uk 
              Dept. of Chemistry
              University of Cambridge
              Cambridge, UK
              
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it
    under the condition of preserving this text
'''


WEBSERVER_MODE=False

import sys,os
import numpy
import sysconfig
#from traceback import print_exc


import PyELM
import s2D
from traceback import print_exc



module_path=str(s2D.__file__)
module_path=os.path.abspath( module_path[:module_path.rfind('/')])
if module_path[-1]!='/' : module_path+='/'

default_c_compiler = sysconfig.get_config_vars('CC')[0]
if default_c_compiler==None : default_c_compiler='cc'

default_parser_executable=module_path+'chkparse'






if WEBSERVER_MODE :
    sys.path=['/home/ps589/.local/lib/python2.7/site-packages/','/home/ps589/.local/lib/python2.7/site-packages/six-1.8.0-py2.7.egg', '/home/ps589/.local/lib/python2.7/site-packages/mock-1.0.1-py2.7.egg','/home/ps589/.local/lib/python2.7/site-packages/distribute-0.6.28-py2.7.egg','/home/ps589/.local/lib/python2.7/site-packages/nose-1.3.4-py2.7.egg','/home/ps589/.local/lib/python2.7/site-packages/matplotlib-1.4.1-py2.7-linux-x86_64.egg']+sys.path
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator,AutoMinorLocator
    matplotlib.use('Agg')
    can_plot=True
else :
    try : # check if plots can be generated
        import matplotlib,distutils
        matplot_vers=  distutils.version.LooseVersion(matplotlib.__version__).version
        if matplot_vers<[1,4] : 
            sys.stderr.write("*Warning*: your matplotilb version is %s but s2D requires >= 1.4 to produce plots.\n  You can still run s2D to make predictions, but you won't be able to use it to make plots.\n" % (str(matplotlib.__version__)))
            can_plot=False
        else : can_plot=True
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator,AutoMinorLocator
        plt.rcdefaults()
        plt.rc('figure',facecolor='white')
            
    except Exception :
        #print_exc(file=sys.stderr)
        sys.stderr.write("\n***WARNING*** matplotlib module not found. attempting to plot the results will generate ERRORS..\n\n")
        sys.stderr.flush()
        can_plot=False
        pass


try : # get default parameters from file
    #from pkg_resources import Requirement, resource_filename
    #filename = resource_filename(Requirement.parse("s2D"),"s2D_parameters.txt")
    #print filename
    #foo_config = resource_string(__name__, 'foo.conf')
    #if os.path.isfile(filename) : #
    if os.path.isfile(module_path+'s2D_parameters.txt') :
        parameter_filename=module_path+'s2D_parameters.txt'
        #imported_parameter_class=s2D.s2D_parameters()
        #imported_parameter_class.read('s2D_parameters.txt')
    else :
        parameter_filename=None 
        #imported_parameter_class=None
except Exception :
    parameter_filename=None
    #imported_parameter_class=None
    #print_exc(file=sys.stderr)
    pass





def convert_to_number(string, force_float=False):
    '''
    this function check if a string is an int or a float and it returns a tuple in the form
    converted_string,bool. Bool is True if the sting has been converted, False if the  string is still in string format.
    the function is quite slow
    '''
    if force_float :
        try :
            return float(string),True
        except ValueError :
            return string,False
    try :
        return int(string),True
    except ValueError :
        try :
            return float(string),True
        except ValueError :
            return string,False
"""
probably to delete:
def get_numbers_from_string(string, force_float=False):
    '''
    from a string (like a file name or such) it extracts all the possible numbers and return them in a list
    it does NOT process negative values (will be read as positives)
    '''
    candidates=[]
    reading=False
    for ch in string:
        if ch.isdigit() or ch=='.' :
            if reading :
                candidates[-1]+=ch 
            else :
                candidates+=[ch]
                reading=True
        else :
            reading=False
    numbers=[]
    for ca in candidates :
        if ca[-1]=='.' : ca=ca[:-1]
        if ca[0]=='.' and len(ca)>1 : ca=ca[1:]
        ca,ok=convert_to_number(ca, force_float=force_float)
        if ok :
            numbers+=[ca]
    return numbers
"""

res_closness_blosum62 = ['V', 'I', 'M', 'L', 'F', 'W', 'Y', 'H', 'N', 'S', 'P', 'T', 'C', 'A', 'G', 'D', 'E', 'Q', 'K', 'R']

# convert residue to neurons without using blast
def residue_to_input(res, numb_neurones=1,out_in_range=(-10.,10.), use_rank=res_closness_blosum62,nres=20):
    if res not in use_rank: 
        raise Exception('*ERROR* residue '+str(res)+' not recognized ')
        return None
    if numb_neurones==1 :
        inp = out_in_range[0] + use_rank.index(res)*(out_in_range[1]-out_in_range[0])/(nres*1.)
        return inp
    else :
        inp=[0]*nres
        inp[ use_rank.index(res) ]=1
        return inp

def linear_correction(out_from_networks,angular_coefficients=numpy.array([1.240177580,1.27546432,1.397699985]) , intercepts=numpy.array([-0.070890685,-0.0577194695,-0.198067935]) ):
    '''
    applies a linear correction. THIS IS NOT USED IN THE CURRENT VERSION
      the default parameters implies an output of 3 number per position of the Seq network.
      helix: y= 1.240177580 *x -0.070890685
      beta : y= 1.27546432  *x -0.057719469498999997
      coil : y= 1.397699985 *x -0.198067935
    '''
    return out_from_networks*angular_coefficients + intercepts


def constrain_in_range(numpy_array,Min=0., Max=1.,shift_by=0.01):
    numpy_array[numpy_array<Min]=Min+shift_by # the addition of shift_by is crucial if you normalize after
    numpy_array[numpy_array>Max]=Max-shift_by
    return numpy_array

def normalize_numpy(numpy_array, scaling_factor=1., constrain_first=True,Min=0., Max=1.,shift_by=0.01):
    '''
    assumes that numpy_array is 2 dimensional, each row is the output of one residue
    and each column represents a property. It normalizes the properties so that they sum to one
    the input array is changed, no copy is made!
    if constrain_first it first put the entries in 0 ,1 
    '''
    if constrain_first : numpy_array=constrain_in_range(numpy_array,Min=Min, Max=Max,shift_by=shift_by)
    den=numpy_array.sum(axis=1)/scaling_factor
    for i in range(numpy_array.shape[1]) :
        numpy_array[:,i]/=den
    del den
    return numpy_array


def linear_correct_and_normalize(out_from_networks):
    '''
    it first applies a linear correction and then normalizes the output
    '''
    return normalize_numpy( linear_correction(out_from_networks) , constrain_first=True)



        
"""
PLOT FUNCTION
  these functions are used to produce a plot of the results. (Besides plot_s2D_results all others are more general and can be used for many purposes).
"""
class cycle_list(list):
    def __init__(self,l=[]):
        list.__init__(self,l)
    def __getitem__(self,y) :
        #x.__getitem__(y) <==> x[y]
        if y>=len(self) : 
            y=y%len(self)
        return list.__getitem__(self,y)

color_palette = cycle_list([(0.25098039215686274, 0.4470588235294118, 0.792156862745098), (0.4235294117647059, 0.6831372549019608, 0.24313725490196078), (0.8470588235294118, 0.7507843137254902, 0.16784313725490197), (0.8196078431372549, 0.5333333333333333, 0.15294117647058825), (0.7764705882352941, 0.29411764705882354, 0.10980392156862745), (0.4549019607843137, 0.3254901960784314, 0.6509803921568628)])


def plot_s2D_results(s2D_output,sequences=None,ss_kind_str=None,seq_names=None,title=None,frame=None, print_all_sequence=True,dont_plot_coil=True,plot_coil_as_bar=False,coil_index=2,coil_color=None,plot_legend=True,y_range=(0,1.00001),figure_and_axts_tuple=None,bar=True,dpi=300,max_res_per_plot=None,color_palette=color_palette \
                     ,xlabel='Residue', ylabel='Secondary structure population',legend_location='upper right',start_rescount=1,y_major_tick_every=0.5,y_minor_tick_every=0.1,x_major_tick_every=None,x_minor_tick_every=None, legend_size=None,figure_size=None,save=True,show=False,**kwargs):
    '''
    plot the results in a graph, can read the output file produced by the s2D class.
    out_tag not used yet
    return fig,axt
    '''
    lw=0.25
    if not bar : lw=3.
    if type(sequences) is not list or type(sequences) is not tuple :
        sequences=[sequences]
    if type(seq_names) is not list or type(seq_names) is not tuple :
        seq_names=[seq_names]
    if type(s2D_output) is str :
        seq_names,sequences,s2D_output,ss_kind_str = read_output(s2D_output)
    elif type(s2D_output) is not list : # One list per sequence # and not ( hasattr(s2D_output,'shape') and len(s2D_output.shape)>1 ): # if it is not a list and not an (at least) bidimensional numpy array
        s2D_output=[s2D_output]
    save_plot=None
    
    for j,out in enumerate(s2D_output) : # loop on sequences
        #if sequences[j]==None :
        #    sequences[j]=' '*len(out)
        legend=['Helix','Strand','Coil','Polyproline-II'][:len(out.T)]
        if seq_names[j]==None : seq_names[j]=''
        if title==None and seq_names[j]!=None and seq_names[j]!='': title ='s2D Prediction '+seq_names[j]
        if dont_plot_coil :
            cols=color_palette[:coil_index]+color_palette[coil_index+1:]
            fig,axt=plot_seq_profile(sequences[j], list(out.T[:coil_index])+list(out.T[coil_index+1:]),annotation_string=ss_kind_str , bar=bar, start_rescount=start_rescount, xlabel=xlabel, ylabel=ylabel, title=title,frame=frame,y_range=y_range,zygg_like_lines=False, print_all_sequence=print_all_sequence, color=cols, max_res_per_plot=max_res_per_plot \
                                             ,y_major_tick_every=y_major_tick_every,y_minor_tick_every=y_minor_tick_every,x_major_tick_every=x_major_tick_every,x_minor_tick_every=x_minor_tick_every, show=False, linewidth=lw,figure_and_axts_tuple=figure_and_axts_tuple,figure_size=figure_size, save=None,**kwargs)
            if plot_legend :
                if bar : add_custom_legend(legend[:coil_index]+legend[coil_index+1:],facecolors=cols,edgecolors=None,figure=fig, legend_location=legend_location, legend_size=legend_size,linewidth=[lw]*len(cols), frame_facecolor=(1.,1.,1.,0.7))
                else : add_custom_legend(legend[:coil_index]+legend[coil_index+1:],facecolors=None,edgecolors=cols,figure=fig, legend_location=legend_location, legend_size=legend_size,linewidth=[lw]*len(cols),frame_facecolor=(1.,1.,1.,0.7))
        else :
            
            if plot_coil_as_bar :
                cols=color_palette[:len(out.T)]
                fig,axt=plot_seq_profile(sequences[j], out.T,annotation_string=ss_kind_str , bar=plot_coil_as_bar, start_rescount=start_rescount, xlabel=xlabel, ylabel=ylabel, title=title,frame=frame,y_range=y_range,zygg_like_lines=False, print_all_sequence=print_all_sequence, color=cols,max_res_per_plot=max_res_per_plot ,y_major_tick_every=y_major_tick_every,y_minor_tick_every=y_minor_tick_every,x_major_tick_every=x_major_tick_every,x_minor_tick_every=x_minor_tick_every, show=False, linewidth=lw,figure_and_axts_tuple=figure_and_axts_tuple,figure_size=figure_size, save=None,**kwargs)
                face_cols,edge_cols,lws = cols,[None]*len(cols),[lw]*len(cols)
            else :
                if 'Coil' in legend : 
                    legend.remove('Coil')
                    legend+=['Coil']
                if coil_color!=None : 
                    coil_col=coil_color
                    cols=color_palette[:len(out.T)-1]
                else : 
                    coil_col=color_palette[coil_index]
                    cols=(color_palette[:coil_index]+color_palette[coil_index+1:])[:len(out.T)-1]
                if bar : lw_coil=lw+2
                else : lw_coil=lw
                fig,axt=plot_seq_profile(sequences[j], list(out.T[:coil_index])+list(out.T[coil_index+1:]),annotation_string=ss_kind_str , bar=bar, start_rescount=start_rescount, xlabel=xlabel, ylabel=ylabel, title=title,frame=frame,y_range=y_range,zygg_like_lines=False, print_all_sequence=print_all_sequence, color=cols,max_res_per_plot=max_res_per_plot ,y_major_tick_every=y_major_tick_every,y_minor_tick_every=y_minor_tick_every,x_major_tick_every=x_major_tick_every,x_minor_tick_every=x_minor_tick_every, show=False, linewidth=lw,figure_and_axts_tuple=figure_and_axts_tuple,figure_size=figure_size, save=None,**kwargs)
                fig,axt=plot_seq_profile(sequences[j], out.T[coil_index],annotation_string=ss_kind_str, bar=False, start_rescount=start_rescount, xlabel=xlabel, ylabel=ylabel, title=title,frame=frame,y_range=y_range,zygg_like_lines=False, print_all_sequence=print_all_sequence, color=coil_col,max_res_per_plot=max_res_per_plot, y_major_tick_every=y_major_tick_every,y_minor_tick_every=y_minor_tick_every,x_major_tick_every=x_major_tick_every,x_minor_tick_every=x_minor_tick_every, show=False, linewidth=lw_coil,figure_and_axts_tuple=(fig,axt),figure_size=figure_size, save=None,**kwargs)
                if bar :
                    face_cols=cols+[None] #color_palette[:coil_index]+[None]+color_palette[coil_index+1:]
                    edge_cols,lws = [None]*(len(legend)-1)+[coil_col],[lw]*(len(legend)-1)+[lw_coil] #[ None if c!=None else color_palette[ji] for ji,c in enumerate(face_cols)], [lw if c!=None else lw_coil for ji,c in enumerate(face_cols) ]
                    #print len(cols),len(face_cols),face_cols
                else : face_cols,edge_cols,lws=[None]*len(legend),cols+[coil_col],[lw]*len(legend)#[None]*len(legend) , color_palette[:len(legend)],[lw]*len(legend)
            if plot_legend :
                add_custom_legend(legend,facecolors=face_cols,edgecolors=edge_cols,figure=fig, legend_location=legend_location, legend_size=legend_size,linewidth=lws,frame_facecolor=(1.,1.,1.,0.7))
            #cols=color_palette[:2]+color_palette[3:4]
            #fig,axt=plot_seq_profile(sequences[j], out.T,annotation_string=ss_kind_str , bar=bar, start_rescount=start_rescount, xlabel=xlabel, ylabel=ylabel,frame=frame, title=title,y_range=y_range,zygg_like_lines=False, print_all_sequence=print_all_sequence, color=cols \
            #                                 ,y_major_tick_every=y_major_tick_every,y_minor_tick_every=y_minor_tick_every,x_major_tick_every=x_major_tick_every,x_minor_tick_every=x_minor_tick_every, show=False, linewidth=lw,figure_and_axts_tuple=figure_and_axts_tuple,figure_size=figure_size, save=None)
            #if plot_legend :
            #    if bar : add_custom_legend(['Helix','Strand','Coil'],facecolors=cols,edgecolors=None,figure=fig, legend_location=legend_location, legend_size=legend_size)
            #    else : add_custom_legend(['Helix','Strand','Coil'],facecolors=None,edgecolors=cols,figure=fig, legend_location=legend_location, legend_size=legend_size)
        if type(save) is str : save_plot=save
        elif save==True : save_plot=seq_names[j].replace('|','_').replace(' ','_').replace('/','_').replace(':','_')+'s2D_plot.png'
        if save!=None and save!=False:
            if 'png' in save_plot : transparent=True
            else : transparent=False 
            fig.savefig(save_plot, dpi=dpi,bbox_inches="tight",transparent=transparent)
    if show : 
        plt.show(block=False)
    return fig,axt


default_parameters = {
    'hgrid':True, 
    'vgrid':True, 
    'frame':True,
    'all_tight':False,
    'seq_max_res_per_plot':200}
default_error_bars = {
    'capsize':4, 
    'capthick':1., 
    'elinewidth':1.}
text_sizes = {
    'value_labels':18, 
    'xlabels':18, 
    'xlabels_many':'small', 
    'xlabel':22, 
    'ylabel':22, 
    'title':24, 
    'legend_size':12}
publication = {
    'value_labels':22, 
    'xlabels':22, 
    'xlabels_many':15, 
    'xlabel':30, 
    'ylabel':30, 
    'title':30, 
    'legend_size':18}
# if you set the size for 'default' all the figures will come out of that size disregarding their type. Otherwise you can change the figure size for each type (key in dictionary)
default_figure_sizes={
'all_tight':False ,
'use_cm':False ,
'dpi':300,
'default':None ,
'sequence':(14.2,8)
}

def set_publish(all_same_figure_size=False, thick_ticks=True, axis_tickness=True, no_grids=True, text_sizes=text_sizes, publication=publication, default_figure_sizes=default_figure_sizes):
    default_error_bars['capsize'] = 8
    default_error_bars['capthick'] = 2
    default_error_bars['elinewidth'] = 2
    for k in publication:
        text_sizes[k] = publication[k] # with an = or with .copy it does not work
    default_parameters['seq_max_res_per_plot']=100
    default_parameters['all_tight']=True
    plt.rc('xtick', labelsize=text_sizes['xlabels'])
    plt.rc('ytick', labelsize=text_sizes['xlabels'])
    plt.rc('ytick.major', width=1.5, size=6)
    plt.rc('ytick.minor', width=1., size=3)
    plt.rc('xtick.major', width=1.5, size=6)
    plt.rc('xtick.minor', width=1., size=3)
    default_figure_sizes['all_tight'] = True
    if no_grids:
        for p in ['vgrid','hgrid'] :
            default_parameters[p] = False
    default_parameters['frame']=True # frame on for s2D
    
    if thick_ticks :
        plt.rc('ytick.major', width=2.5, size=10)
        plt.rc('ytick.minor', width=2, size=6)
        plt.rc('xtick.major', width=2.5, size=10)
        plt.rc('xtick.minor', width=2, size=6)
    if axis_tickness:
        plt.rc('axes', linewidth=2, edgecolor='black')
    if all_same_figure_size != False: # I could just change default_figure_sizes['default'] but this is safer
        if type(all_same_figure_size) is tuple:
            for s in default_figure_sizes:
                if s == 'all_tight' or s == 'use_cm' or s == 'dpi':
                    continue
                default_figure_sizes[s] = all_same_figure_size
        else:
            for s in default_figure_sizes:
                if s == 'all_tight' or s == 'use_cm' or s == 'dpi':
                    continue
                default_figure_sizes[s] = 10, 10
    return


def plot_seq_profile(sequence, profile, annotation_string=None,use_subplot=True,bar=False,bar_sep=0.2,log_scale=False,avoid_scientific_notation=True,max_res_per_plot=None,stacked=False,start_rescount=1, label='',xlabel='Residue',ylabel='Score (a.u.)',title=None,zygg_like_lines=True,hline=0,vgrid=None,frame=None ,print_all_sequence=True,color=None, show=True\
                             ,yerr=None,y_range=None,y_major_tick_every=None,y_minor_tick_every=None,x_major_tick_every=None,x_minor_tick_every=None,ls='-', linewidth=1,marker='',markerfacecolor=True,markeredgecolor=True, markersize=18,upper_label_rotation='horizontal',legend_location='upper right',legend_size=None,figure_size=None, figure_and_axts_tuple=None, save='') :
    '''
    figure_and_axts_tuple can be given to superimpose
    '''
    if frame==None : frame=default_parameters['frame']
    #if hgrid==None : hgrid=default_parameters['hgrid']
    if vgrid==None : vgrid=default_parameters['vgrid']
    if max_res_per_plot==None : max_res_per_plot=default_parameters['seq_max_res_per_plot']
    if figure_size==None : 
        if default_figure_sizes['default']==None : 
            figure_size = default_figure_sizes['sequence']
        else : figure_size = default_figure_sizes['default']
        if type(sequence) is str and len(sequence)>900 :
            figure_size=(figure_size[0],min([20,int(len(sequence)/10)]))
            #print 'DEB figure_size',figure_size
    #if default_figure_sizes['use_cm'] : figure_size= ( cmToinch(figure_size[0]),  cmToinch(figure_size[1]))
    if legend_size==None : legend_size=text_sizes['legend_size']
    plt.rc('xtick', labelsize=text_sizes['xlabels'])
    plt.rc('ytick', labelsize=text_sizes['xlabels'])
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
    
    #if default_figure_sizes['all_tight'] : 
    #    #if y_major_tick_every==None and x_major_tick_every==None :
    #    plt.locator_params(axis='both', tight=None, nbins=5)
    if hasattr(profile[0],'__len__') : 
        ismulti=True
        if label!='' and label!=None and type(label) is not list and type(label) is not tuple : label=len(profile)*[label]
        if type(ls) is str or ls==None : ls=[ls]*len(profile)
        if type(linewidth) is int or type(linewidth) is float or linewidth==None : linewidth=[linewidth]*len(profile)
        lengths=[]
        for p in profile :
            if len(p) not in lengths : lengths+=[len(p)]
        maxlength=max(lengths)
        if len(lengths)>1 :
                sys.stderr.write("**WARNING in plot_seq_profile() given profiles of different lengths (found %s). Appending zeros at end of shorter profiles!\n\n" % (str(lengths)) )
                for j,p in enumerate(profile) :
                    if len(p)< maxlength :
                        profile[j]=list(p)+[0.]*(maxlength-len(p))
        profile=numpy.array(profile)
        if yerr!=None : yerr=numpy.array(yerr)
        prof=profile[0]
        if y_range == None :
            Min,Max= profile.min(),profile.max()
            while hasattr(Min,'__len__') : Min=min(Min)
            while hasattr(Max,'__len__') : Max=max(Max)
    else : 
        ismulti=False
        prof=profile
        if y_range == None :
            Min=min(profile)
            Max=max(profile)
    if y_range == None :
        ymin=int(Min -1.)
        ymax=int(Max +1.)
    else :
        ymin,ymax=y_range
    if sequence!=None and len(sequence)!=len(prof) :
        sys.stderr.write('**WARNING** in plot_seq_profile() len(sequence)!=len(profile) %d!=%d\n' % (len(sequence),len(prof)))
    #if ecolor==None : ecolor=color
    if type(markeredgecolor) is bool and markeredgecolor==True : markeredgecolor=color
    if type(markerfacecolor) is bool and markerfacecolor==True : markerfacecolor=color
    if type(color) is list or isinstance(color,cycle_list) :
        if type(markerfacecolor) is not list and not isinstance(markerfacecolor,cycle_list) : markerfacecolor= cycle_list([markerfacecolor]*len(color))
        if type(markeredgecolor) is not list and not isinstance(markeredgecolor,cycle_list) : markeredgecolor= cycle_list([markeredgecolor]*len(color))
        #if type(ecolor) is not list and not isinstance(ecolor,cycle_list) : ecolor=cycle_list([ecolor]*len(color))
        
    if use_subplot  :
        
        if figure_and_axts_tuple!=None : 
            fig,axt = figure_and_axts_tuple
            n_profs=len(axt)
            do_tight=False
        else :
            if len(prof)%int(max_res_per_plot) > 0 : add=1
            else : add=0.001 # in case of rounding errors from the float conversion
            n_profs=max([1 , int(len(prof)/float(max_res_per_plot)+add)]) # up to 199 residues per plot 
            fig,axt = plt.subplots(n_profs, sharey=True,figsize=figure_size) # axt is a tuple of n_profs size Use , sharex=True to share x axis
            if n_profs==1 : axt=(axt,)
            do_tight=True
        # determine the number of residues per subplot
        res_per_plot=len(prof)/n_profs
        rem=len(prof)%n_profs
        pzise=[res_per_plot]*n_profs
        j=0
        while rem>0 :
            pzise[j]+=1
            rem -=1
            j +=1
        start=0
        
        line_styles=[]
        for j,nres in enumerate(pzise) :
            xlim_m=start+start_rescount-0.5
            xlim_M=start+nres+start_rescount-0.5
            
            if default_figure_sizes['all_tight'] : 
                if x_minor_tick_every==None : axt[j].xaxis.set_minor_locator(AutoMinorLocator(n=2))
                if y_minor_tick_every==None : axt[j].yaxis.set_minor_locator(AutoMinorLocator(n=2))
            
            if ismulti : 
                to_plot=profile[:,start:start+nres]
                if yerr!=None : pl_yerr=[ a[start:start+nres] if a!=None else None for a in yerr ]
                else : pl_yerr=[None]*len(to_plot)
            else : 
                to_plot=profile[start:start+nres]
                if yerr!=None : pl_yerr=yerr[start:start+nres]
                else : pl_yerr=None
            
            #if log_scale :
            #    if j==0 and y_range!=None : print "WARNING log_scale is overwriting y_range"
            #    _,to_plot,y_range =logscale(axt[j], entries=to_plot,add_one=True,add_zero=True)        
                
            if bar :
                if ismulti :
                    if stacked :   
                        bottom=numpy.zeros(nres)
                        sep,bar_width = bar_sep, (1.-bar_sep)
                        left=numpy.array(range(start+start_rescount,start+nres+start_rescount))-0.5+sep/2.
                    else :
                        sep, bar_width= float(bar_sep)/(len(profile)+1), (1.-bar_sep)/len(profile)# +1 is there so that bar groups will be separated by 2*sep 
                        bottom=None
                        left = numpy.array(range(start+start_rescount,start+nres+start_rescount))-0.5+sep 
                    for i,prof in enumerate(to_plot) :
                        if type(color) is list or isinstance(color,cycle_list)  : l= axt[j].bar(left, prof,bottom=bottom,yerr=pl_yerr[i] , width=bar_width,linewidth=linewidth[i],color=color[i])
                        else : l=axt[j].bar(left, prof, width=bar_width,yerr=pl_yerr[i], linewidth=linewidth[i],color=color)
                        if stacked : bottom+=numpy.array(prof)
                        else : left+=sep+bar_width
                        if start==0 : line_styles+=[l]
                        del l
                    
                else :    
                    sep,bar_width = bar_sep, (1.-bar_sep)
                    left=numpy.array(range(start+start_rescount,start+nres+start_rescount))-0.5+sep/2.
                    l=axt[j].bar(left,to_plot,width=bar_width,yerr=pl_yerr, linewidth=linewidth,color=color)
                    if start==0 : line_styles+=[l]
                    del l
            else :
                if ismulti :
                    x_pos=range(start+start_rescount,start+nres+start_rescount)
                    for i,prof in enumerate(to_plot) :
                        if color==None : l=axt[j].errorbar(x_pos, prof, linewidth=linewidth[i],yerr=pl_yerr[i],ls=ls[i],elinewidth=default_error_bars['elinewidth'], capsize=default_error_bars['capsize'], capthick=default_error_bars['capthick']\
                                                                      , marker=marker,markersize=markersize,markeredgecolor=markeredgecolor,markerfacecolor=markerfacecolor)
                        elif type(color) is list or isinstance(color,cycle_list)  :
                            l=axt[j].errorbar(x_pos, prof,yerr=pl_yerr[i], linewidth=linewidth[i],ls=ls[i],color=color[i],elinewidth=default_error_bars['elinewidth'], capsize=default_error_bars['capsize'], capthick=default_error_bars['capthick']\
                                                         , marker=marker,markersize=markersize,markeredgecolor=markeredgecolor[i],markerfacecolor=markerfacecolor[i])
                        else : l=axt[j].errorbar(x_pos, prof, linewidth=linewidth[i],yerr=pl_yerr[i],ls=ls[i],color=color,elinewidth=default_error_bars['elinewidth'], capsize=default_error_bars['capsize'], capthick=default_error_bars['capthick']\
                                                            , marker=marker,markersize=markersize,markeredgecolor=markeredgecolor,markerfacecolor=markerfacecolor)
                        if start==0 : line_styles+=[l]
                        del l
                else :
                    if color==None :l=axt[j].errorbar(range(start+start_rescount,start+nres+start_rescount),to_plot, linewidth=linewidth,yerr=pl_yerr,ls=ls,elinewidth=default_error_bars['elinewidth'], capsize=default_error_bars['capsize'], capthick=default_error_bars['capthick'], marker=marker,markersize=markersize,markeredgecolor=markeredgecolor,markerfacecolor=markerfacecolor)
                    else :l=axt[j].errorbar(range(start+start_rescount,start+nres+start_rescount),to_plot, linewidth=linewidth,yerr=pl_yerr,ls=ls,color=color,elinewidth=default_error_bars['elinewidth'], capsize=default_error_bars['capsize'], capthick=default_error_bars['capthick'], marker=marker,markersize=markersize,markeredgecolor=markeredgecolor,markerfacecolor=markerfacecolor)
                    if start==0 : line_styles+=[l]
                    del l
            #axt[j].set_xlim(start+start_rescount,start+nres+start_rescount-1)
            
            if hline!=None :
                if type(hline) is not list or type(hline) is not tuple : hline=[hline]
                for ypos in hline :
                    axt[j].axhline(ypos,color='black',ls='-') #use thick line in this case, it represent the axis)
            if zygg_like_lines!=False :
                if type(zygg_like_lines) is not list or type(zygg_like_lines) is not tuple :
                    zygg_like_lines=(-1,1)
                for ypos in zygg_like_lines :
                    axt[j].axhline(ypos,color='black',ls='--',lw=0.5) #use thick line in this case, it represent the axis)
                        
            if log_scale : 
                axt[j].set_yscale('symlog',basey=10)
                if avoid_scientific_notation:
                    yticks=axt[j].yaxis.get_majorticklocs()
                    xlab=[ 10**i for i in xrange(len(yticks))]
                    axt[j].set_yticklabels(xlab,rotation='horizontal',verticalalignment='center',horizontalalignment='right',fontsize=text_sizes['xlabels'])
            
            
            if y_minor_tick_every!=None :
                yminorLocator   = MultipleLocator(y_minor_tick_every)
                #for the minor ticks, use no labels; default NullFormatter
                axt[j].yaxis.set_minor_locator(yminorLocator)
            if y_major_tick_every!=None :
                ymajorLocator   = MultipleLocator(y_major_tick_every)
                axt[j].yaxis.set_major_locator(ymajorLocator)
            if x_minor_tick_every!=None :
                xminorLocator   = MultipleLocator(x_minor_tick_every)
                #for the minor ticks, use no labels; default NullFormatter
                axt[j].xaxis.set_minor_locator(xminorLocator)
            if x_major_tick_every!=None :
                xmajorLocator   = MultipleLocator(x_major_tick_every)
                axt[j].xaxis.set_major_locator(xmajorLocator)
            #axt[j].set_ylim(ymin, ymax )
            
            xticks=axt[j].xaxis.get_majorticklocs()
            xticks=map(float,xticks)
            sp=1.*nres
            to_remove=[]
            for x in xticks :
                if abs((x-start+start_rescount)/sp)<0.1 : to_remove.append(x)
                elif abs((start+nres+start_rescount-1-x)/sp)<0.1 : to_remove.append(x)
            for x in to_remove :
                xticks.remove(x)
            xticks+=[start+start_rescount,start+nres+start_rescount-1]
            axt[j].set_xticks(xticks)
            axt[j].set_xlim(xlim_m,xlim_M)
            #handle_grid( axt[j] , vgrid=False , hgrid=hgrid ) # custom vgrid for this plot
            if vgrid :
                if type(vgrid) is list or type(vgrid) is tuple :
                    for vl in vgrid : axt[j].axvline(vl,color='grey',ls=':')
                else :
                    for count in range(start+start_rescount,start+nres+start_rescount) :
                        if type(vgrid) is int and count%vgrid==0:
                            axt[j].axvline(count,color='grey',ls=':')
                        elif count%10==0:
                            axt[j].axvline(count,color='grey',ls=':')
                    # axt[j].annotate(sequence[count-1], xy=(count,ymin), xytext=(0, -5),rotation=rotation, textcoords='offset points', va='top', ha='center',size='small')
            ax2=None
            if ( sequence!=None or annotation_string!=None ) and frame!=False:
                ax2=axt[j].twiny()
                ax2.set_xlim(axt[j].get_xlim())
                if print_all_sequence==True : ju=1
                elif type(print_all_sequence) is int : ju=print_all_sequence
                else : ju=3
                ax2.set_xticks(range(start+start_rescount,start+nres+start_rescount),minor=True) #
                ax2.set_xticks(range(start+start_rescount,start+nres+start_rescount,ju))
                ax2_ticks=ax2.get_xticks()
                if sequence!=None : 
                    ax2.set_xticklabels(sequence[start:start+nres:ju],rotation=upper_label_rotation,verticalalignment='bottom',fontsize=text_sizes['xlabels_many'])
                    #an=list(sequence[start:start+nres:ju])
                    #for ja,xt in enumerate(ax2_ticks[::ju]) :    
                    #    axt[j].annotate( an[ja], (xt,ymax),(0,5), xycoords='data' \
                    #    , size=text_sizes['xlabels_many'],textcoords = 'offset points', ha = 'center', va = 'bottom' )
                if annotation_string!=None :
                    an=list(annotation_string[start:start+nres:ju])
                    if len(an)>=len(ax2_ticks) :
                        for ja,xt in enumerate(ax2_ticks) :
                            axt[j].annotate( an[ja], (xt,ymax),(0,-5), xycoords='data' \
                            , size=text_sizes['xlabels_many'],textcoords = 'offset points', ha = 'center', va = 'top' )
                    else :
                        sys.stderr.write("Warn in plot_seq_profile. Not plotting annotation_string as length of processed str is larger than number of ticks on top axis. [%d %d]\n" % (len(an),len(ax2_ticks)))
            if not frame :
                # this remove top and right axis
                print 'Removing frame'
                axt[j].spines["right"].set_visible(False)
                axt[j].spines["top"].set_visible(False)
                axt[j].get_xaxis().tick_bottom() # ticks only on bottom axis
                axt[j].get_yaxis().tick_left() # ticks only on left axis
                if ax2!=None :
                    ax2.spines["right"].set_visible(False)
                    ax2.spines["top"].set_visible(False)
                    ax2.get_xaxis().tick_bottom() # ticks only on bottom axis
                    ax2.get_yaxis().tick_left() # ticks only on left axis
            start=start+nres
        #yticks=axt[0].yaxis.get_majorticklocs() # the y axis is shared
        #yticks=map(float,yticks)
        #yticks.remove(min(yticks))
        #axt[0].set_yticks(yticks)
        axt[0].set_ylim(ymin, ymax )
        if xlabel!=None : fig.text(0.5, 0.03, xlabel,fontsize=text_sizes['xlabel'], ha='center', va='center')
        if ylabel!=None :fig.text(0.015, 0.5, ylabel,fontsize=text_sizes['ylabel'],rotation='vertical', ha='center', va='center')
        if title!=None : fig.text(0.5, 0.97, title,horizontalalignment='center',fontsize=text_sizes['title']) 
        
        if label!='' and label!=None :
            if type(label) is not list and type(label) is not tuple : label=[label]
            legend=fig.legend(line_styles, label, loc=legend_location,prop={'size':legend_size},frameon=True,framealpha=0.5)
            legendframe = legend.get_frame()
            legendframe.set_facecolor((1.,1.,1.,0.7))
        if do_tight : fig.tight_layout(pad=3.5,h_pad=1.08,w_pad=1.08,rect=(0, 0, 1, 1))
        #if default_figure_sizes['all_tight'] : figure.tight_layout()
    else :
        raise Exception("Not implemented, set use_subplot to True")
    plt.draw() 
    if save!=None and save!='' :
        if '.' not in save : save+='.pdf'
        fig.savefig(save, dpi=default_figure_sizes['dpi'],bbox_inches="tight",transparent=True) #  bbox_inches=0 remove white space around the figure.. ,
    if show :
        plt.show(block=False)
    return fig,axt


def add_custom_legend(labels,facecolors,edgecolors=None,marker_types=None,markersize=None,linewidth=None,figure=None,proxy_point=1,frame=True, legend_location='upper right', legend_size=None,frame_facecolor=None,shadow=False,framealpha=None):
    '''
    this will draw a custom legend to the existing figure (or to figure if given)
    if you want to represent a line for label j give facecolor[j]=None and the desired edgecolor,
    if you wish to represent a marker give marker_types[j] != None
    frame_facecolor=(1.,1.,1.,0.7) give alpha of 0.7 to a white background
    '''
    if legend_size==None : legend_size=text_sizes['legend_size']
    proxy_point=int(proxy_point)
    if type(labels) is not list : labels=[labels]
    if type(facecolors) is not list : facecolors=[facecolors]*len(labels)
    if type(edgecolors) is not list :
        if  edgecolors==None : edgecolors=[None]*len(labels)
        else : edgecolors=[edgecolors]*len(labels)
    if type(marker_types) is not list : 
        if  marker_types==None : marker_types=[None]*len(labels)
        else : marker_types=[marker_types]*len(labels)
    if type(markersize) is not list :
        markersize=[markersize]*len(labels)
    if type(linewidth) is not list :
        if linewidth==None : linewidth=[None]*len(labels)
        else : linewidth=[linewidth]*len(labels)
        
    proxy_artists=[]
    for j in xrange(len(labels)) :
        if marker_types[j]!=None :
            pro = plt.Line2D(range(proxy_point), range(proxy_point), color="white", marker=marker_types[j],markersize=markersize[j], markerfacecolor=facecolors[j],markeredgecolor=edgecolors[j],linewidth=linewidth[j])
        elif facecolors[j]==None :
            pro = plt.hlines(y=proxy_point, xmin=proxy_point, xmax=proxy_point, color=edgecolors[j],linewidth=linewidth[j])
        else :
            pro = plt.Rectangle((proxy_point, proxy_point), 0, 0, facecolor=facecolors[j],edgecolor=edgecolors[j],linewidth=linewidth[j]) #,linewidth=2)
        proxy_artists+=[pro]
    if figure!=None :
        if type(legend_location) is tuple :leg=figure.legend( proxy_artists , labels,frameon=frame,numpoints=1, bbox_to_anchor=legend_location,prop={'size':legend_size}, shadow=shadow,framealpha=framealpha) 
        else :leg=figure.legend( proxy_artists , labels,frameon=frame,numpoints=1, loc=legend_location, prop={'size':legend_size}, shadow=shadow,framealpha=framealpha)
    else :
        if type(legend_location) is tuple :leg=plt.legend( proxy_artists, labels,frameon=frame,numpoints=1, bbox_to_anchor=legend_location, prop={'size':legend_size}, shadow=shadow,framealpha=framealpha) 
        else :leg=plt.legend( proxy_artists , labels,frameon=frame,numpoints=1, loc=legend_location,prop={'size':legend_size}, shadow=shadow,framealpha=framealpha)
    if frame_facecolor!=None :
        #leg=fig.legend(line_styles, label, loc=legend_location,prop={'size':legend_size},frameon=True,framealpha=0.5)
        legendframe = leg.get_frame()
        legendframe.set_facecolor(frame_facecolor)
    plt.draw()
    return leg
















    
        

def s2D_profile_to_color(output_matrices,color_rgb_list=[(0,0,1),(0,1,0),(1,0,0),(1.0, 1.0, 0.0)],add_PPII_to_coil=False):
    try :
        import chroma
    except ImportError :
        class chroma :
            def Color(self,*args,**kwargs):
                raise Exception("chroma MODULE NOT AVAILABLE. Cannot run s2D_profile_to_color\n in terminal try running 'pip install chroma'\n")
    if hasattr(output_matrices, 'shape') : # only one sequence
        output_matrices=[output_matrices] # mimic output from file with more sequences
    elif type(output_matrices) is str :
        seq_names, sequences, output_matrices, annotation_str =read_output(output_matrices, add_PPII_to_coil=add_PPII_to_coil)
    
    starting_cols=[ chroma.Color(c,format='RGB') for c in color_rgb_list]
    #print starting_cols
    out_color_list=[]
    for i,mat in enumerate(output_matrices) :
        out_color_list+=[[]]
        for res_out in mat :
            rc=None
            for j,x in enumerate(res_out) :
                nc=list(starting_cols[j].hsv)
                nc[1]=x
                #print nc,starting_cols[j].hsv,res_out
                if not hasattr(rc, 'rgb') : rc=chroma.Color( nc, format='HSV')
                else : rc -= chroma.Color( nc, format='HSV')
            out_color_list[-1]+=[rc.rgb]
    return out_color_list










def read_output(filename,force_float=True, add_PPII_to_coil=False, verbose=True):
    '''
    reads an output file of the s2D class
    return seq_names,sequences,output_matrices,annotation_str
      annotation_str will contain the content of eventual str columns (such the ss kind for s2D)
      output_matrices is a list of numpy 3D arrays with all the float output, one column per neuron (secondary structure type). 
      one array per sequence.
    '''
    #content=file(filename).read().splitlines()
    sequences=[]
    seq_names=[]
    annotation_str=[]
    output_matrices=[  ]
    for line in open(filename) :
        if len(line)<1 : continue
        if line[0]=='>' :
            seq_names.append( line[1:].strip() )
            sequences.append('')
            annotation_str.append('')
            output_matrices.append( [] )
        elif line[0]!='#' :
            if sequences==[] : # probably single-sequence output file
                sequences.append('')
                annotation_str.append('')
                output_matrices.append( [] )
                if verbose : print "s2D_class.read_output() -> Reading %s as a single-sequence output file!"
            line=line.split()
            if line!=[] :
                sequences[-1]+=line[1]
                output=[]
                ann=0
                for el in  line[2:] :
                    el,isnumber = convert_to_number(el, force_float=force_float)
                    if isnumber :
                        output.append(el)
                    elif ann==0 :
                        annotation_str[-1]+=el
                        ann+=1
                    else :
                        print '**W** ii s2D_class.read_output() Too many annotation strings in  line %s' % (str(line))
                output_matrices[-1].append(output)
    for j,out in enumerate(output_matrices) :
        if add_PPII_to_coil and len(out[0])>3 : output_matrices[j]=numpy.hstack(  ( numpy.array(out)[:,:2], numpy.array(out)[:,2:].sum(axis=1)[:,numpy.newaxis ] ) )
        else : output_matrices[j]=numpy.array(out)
    return seq_names,sequences,output_matrices,annotation_str

def average_multinetwork_prediction(list_of_network, seq_input_vec,postprocess=None):
    '''
    designed to average the results of multiple Seq network in case
    different sliding window sizes are employed (or any other difference in the network parameters)
    ''' 
    n=float(len(list_of_network))
    Pred=list_of_network[0].predict(seq_input_vec)
    for net in list_of_network[1:] :
        Pred += net.predict(seq_input_vec) # sum and average
    Pred/=n
    if postprocess!=None and hasattr(postprocess, '__call__') :
        Pred=postprocess(Pred)
    return Pred


def run_net(net,input_vec, postprocess=None):
    #print input_vec.shape
    if len(input_vec)==0 :
        raise Exception("**ERROR** can't perform prediction, EMPTY input vector")
    Pred=net.predict(input_vec)
    if postprocess!=None  :
        Pred=postprocess(Pred)
    return Pred







"""
BLAST FUNCTIONS
"""
chkparseC="""
/* On the first usage this gets printed in a .c file which is authomatically compiled */
/* content of the chkparse C script used to make the checkpoint file from psiblast more user friendly */
/* chkparse - generate PSIPRED compatible mtx file from BLAST+ checkpoint file */
/* V0.3 */
/* Copyright (C) 2010 D.T. Jones */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#define MAXSEQLEN 65536
#define EPSILON 1e-6
#define FALSE 0
#define TRUE 1
#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
const char *ncbicodes = "*A*CDEFGHIKLMNPQRSTVWXY*****";
/*  BLOSUM 62 */
const short           aamat[23][23] =
{
    {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0},
    {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1},
    {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1},
    {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1},
    {0, -3, -3, -3,10, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2},
    {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1},
    {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1},
    {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1},
    {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1},
    {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1},
    {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1},
    {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1},
    {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1},
    {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1},
    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2},
    {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0},
    {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0},
    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2},
    {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1},
    {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1},
    {-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1},
    {-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1},
    {0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, 4}
};
/* Standard BLAST+ a.a. frequencies */
float aafreq[26] =
{
    0.00000, 0.07805, 0.00000, 0.01925, 0.05364, 0.06295, 0.03856, 0.07377, 0.02199, 0.05142, 0.05744, 0.09019,
    0.02243, 0.04487, 0.05203, 0.04264,    0.05129, 0.07120, 0.05841, 0.06441, 0.01330, 0.00000, 0.03216, 0.00000,
    0.00000, 0.00000
};
/* PSSM arrays */
float fratio[MAXSEQLEN][28], pssm[MAXSEQLEN][28];
/* Dump a rude message to standard error and exit */
void
  fail(char *errstr)
{
    fprintf(stderr, "\\n*** %s\\n\\n", errstr);
    exit(-1);
}
/* Convert AA letter to numeric code (0-22 in 3-letter code order) */
int aanum(int ch)
{
    static const int aacvs[] =
    {
    999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
    22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };
    return (isalpha(ch) ? aacvs[ch & 31] : 22);
}
/* Scan ahead for file tokens */
void findtoken(char *buf, char *token, FILE *ifp)
{
    for (;;)
    {
    if (fscanf(ifp, "%s", buf) != 1)
        fail("Cannot find token in checkpoint file!");
    if (!token[0] || !strcmp(buf, token))
        break;
    }
}
/* Read hex sequence string */
int readhex(char *seq, FILE *ifp)
{
    int ch, aa, nres=0;
    while ((ch = fgetc(ifp)) != EOF)
    if (ch == '\\'')
        break;
    if (ch == EOF)
    fail("Bad sequence record in checkpoint file!");
    for (;;)
    {
    ch = fgetc(ifp);
    if (ch == '\\'')
        break;
    if (isspace(ch))
        continue;
    if (!isxdigit(ch))
        fail("Bad sequence record in checkpoint file!");
    if (ch >= 'A')
        aa = 16 * (10 + ch - 'A');
    else
        aa = 16 * (ch - '0');
    ch = fgetc(ifp);
    if (!isxdigit(ch))
        fail("Bad sequence record in checkpoint file!");
    if (ch >= 'A')
        aa += 10 + ch - 'A';
    else
        aa += ch - '0';
    if (nres > MAXSEQLEN)
        break;
    seq[nres++] = aa;
    }
    return nres;
}
/* This routine will extract PSSM data from a BLAST+ checkpoint file */
int getpssm(char *dseq, FILE *ifp)
{
    int i, j, len;
    float pssmrow[28], val, base, power;
    char buf[4096];
    findtoken(buf, "", ifp);
    if (strcmp(buf, "PssmWithParameters"))
    fail("Unknown checkpoint file format!");
    findtoken(buf, "numColumns", ifp);
    if (fscanf(ifp, "%d", &len) != 1)
    fail("Unknown checkpoint file format!");
    findtoken(buf, "ncbistdaa", ifp);
    if (len != readhex(dseq, ifp))
    fail("Mismatching sequence length in checkpoint file!");
    findtoken(buf, "freqRatios", ifp);
    findtoken(buf, "", ifp);
    for (i=0; i<len; i++)
    for (j=0; j<28; j++)
    {
        findtoken(buf, "", ifp);
        findtoken(buf, "", ifp);
        if (sscanf(buf, "%f", &val) != 1)
        fail("Unknown checkpoint file format!");
        findtoken(buf, "", ifp);
        if (sscanf(buf, "%f", &base) != 1)
        fail("Unknown checkpoint file format!");
        findtoken(buf, "", ifp);
        if (sscanf(buf, "%f", &power) != 1)
        fail("Unknown checkpoint file format!");
        findtoken(buf, "", ifp);

        fratio[i][j] = val * pow(base, power);
    }
    findtoken(buf, "scores", ifp);
    findtoken(buf, "", ifp);
    for (i=0; i<len; i++)
    for (j=0; j<28; j++)
    {
        findtoken(buf, "", ifp);
        if (sscanf(buf, "%f", &val) != 1)
        fail("Unknown checkpoint file format!");
        pssm[i][j] = val;
    }
    return len;
}
int roundint(double x)
{
    x += (x >= 0.0 ? 0.5 : -0.5);

    return (int)x;
}
int main(int argc, char **argv)
{
    int i, j, seqlen=0, nf;
    char seq[MAXSEQLEN];
    double scale, x, y, sxx, sxy;
    FILE *ifp;
    int use_psipred_format=0;
    if (argc != 2)
    fail("Usage: chkparse chk-file");
    ifp = fopen(argv[1], "r");
    if (!ifp)
    fail("Unable to open checkpoint file!");
    seqlen = getpssm(seq, ifp); // read the sequence from input file and save its length
    if (seqlen < 5 || seqlen >= MAXSEQLEN)
    fail("Sequence length error!");
    /* Estimate original scaling factor by weighted least squares regression */
    for (sxx=sxy=i=0; i<seqlen; i++)
    for (j=0; j<26; j++)
        if (fratio[i][j] > EPSILON && aafreq[j] > EPSILON)
        {
        x = log(fratio[i][j] / aafreq[j]);
        y = pssm[i][j];
        sxx += (y*y) * x * x; /* Weight by y^2 */
        sxy += (y*y) * x * y;
        }
    scale = 100.0 * sxy / sxx;
    if(use_psipred_format) 
    {
       printf("%d\\n", seqlen); // print sequence length
       for (i=0; i<seqlen; i++)  // print actual sequence
          putchar(ncbicodes[seq[i]]);
       printf("\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n");
       for (i=0; i<seqlen; i++)
       {
         for (j=0; j<28; j++)
               if (ncbicodes[j] != '*')
               {
             if (fratio[i][j] > EPSILON)
                   printf("%d  ", roundint(scale * log(fratio[i][j] / aafreq[j])));
             else
                      printf("%d  ", 100*aamat[aanum(ncbicodes[seq[i]])][aanum(ncbicodes[j])]);
             }
               else
               printf("-32768  ");
         putchar('\\n');
        }    
    }else
    {
      //print header
      printf("        ");
      for (j=0; j<28; j++)
          if (ncbicodes[j] != '*')
              printf("   %c   ",ncbicodes[j]);
      putchar('\\n'); 
  
      for (i=0; i<seqlen; i++)
      {
         printf("%5d %c ",i+1,ncbicodes[seq[i]]);
         for (j=0; j<28; j++)
            if (ncbicodes[j] != '*')
            {
             if (fratio[i][j] > EPSILON)
                   printf(" %5.2lf ", roundint(scale * log(fratio[i][j] / aafreq[j]))*0.01 );
             else
                      printf(" %5.2lf ", 1.*aamat[aanum(ncbicodes[seq[i]])][aanum(ncbicodes[j])]);
            }
       putchar('\\n');
       }
    }
    return 0;
}
"""

def psiblast_checkpoint(sequence, psi_blast_database_no_last_extension, BLAST_PATH='',sequence_name='to_blast',BLASTFILE=None,parser_executable=default_parser_executable, c_compiler=default_c_compiler,temporary_file_directory='',num_iterations=3,ncpu=2,str_content_of_chkparse=chkparseC,psi_blast_file=None):
    '''
    run psiblast but reads results from a checkpoint file so that the scoring_matrix results are not rounded down. 
    '''
    if BLASTFILE==None : BLASTFILE=temporary_file_directory+'blast_tmp.txt'
    database=psi_blast_database_no_last_extension
    if not os.path.isfile(parser_executable) :
        if os.path.isfile('chkparse.c') :
            os.system(c_compiler+' -O chkparse.c -lm -o '+parser_executable)
        elif str_content_of_chkparse!=None :
            out=open('chkparse.c','w')
            out.write(str_content_of_chkparse)
            out.close()
            os.system(c_compiler+' -O chkparse.c -lm -o '+parser_executable)
        else :
            raise IOError('***ERROR*** in psiblast_checkpoint() cant find chkparse (nor chkparse.c) in current folder and neither can find %s (its a compiled c file and its required)\n' % (parser_executable) )
    if BLAST_PATH!='' and not os.path.isfile(BLAST_PATH+'makeblastdb') :    #check if the BLAST_PATH is correct
        raise IOError('***ERROR*** in psiblast_checkpoint() path %s doesnt lead to blast directory where makeblastdb should be located' % (BLAST_PATH) )
    if psi_blast_file==None :
        psi_blast_file='psiblast_'+sequence_name.replace('|','').replace(' ','_').replace('/','').replace(':','_')+'.txt'
    pid=str(os.getpid())
    seq_file = file(temporary_file_directory+'sequence_tmp'+pid+'.fa','w')
    seq_file.write('> '+sequence_name+'\n'+sequence+'\n')
    seq_file.close()
    #check if the blast database has already been built, if not build it
    if (not ( os.path.isfile(database+'.phr') and os.path.isfile(database+'.pin') and os.path.isfile(database+'.psq')) and (not os.path.isfile(database+'.pal'))  ) : # .pal is for large datasets
        sys.stderr.write('\n********  WARNING  ********\n==> the blast database provided (%s) has not yet been compiled by makeblastdb\n    Running makeblastdb... (this will take a VERY LONG TIME, but needs to be done only once unless the database is deleted/renamed).\n    this may also print some ERROR messages "Error: (1431.1) FASTA-Reader:..." which can safely be ignored\n********           ********\n' % (database))
        sys.stderr.flush()
        try :
            os.system(BLAST_PATH+'makeblastdb -dbtype prot -in '+database)
        except :
            print '***ERROR*** in psiblast_checkpoint() cannot build blast database %s maybe you whish to  set BLAST_PATH to correct directory' % (database)
            raise
    # psipred command: "psiblast -db $dbname -query $tmproot.fasta -inclusion_ethresh 0.001 -out_pssm $tmproot.chk -num_iterations 3 -num_alignments 0 >& $tmproot.blast"
    if WEBSERVER_MODE : add=' 2> /dev/null'
    else : add=''
    os.system(BLAST_PATH+'psiblast -query '+temporary_file_directory+'sequence_tmp'+pid+'.fa'+' -db '+database+' -num_iterations '+str(num_iterations)+'  -inclusion_ethresh 0.001 -out_pssm '+temporary_file_directory+'tmproot'+pid+'.chk -num_alignments 0 -num_threads '+str(ncpu)+' > '+BLASTFILE+add)
    
    os.system(parser_executable+' '+temporary_file_directory+'tmproot'+pid+'.chk > '+psi_blast_file)
    os.system('rm -f '+temporary_file_directory+'tmproot'+pid+'.chk '+temporary_file_directory+'sequence_temp'+pid+'.fa')
    return psi_blast_file

def parse_psiblast_checkpoint_output(psi_blast_filename):
    '''
    # parse files generated by psi_blast() with psipred like options
    # it returns  entries,aa_header
    #  entries=[] # list of dictionary, one per position in the sequence, keys are 'id','aa', and 'scoring_matrix'
    #  aa_header list of the 20 residues in the order they appear in entry (both in 'occurrence' and in 'scoring_matrix')
    '''
    full_file=file(psi_blast_filename).read().splitlines()
    entries=[] # list of dictionary, one per position in the sequence
    aa_header=[]
    try :
        for i,line in enumerate(full_file) :
            #line=line.split()
            if line!='' and len(line)>6 :
                
                if line[:6].strip().isdigit() :
                    line=line.split()
                    entries+=[ {} ]
                    entries[-1]['id']=int(line[0])
                    entries[-1]['aa']=line[1]
                    entries[-1]['scoring_matrix']=map(float,line[2:])
                else :
                    line=line.split()
                    if line[0]=='A' :  # read the header
                        aa_header=line
    except Exception :
        sys.stderr.write("**ERROR** in parse_psiblast_checkpoint_output() while parsing %s at line %d\n" % (psi_blast_filename,i))
        raise
    return entries,aa_header


# just returns the psi_blast_file for analyisis (it does not parse it!)
# standard database is composed of three files Fragment.fasta.pin Fragment.fasta.psq Fragment.fasta.phr, you should only input Fragment.fasta. If the three files are not
#   found the database is generated
def psi_blast(sequence, psi_blast_database_no_last_extension, BLAST_PATH='/usr/local/ncbi/blast/bin/',sequence_name='to_blast',BLASTFILE='blast_tmp.txt',num_iterations=3,ncpu=1,psi_blast_file=None):
    database=psi_blast_database_no_last_extension
    if BLAST_PATH!='' and not os.path.isfile(BLAST_PATH+'makeblastdb') :    #check if the BLAST_PATH is correct
        raise IOError('***ERROR*** in psi_blast() path %s doesnt lead to blast directory where makeblastdb should be located' % (BLAST_PATH) )
    if psi_blast_file==None :
        psi_blast_file='psiblast_'+sequence_name.replace('|','').replace(' ','_').replace('/','').replace(':','_')+'.txt'
    
    seq_file = file('sequence_temp.txt','w')
    seq_file.write('> '+sequence_name+'\n'+sequence+'\n')
    seq_file.close()
    #check if the blast database has already been built, if not build it
    if (not ( os.path.isfile(database+'.phr') and os.path.isfile(database+'.pin') and os.path.isfile(database+'.psq')) and (not os.path.isfile(database+'.pal'))  ) : # .pal is for large datasets
        try :
            sys.stderr.write('\n==> the blast database provided (%s) has not yet been compiled by makeblastdb\n    Running makeblastdb... (this will take a LONG TIME, but needs to be done only once unless the database is deleted/renamed).\n' % (database))
            sys.stderr.flush()
            os.system(BLAST_PATH+'makeblastdb -dbtype prot -in '+database)
        except :
            sys.stderr.write( '\n\n***ERROR*** in psi_blast() cannot build blast database %s maybe you whish to  set BLAST_PATH to correct directory\n' % (database))
            raise
    if WEBSERVER_MODE : add=' 2> /dev/null'
    else : add=''
    os.system(BLAST_PATH+'psiblast -query sequence_temp.txt -db '+database+' -num_iterations '+str(num_iterations)+' -out_ascii_pssm '+psi_blast_file+' -out '+BLASTFILE+' -num_threads '+str(ncpu)+add)
    os.system('rm -f sequence_temp.txt')
    return psi_blast_file
# parse files generated by psi_blast() with -out_ascii_pssm option.
# it returns  entries,aa_header,general_results
#  entries=[] # list of dictionary, one per position in the sequence
#  aa_header list of the 20 residues in the order they appear in entry (both in 'occurrence' and in 'scoring_matrix')
#  general_results contains general results on the alignment such as 'Standard_UngappedK/L'
def parse_psiblast_output(psi_blast_filename, percentage_to_fractions=False, calculate_entropy_of_profile=True):
    '''
    # parse files generated by psi_blast() with -out_ascii_pssm option.
    # it returns  entries,aa_header,general_results
    #  entries=[] # list of dictionary, one per position in the sequence, keys are 'id','aa', 'occurrence' and 'scoring_matrix'
    #  aa_header list of the 20 residues in the order they appear in entry (both in 'occurrence' and in 'scoring_matrix')
    #  general_results contains general results on the alignment such as 'Standard_UngappedK/L'
    #   if calculate_entropy_of_profile then it contains also the entropy exp(sum(pi log(pi))) where pi is the number of different amino acids that appear in the alignment
    '''
    full_file=file(psi_blast_filename).read().splitlines()
    entries=[] # list of dictionary, one per position in the sequence
    aa_header=[]
    general_results={}
    try :
        if calculate_entropy_of_profile :
            entropy=0.
        for i,line in enumerate(full_file) :
            #line=line.split()
            if line!='' and len(line)>6 :
                
                if line[:6].strip().isdigit() :
                    if 161<=len(line)<=165 : # psiblat version quite old, 165 shoudl not be the case, it happens sometimes when the occurrences are all zeros
                        matrix_st, matrix_end,matrix_every,occ_start,occ_end = 9,69,3,70,150
                    elif 181<=len(line)<=185 : # recent psiblat version 
                        matrix_st, matrix_end,matrix_every,occ_start,occ_end = 9,89,4,90,170
                    else :
                        raise Exception("Line length of %d not recognized" % (len(line)))
                    entries+=[ {} ]
                    entries[-1]['id']=int(line[:6].strip())
                    entries[-1]['aa']=line[6]
                    entries[-1]['scoring_matrix']=[int(score) for score in split_every(line[matrix_st:matrix_end],matrix_every) ]
                    if percentage_to_fractions : entries[-1]['occurrence']=[float(score)/100. for score in line[occ_start:occ_end].split() ]
                    else : entries[-1]['occurrence']=[int(score) for score in line[occ_start:occ_end].split() ]
                    if calculate_entropy_of_profile :
                        if percentage_to_fractions : tmp=numpy.array( entries[-1]['occurrence'] )
                        else : tmp=numpy.array( entries[-1]['occurrence'] )/100.
                        
                        if all(tmp<=0.000001) : 
                            entropy += 1.
                        else : 
                            C_i = sum(tmp*numpy.log(tmp+0.00000001))# sum_over_20_aa( frequency * log(frequency))
                            entropy += numpy.exp( - C_i )  # exp( -sum )
                        
                        
                    entries[-1]['information']=float(line[occ_end:occ_end+7])
                    entries[-1]['gapless_match_to_pseudocount']=float(line[occ_end+7:])
                else :
                    line=line.split()
                    if line[0]=='A' :  # read the header
                        aa_header=line[:20]
                    elif 'Standard'==line[0] and 'Ungapped'==line[1] :
                        general_results['Standard_UngappedK/L']=[float(l) for l in line[2:]]
                    elif 'Standard'==line[0] and 'Gapped'==line[1] :
                        general_results['Standard_GappedK/L']=[float(l) for l in line[2:]]
                    elif 'PSI'==line[0] and 'Ungapped'==line[1] :
                        general_results['PSI_UngappedK/L']=[float(l) for l in line[2:]]
                    elif 'PSI'==line[0] and 'Gapped'==line[1] :
                        general_results['PSI_GappedK/L']=[float(l) for l in line[2:]]
        if calculate_entropy_of_profile :
            entropy/=(1.*len(entries))
            general_results['entropy']=entropy
    except Exception :
        sys.stderr.write("\n**ERROR** in parse_psiblast_output() while parsing %s at line %d\n" % (psi_blast_filename,i))
        raise
    return entries,aa_header,general_results



def split_every(staff, num):
    '''
    split a string a list or a tuple every num elements
    '''
    return [ staff[start:start+num] for start in range(0, len(staff), num) ]

def loose_compare_keyword(keyword, listt , begin_with=False, end_with=False):
    '''
    see if keyword is contained in any of the elements in listt.
    Only the first element satisfying this condition is returned
    
     if begin_with keyword has to be found at the beginning of the element, 
     if end_with keyword has to be found at the end of the element.
     if both are True keyword has to be either at the end or at the beginning of the element
    it returns True and the first matching element or False and None
    '''
    le=len(keyword)
    if begin_with and end_with :
        for k in listt :
            if len(k) >= le :
                if keyword==k[:le] or keyword==k[-le:] :
                    return True,k
        return False,None
    if begin_with :
        for k in listt :
            if len(k)>= le :
                if keyword==k[:le] :
                    return True,k
        return False,None
    if end_with :
        for k in listt :
            if len(k)>= le :
                if keyword==k[-le:] :
                    return True,k
        return False,None
    for k in listt :
        if keyword in k : return True,k
    return False,None




class s2D_method :
    def __init__(self, parameter_class=parameter_filename ,network_files=None, deep_learning_net_file=None,Nto1_netfile=None, postprocess=normalize_numpy,Nto1Postprocess=None,multiply_before_DL=10., use_psiblast=True, net_trained_on_scoring_matrix=True 
                 ,psiblast_ncpu=2,uniref90_psi_blast_database=None,use_psiblast_checkpoint=True,temporary_file_directory='/tmp/',keep_pssm_profile_files=True,folder_with_psiblast_files='psiblast_files/' \
                 ,out_tags=['Helix','Beta','Coil','Polyproline-II']):
        
        if parameter_class!=None and network_files==None :
            if type(parameter_class) is str :
                tmp=parameter_class
                parameter_class=s2D.s2D_parameters()
                parameter_class.read(tmp)
            network_files=parameter_class.networks
            deep_learning_net_file=parameter_class.DL_network
            temporary_file_directory=parameter_class.temporary_file_directory
            folder_with_psiblast_files=parameter_class.psiblast_files_folder
            uniref90_psi_blast_database=parameter_class.psiblast_database
            use_psiblast_checkpoint=parameter_class.use_psiblast_checkpoint
            psiblast_ncpu=parameter_class.psiblast_ncpu
            keep_pssm_profile_files=parameter_class.keep_pssm_profile_files
        elif network_files==None  :
            sys.stderr.write("**WARNING** parameter file not found or not declared\n")
        self.network_files=network_files
        self.Nto1_net_list=None
        self.multiply_before_DL=multiply_before_DL 
        self.max_window_size=0
        if self.multiply_before_DL==None : self.multiply_before_DL=1.
        if type(self.network_files) is list or type(self.network_files) is tuple :
            self.net=[]
            self.Nto1_net_list=[]
            try :
                for j,net in enumerate(self.network_files) :
                    self.net += [ PyELM.loadModel(net,verbose=not WEBSERVER_MODE) ]
                    self.net[-1].name=net
                    if isinstance(self.net[-1],PyELM.WNto1) or isinstance(self.net[-1],PyELM.Nto1) : self.Nto1_net_list+=[j]
                    elif self.net[-1].win>self.max_window_size : self.max_window_size=self.net[-1].win
                if not WEBSERVER_MODE: print 'loaded %d networks' % (len(self.net))
                if self.Nto1_net_list!=[] :
                    remove=sorted(self.Nto1_net_list,reverse=True)
                    self.Nto1_net_list=[ self.net[j] for j in remove ]
                    for j in remove : del self.net[j]
                    if not WEBSERVER_MODE: print ' of which %d Nto1 networks: %s' % (len(self.Nto1_net_list),' '.join([n.name for n in self.Nto1_net_list]))
                    self.Nto1_net_list.sort(key=lambda x : x.name )
                else : 
                    self.Nto1_net_list=None
                self.net.sort(key=lambda x : x.name )
                if not WEBSERVER_MODE: print 'Sorted:',[n.name for n in self.net]
            except Exception :
                sys.stderr.write('\n**ERROR** probably when attempting to load network %s,\n   maybe the file is not in the correct directory\n\n' % (net))
                raise
        elif network_files!=None :
            self.net=PyELM.loadModel(self.network_files,verbose=not WEBSERVER_MODE)
        self.deep_learning_net_file=deep_learning_net_file
        if self.deep_learning_net_file!=None :
            self.DL=True
            self.DL_net=PyELM.loadModel(self.deep_learning_net_file,verbose=not WEBSERVER_MODE)
            self.DL_net.name=self.deep_learning_net_file
            if 'win' in dir(self.DL_net) and self.DL_net.win>self.max_window_size : self.max_window_size=self.DL_net.win
            elif 'window' in dir(self.DL_net) and self.DL_net.window>self.max_window_size : self.max_window_size=self.DL_net.window
        else :
            self.DL=False
        self.Nto1Postprocess=Nto1Postprocess
        self.Nto1_net=None
        self.Nto1_netfile=Nto1_netfile
        self.run_Nto1_first=False # this is an option never tested that allows to run an Nto1 network before anything else, and then use it as a deep learning
        if self.Nto1_net==None and self.Nto1_netfile!=None :
            self.Nto1_net=PyELM.loadModel(self.Nto1_netfile,verbose=not WEBSERVER_MODE)
            self.run_Nto1_first=True
        self.postprocess=postprocess
        self.out_tags=out_tags
        self.psiblast_ncpu=psiblast_ncpu
        self.uniref90_psi_blast_database=uniref90_psi_blast_database
        self.psi_blast_path=''
        self.use_psiblast_checkpoint=use_psiblast_checkpoint
        self.keep_pssm_profile_files=keep_pssm_profile_files
        if not os.path.isdir(temporary_file_directory) :
            sys.stderr.write("**Warning declared temporary_file_directory %s does not seem to exist, using the current directory.\n" % (str(temporary_file_directory)))
            temporary_file_directory=''
        elif temporary_file_directory[-1]!='/' : temporary_file_directory+='/'
        self.temporary_file_directory=temporary_file_directory
        self.folder_with_psiblast_files=folder_with_psiblast_files
        self.function_residue_to_input=residue_to_input 
        self.use_psiblast=use_psiblast
        self.seq_name=None
        self.net_trained_on_scoring_matrix=net_trained_on_scoring_matrix
    def run_on_psiblast_file(self, psiblast_file,input_sequence=None):
        inp,sequence=self.input_from_psiblast_file(psiblast_file,input_sequence=input_sequence, delete_pssm_file=False)
        if inp==None : return None,None,None
        if self.run_Nto1_first :
            Nto1_out,_=self.run_Nto1(sequence, read_staff_from_psiblastfile=None, seq_name=self.seq_name)
            numpy.vstack((Nto1_out,inp))
        if not self.DL : # just run the networks if no Deep learning is requested
            if type(self.net) is list :
                output=average_multinetwork_prediction(self.net, inp, postprocess=self.postprocess)
            else :
                output=run_net(self.net, inp, postprocess=self.postprocess)
        else :
            # run networks
            self.net_output=run_net(self.net[0], inp, postprocess=None) # do not postprocess before DL!
            for ne in self.net[1:] :
                self.net_output= numpy.hstack((self.net_output, run_net( ne, inp, postprocess=None) ) ) # do not postprocess before DL!
            self.dl_input=numpy.hstack((inp, self.multiply_before_DL* self.net_output))
            # run N-to-1 network
            if self.Nto1_net_list!=None :
                self.n1_out=run_net(self.Nto1_net_list[0],[self.dl_input[:,:self.Nto1_net_list[0].dimInput]], postprocess=self.Nto1Postprocess)[0] # run the Nto1, This won't have one vertical output per amino acid, but only one for the whole sequence
                for Nto1 in self.Nto1_net_list[1:] :
                    self.n1_out= numpy.hstack( (self.n1_out, run_net(Nto1,[self.dl_input[:,:Nto1.dimInput]], postprocess=self.Nto1Postprocess)[0] )) # run the Nto1
                self.dl_input=numpy.hstack((self.dl_input, self.multiply_before_DL*numpy.array(self.n1_out)*numpy.ones( (inp.shape[0],1) ) ))
            # run DL network
            output=run_net(self.DL_net, self.dl_input, postprocess=self.postprocess)
        return inp,sequence, output
    
    def input_from_psiblast_file(self, psiblast_file,input_sequence=None,delete_pssm_file=False):
        if self.use_psiblast_checkpoint :
            entries,_ = parse_psiblast_checkpoint_output(psiblast_file )
            #print 'entries',len(entries),len(entries[0]['scoring_matrix'])
        else :
            entries,_,_ = parse_psiblast_output(psiblast_file , percentage_to_fractions=True,calculate_entropy_of_profile=False)
        add_aa_specific_neurone=False # determine whether we should add a specific neurone or not
        if entries==[] :
            sys.stderr.write('\n***WARNING*** blast file %s is probably empty.. reblasting...\n' % (psiblast_file.split('/')[-1]))
            sys.stderr.flush()
            return None,None
        elif len(entries[0]['scoring_matrix'])==self.net[0].numID-1:
            add_aa_specific_neurone=True
        sequence=''
        if self.net_trained_on_scoring_matrix or self.use_psiblast_checkpoint :
            inp=numpy.array(entries[0]['scoring_matrix'])
            sequence+=entries[0]['aa']
            for en in entries[1:] :
                inp=numpy.vstack((inp, numpy.array(en['scoring_matrix'])))
                sequence+=en['aa']
            #print 'inp',inp.shape
        else :
            inp=numpy.array(entries[0]['occurrence'])
            sequence+=entries[0]['aa']
            for en in entries[1:] :
                inp=numpy.vstack((inp, numpy.array(en['occurrence'])))
                sequence+=en['aa']
        if input_sequence!=None and input_sequence!=sequence :
            sys.stderr.write('*WARNING* sequence read from blast file %s is different from the one given as input.. reblasting...' % (psiblast_file))
            sys.stderr.flush()
            if delete_pssm_file :os.system('rm -f %s' % (psiblast_file))
            return None,None
        if add_aa_specific_neurone :
            sp=[]
            for s in sequence :
                sp.append(residue_to_input(s, numb_neurones=1,out_in_range=(-10.,10.), use_rank=res_closness_blosum62,nres=20))
            #print len(sp),inp.shape
            inp=numpy.hstack( (numpy.array(sp,ndmin=2).T,inp))
        if delete_pssm_file :os.system('rm -f %s' % (psiblast_file))
        return inp,sequence
    def psiblast(self,sequence):
        '''
        run psiblast on the sequence using all the default parameters. It returns the path to the output psiblast file.
        '''
        psiblast_file,_=self.psiblast_sequence(sequence, sequence_name=self.seq_name, uniref90_psi_blast_database=self.uniref90_psi_blast_database, psi_blast_path=self.psi_blast_path, folder_with_psiblast_files=self.folder_with_psiblast_files, keep_psiblast_alignments=False)
        return psiblast_file
    def psiblast_sequence(self,sequence,sequence_name='myseq', uniref90_psi_blast_database=None, psi_blast_path='', folder_with_psiblast_files='psiblast_files/', keep_psiblast_alignments=False):
        '''
        used internally, use psiblast() as a user
        This run psiblast on the sequenc unless the corresponding psiblast output file is found in the folder folder_with_psiblast_files
        return psiblast_file,gzip (the path to the file (with file name) and gzip=True if the file was found gzipped (however the actual file has been gunzipped by the function when the function returns).
        psi_blast_path=='' means psiblast is added to the system path and can be called from anywhere. (make sure the path ends with '/' if psiblast is not installed globally
        '''
        if type(sequence) is not str : # we assume is a Seq object
            sequence_name=sequence.id
            sequence=str(sequence.seq)
        
        if uniref90_psi_blast_database==None or not os.path.isfile(uniref90_psi_blast_database) :
            sys.stderr.write("\n**ERROR** psiblast database not found or not specified in run_on_sequence()!!\n  (default should be uniref90 filtered from low-complexity regions) declared %s \n" % (str(uniref90_psi_blast_database)) )
            sys.stderr.flush()
            return 1
        if folder_with_psiblast_files[-1]!='/' : folder_with_psiblast_files+='/'
        if not os.path.isdir(folder_with_psiblast_files) :
            os.system('mkdir '+folder_with_psiblast_files)
        
        psi_blast_file=sequence_name.replace('|','').replace(' ','_').replace('/','').replace(':','_')+'_psi_blast.txt'        
        gzip=False
        blast_anyway=True

        #print folder_with_psiblast_files+psi_blast_file
        # check if the psiblast file for this sequence name already exists (it is kind of useless to have it saved according to the sequence name)
        if os.path.isfile(folder_with_psiblast_files+psi_blast_file+'.gz') :
            psiblast_file=folder_with_psiblast_files+psi_blast_file
            gzip=True
            os.system('gunzip '+folder_with_psiblast_files+psi_blast_file+'.gz')
            psiblast_file=folder_with_psiblast_files+psi_blast_file
            _, sequence2 = self.input_from_psiblast_file(psiblast_file, input_sequence=None, delete_pssm_file=False)
            if sequence==sequence2 : blast_anyway=False
        elif os.path.isfile(folder_with_psiblast_files+psi_blast_file) :
            psiblast_file=folder_with_psiblast_files+psi_blast_file
            _, sequence2 = self.input_from_psiblast_file(psiblast_file, input_sequence=None, delete_pssm_file=False)
            if sequence==sequence2 : blast_anyway=False 
        if  blast_anyway :
            psiblast_alignments=sequence_name.replace('|','').replace(' ','_').replace('/','').replace(':','_')+'_blast.txt'
            if self.use_psiblast_checkpoint : psiblast_file= psiblast_checkpoint(sequence, uniref90_psi_blast_database, BLAST_PATH=psi_blast_path, sequence_name=sequence_name, BLASTFILE=psiblast_alignments,temporary_file_directory=self.temporary_file_directory, num_iterations=3,ncpu=self.psiblast_ncpu, psi_blast_file=folder_with_psiblast_files+psi_blast_file)
            else : psiblast_file=psi_blast(sequence, uniref90_psi_blast_database, BLAST_PATH=psi_blast_path, sequence_name=sequence_name, BLASTFILE=psiblast_alignments, num_iterations=3,ncpu=self.psiblast_ncpu, psi_blast_file=folder_with_psiblast_files+psi_blast_file)
            if not keep_psiblast_alignments : 
                os.system('rm -f %s' % (psiblast_alignments)) # remove blast file with the alignment, we need only the psiblast one
        return psiblast_file,gzip
    def run_on_sequence(self,sequence, uniref90_psi_blast_database=None,sequence_name='myseq', psi_blast_path='', folder_with_psiblast_files='psiblast_files/', keep_psiblast_alignments=False):
        '''
        psi_blast_path=='' means psiblast is added to the system path and can be called from anywhere. (make sure the path ends with '/' if psiblast is not installed globally)
        '''
        if type(sequence) is not str : # we assume is a Seq object
            sequence_name=sequence.id
            sequence=str(sequence.seq)
        if self.use_psiblast :
            psiblast_file,gzip= self.psiblast_sequence(sequence,sequence_name=sequence_name, uniref90_psi_blast_database=uniref90_psi_blast_database, psi_blast_path=psi_blast_path, folder_with_psiblast_files=folder_with_psiblast_files, keep_psiblast_alignments=keep_psiblast_alignments)
            inp,sequence2, output = self.run_on_psiblast_file(psiblast_file,input_sequence=sequence)
            if gzip :
                os.system('gzip '+psiblast_file)
            if sequence!=sequence2 :
                if inp!=None : sys.stderr.write("**ERROR** after psiblast. Returned sequence is different from original one. Maybe in folder %s there is already a file named %s. Remove it or give a different sequence !!\n\n" % (folder_with_psiblast_files,psiblast_file))
                else : sys.stderr.write("psiblast failed on sequence:\n%s\n" % (sequence))
                return None , None
            if not self.keep_pssm_profile_files : os.system('rm -f %s*' % (psiblast_file))
        return inp, output
    def run_Nto1(self,sequence,read_staff_from_psiblastfile=None,seq_name='seq_for_s2D'):
        '''
        it runs an Nto1 network
        '''
        if type(sequence) is str :
            self.sequence=sequence
            if self.seq_name==None : self.seq_name=seq_name
        else :
            self.sequence=str(sequence.seq)
            self.seq_name=str(sequence.id)
        if read_staff_from_psiblastfile!=None :
            if read_staff_from_psiblastfile==True :
                read_staff_from_psiblastfile=self.psiblast(sequence)
            inp,_=self.input_from_psiblast_file(read_staff_from_psiblastfile, input_sequence=sequence, delete_pssm_file=False)
        else :
            inp=[]
            for s in sequence :
                inp.append(residue_to_input(s, numb_neurones=20, use_rank=res_closness_blosum62,nres=20))
        if self.Nto1_net==None and self.Nto1_netfile!=None :
            self.Nto1_net=PyELM.loadModel(self.Nto1_netfile,verbose=not WEBSERVER_MODE)
        Pred=run_net(self.Nto1_net, [inp], postprocess=self.postprocess)
        return Pred, inp
    def run(self,sequence, seq_name='seq_for_s2D',keep_psiblast_alignments=False):
        if type(sequence) is str :
            self.sequence=sequence
            self.seq_name=seq_name
        else :
            self.sequence=str(sequence.seq)
            sequence.id=sequence.id #.replace('|','').replace(' ','_').replace('/','').replace(':','_')
            self.seq_name=str(sequence.id)
        #if self.use_psiblast : print 'Using psiblastDB %s' % (str(self.uniref90_psi_blast_database))
###
        #print self.sequence,self.seq_name
        self.input, self.output= self.run_on_sequence(self.sequence,sequence_name=self.seq_name,uniref90_psi_blast_database=self.uniref90_psi_blast_database, folder_with_psiblast_files=self.folder_with_psiblast_files,keep_psiblast_alignments=keep_psiblast_alignments)
    def get_ss_kind(self,out_per_res):
        f=out_per_res.argmax()
        if f==0 : return 'H'
        if f==1 : return 'E'
        if f==2 : return 'C'
        if f==3 : return 'P'
    def print_results(self, out=sys.stdout):
        close_file=False
        if type(out) is str :
            close_file=True
            print '\nSaving output to file %s...' % (out)
            out=open(out,'w')
        if len(self.output)!=len(self.sequence) :
            sys.stderr.write("ERROR in print_results() len output!= len sequence %d != %d\n\n" % (len(self.output),len(self.sequence)))
        self.ss_string=''
        self.ss_profiles=[]
        n_out=len(self.output[0])
        out.write('> %s\n' % (self.seq_name))
        tmp='#rid\tresname'
        for i in xrange(n_out) :
            self.ss_profiles.append([]) 
            tmp+='\t%s' % (self.out_tags[i])
        tmp+='\tss_kind'
        out.write(tmp+'\n')
        for j,r in enumerate(self.sequence) :
            tmp=''
            for i,float_n in enumerate(self.output[j]) : 
                tmp+='\t%lf' % (float_n)
                self.ss_profiles[i]+=[ float_n ]
            k=self.get_ss_kind(self.output[j])
            out.write('%d\t%s%s\t%s\n' % (j+1,r, tmp  ,k))
            self.ss_string+=k
        if close_file :
            out.close()
    def plot_results(self,save=True,show=False,dpi=250, plotCoil=False,savefolder='',**kwargs):
        if not can_plot :
            raise Exception("When first importing the module there was a problem importing matplotlib (maybe you don't have version >= 1.4?)")
        if type(save) is bool and save==True :
            save =self.seq_name.replace('|','_').replace(' ','_').replace('/','_').replace(':','_')+'_s2D_plot.png'
        if savefolder!='' and savefolder!=None :
            if savefolder[-1]!='/' : savefolder+='/'
            save=savefolder+save
        plot_s2D_results(self.output,self.sequence,self.ss_string,seq_names=self.seq_name,bar=True,dont_plot_coil=not plotCoil,y_range=(0,1),dpi=dpi,save=save,show=show,**kwargs)

