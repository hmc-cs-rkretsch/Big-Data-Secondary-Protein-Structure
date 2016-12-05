"""
    Copyright (C) 2014 Pietro Sormanni (for s2D) and Piero Fariselli (for PyELM) (Also David T. Jones for pfilt and chkparse included in 'extra/')

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
        Please contact me if you have problems running or installing s2D. 
        If you have problems installing python, numpy or blast+ please 
        contact the relevant support.
              
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it
    under the condition of preserving this text

"""


import os,sys
import shutil
from distutils.sysconfig import get_python_lib
import sysconfig
import site
#site.USER_BASE

PARAMETER_FILE='s2D_parameters.txt'
CREATE_BIN=True
default_c_compiler = sysconfig.get_config_vars('CC')[0]
if default_c_compiler==None : default_c_compiler='cc'


param_file_cont=None
try :
    if not os.path.isfile(PARAMETER_FILE) :  
        raise Exception("\n**ERROR** parameter file not found. setup only works with the paths to the network files and the psiblast database in a file named %s" % (PARAMETER_FILE))

    from s2D import s2D_parameters
    parms=s2D_parameters()
    parms.read(PARAMETER_FILE)
    param_file_cont=file(PARAMETER_FILE).read()    
    #make paths absolute
    if len(sys.argv)>1 :
       if os.path.isdir(sys.argv[1]) : destination_folder=sys.argv[1]
       elif 'user' in sys.argv[1] or '--user' in sys.argv : destination_folder=site.getusersitepackages()
       elif 'install' in sys.argv[1] : destination_folder=get_python_lib()
    else : destination_folder=get_python_lib()
    if destination_folder[-1]!='/' : destination_folder+='/'
    print "\n==> Installing in folder '%s' \n" % (destination_folder)
    parms.psiblast_database=os.path.abspath(parms.psiblast_database)
##
    print parms.psiblast_files_folder
    if parms.psiblast_files_folder!='.' and  parms.psiblast_files_folder!='./' : parms.psiblast_files_folder=os.path.abspath(parms.psiblast_files_folder)
    print parms.psiblast_files_folder
    parms.temporary_file_directory=True # in this way each time the tempfile module will authomatically locate the current temporary directory
   
    # check whether there is the psiblast database
    if not os.path.isfile(parms.psiblast_database) :
        sys.stderr.write('\n*** WARNING*** blast database not found.\n* Cannot locate the psiblast database %s, declared in the s2D parameters file %s (with keyword psiblast_database).\n* We reccomend that you follow the instructions in README.txt to install BLAST+ and create the relevant psiblast database (Uniref90 filtered with pfilt).\n* Then you can update the keyword psiblast_database in the file %s so that it points to the database before running setup.py to install s2D.\n ' % (str(parms.psiblast_database),PARAMETER_FILE,PARAMETER_FILE))
        c=raw_input("If you already have the correct blast database on your computer then type the path to it (e.g. /Users/Shared/blastDB/uniref90filt.fasta), otherwise decide if you wish to continue installing s2D anyway\ntype [ path / y / n] and press return\n")
        while True :
            if len(c.strip())==1 : 
                if c.strip().lower()=='n' : sys.exit(0)
                elif c.strip().lower()=='y' : break
            elif os.path.isfile(c) : 
                parms.psiblast_database=os.path.abspath(c)
                break
            else : 
                sys.stderr.write("\n*** WARNING*** no file found at:\n'%s'\n" % (c))
                c= raw_input("TYPE: [ path / y / n] and press return\n  n= abort installation;\n  y= continue installation (it won't work without a blast database)\n  path= path to the blastdb file\n")
    print '\nUsing psiblast database %s\n' % (str(parms.psiblast_database))
    # update network folder name
    all_nets=parms.networks[:]+[ parms.DL_network ] #dataFiles= [ ('s2D_networks', parms.networks[:]+[ parms.DL_network ] ),('.',['s2D_parameters.txt']) ]
    parms.networks=[ destination_folder+'s2D_networks/'+n.split('/')[-1] for n in parms.networks ]
    parms.DL_network = destination_folder+'s2D_networks/'+parms.DL_network.split('/')[-1] 
    parms.write('s2D_parameters.txt')
    
    # INSTALL
    from s2D_class import chkparseC,default_parser_executable
    default_parser_executable=default_parser_executable.split('/')[-1]
    out=open('chkparse.c','w')
    out.write(chkparseC)
    out.close()
    # compile chkparse to parse the psiblast checkpoint file with PSSM matrix with float numbers (not rounded down).
    os.system(default_c_compiler+' -O chkparse.c -lm -o '+default_parser_executable)
    shutil.copy(default_parser_executable,destination_folder+default_parser_executable)
    shutil.copy('s2D_parameters.txt',destination_folder+'s2D_parameters.txt')
    if not os.path.exists(destination_folder+'s2D_networks'):
        os.makedirs(destination_folder+'s2D_networks')
    for n in all_nets :
        shutil.copy(n,destination_folder+'s2D_networks/.')
    shutil.copy('s2D.py',destination_folder+'s2D.py')
    shutil.copy('s2D_class.py',destination_folder+'s2D_class.py')
    shutil.copy('PyELM.py',destination_folder+'PyELM.py')


    
    sys.stdout.write(" Python module installation completed. Creating executable\n")
    sys.stdout.flush()
    # make executable (can be moved in bin directory)
    if CREATE_BIN :
        bin_str='echo "python %ss2D.py \$*" > s2D\nchmod 777 s2D'  % (destination_folder)
        os.system(bin_str)
        sys.stdout.write("\n====> Created an executable file named 's2D'.\nYou may copy this file into a bin directory included in your $PATH (e.g. '/usr/local/bin' or '/usr/bin' on Unix) to be able to call s2D from any folder\n\n") 
    #sh_str='if [ $# -eq 1 ]; then\n      PA=$1\nelse\n      PA="/usr/bin/"\nfi\necho "Creating an executable link to s2D in directory $PA" \necho "python %ss2D.py \$*" > s2D\nchmod 777 s2D\nmv s2D $PA/.\nexit\n' % (destination_folder)
    #out=open('install_s2D_bin.sh','w')
    #out.write(sh_str)
    #out.close()


except Exception :
  if param_file_cont!=None :
    out=open(PARAMETER_FILE,'w')
    out.write(param_file_cont)
    out.close()
  raise

if param_file_cont!=None :
    out=open(PARAMETER_FILE,'w')
    out.write(param_file_cont)
    out.close()
    try :
        os.system("rm -f *.pyc chkparse*")
    except Exception : pass
    print "DONE!  ==> s2D is now installed as a python module. See README.txt for usage.\n"




"""


  from distutils.core import setup, Extension



#module1 = Extension('s2D',define_macros = [('MAJOR_VERSION', '1'),('MINOR_VERSION', '0')],sources = ['s2D.py'])
#module2 = Extension('s2D_class', sources = ['s2D_class.py'])
#module3 = Extension('PyELM', sources = ['PyELM.py'])


  setup (name = 's2D',
       version = '2.0',
       description = 'This is a module to run s2D to predict protein secondary structure populations from the amino acid sequence',
       author = 'P. Sormanni',
       author_email = 'ps589@cam.ac.uk',
       url = 'http://docs.python.org/extending/building',
       long_description = '''
''',   py_modules=['s2D', 's2D_class','PyELM'],
       include_package_data = True,
       package_data = {'s2D' : dataFiles, 's2D_class':dataFiles },
       data_files=dataFiles ) # ,   ext_modules = [module1,module2,module3]

"""
