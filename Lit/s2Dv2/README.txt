    Copyright (C) 2014-2016 Pietro Sormanni (for s2D) and Piero Fariselli (for PyELM) 
    (Also David T. Jones for pfilt and chkparse included in 'extra/')

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
     If you have problems installing python, numpy or blast+ please contact the 
     relevant support.
              
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it
    under the condition of preserving this text





PLEASE REPORT BUGS/ISSUES TO: Pietro Sormanni ps589@cam.ac.uk


REQUIREMENTS

In order to function s2D needs the following software:

python 2.7   (should work with other python versions as well, 2.7 has been tested)
numpy        (generally included with python as a module)
psiblast     (distributed with BLAST+, tested versions are 2.2.28 and 2.2.29 but should work with newer versions also. See section "obtain blast+" below)

In addition it needs the filtered UniRef90 Database, against which the psiblast search is carried out (see section "the psiblast database" below).

the s2D_parameters.txt file must be edited to point to the location of this database (see “very important” below).


OBTAIN BLAST+

Blast+ executables and source code can be downloaded from the NCBI web site at http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download for most operating systems. To retrieve the installer login as guest to the ftp server ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ (just copy and paste url in your web browser) or to ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ to install the version tested with s2D (the latest should work anyway).
You should isntall psiblast so that it can be called from any folder in the terminal by simply typing 'psiblast' (this may require adding the blast+ directory to your path depending on your system or how you installed blast+)



THE PSIBLAST DATABASE

Two options are possible:
1- (preferred) The ready-to-be-compiled database can be downloaded together with the s2D source code at:
   http://www-mvsoftware.ch.cam.ac.uk/s2D (large file! > 3Gb). Unzip the file in the destination folder of your choice (see section 'very important' below)

2- Alternatively one can download the file uniref90.fasta.gz database from ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/ 
   In this case then the file pfilt.c in the folder extra/ needs to be compiled (e.g. give the command 'cc  -O pfilt.c -lm -o pfilt') and the database
   filtered to exclude low complexity regions. 'pfilt uniref90.fasta > uniref90filt.fasta' (you need to unzip first 'gunzip -v uniref90.fasta.gz')


Once the database is downloaded it should be placed in a folder and its path should be declared by editing the s2D_parameters.txt file




VERY IMPORTANT *** ***

One should tell s2D where the psiblast database is located. 
To this end edit the s2D_parameters.txt file and replace the path found in the line beginning with the keyword 'psiblast_database' (do not change the keyword!) with the path of of your choice. The path can be either absolute (preferred) or relative to the folder where s2D.py and setup.py are located. The edited s2D_parameters.txt line should look like 'psiblast_database:   /My/Path/uniref90filt.fasta' (without the apostrophes ').

If you plan to install s2D for all users of the computer, the database should be located in a folder accessible to every users (MacOSX default '/Users/Shared' Linux default depends on your configuration (maybe '/home/shared/') and Windows default should be 'Shared Documents/')
Allocating more cpu will improve speed (see below).



THE s2D_parameters.txt FILE 

this file contains some parameters that can be changed to customize the performance of the s2D method.
The parameters that can be customized are all given as 'keyword    parameter_value'. The most important one is the psiblast_database discussed above in section 'very important'. Other parameters that can be edited are:

ncpu:  2	# the number of cpus used by psiblast. In general the more you can allocate (within the limit or your computer) the faster the psiblast search will be (but check the blast+ documentation for the keyword num_threads therein)
keep_pssm_profile_files:  False  # If True the pssm files generated by psiblast are kept in the folder given by the keyword psiblast_files_folder (otherwise are deleted)
psiblast_files_folder: .	 # the path to a folder where the pssm profile files are stored if keep_pssm_profile_files is set to True (otherwise it is used as a temporary files folder for the pssm profiles) '.' means current directory (that will be the directory from which s2D is called).

RUNNING WHOLE PROTEOMES?
Note that the speed-limiting step of s2D is the psiblast search. The performance of psiblast does not increase linearly with the number of cpus used by it (see blast+ documentation). If your system has e.g. 4 physical cores corresponding to 8 virtual ones (hyperthreading enabled) it will be faster to split your input proteome file in 2 files with approximately the same number of sequences each and call s2D twice (once per file) allocating 4 cores per run. Doing one s2D call on the whole proteome with 8 allocated cores will be slower.



TEST s2D

once you have psiblast installed (typing 'psiblast' in the terminal should yield "BLAST query/options error: ... ")
and the psiblast database correctly declared in the file s2D_parameters.txt
you can test s2D by typing (in the current folder where s2D.py is located)

python s2D.py test_sequence.fasta -plotCoil

the option -plotCoil requires matplotlib version >= 1.4 to be installed, if you do not have it omit the option (you will get the text output and not the plot).
If you have never compiled the psiblast database (with makeblastdb) s2D will compile it authomatically, this takes a very LONG TIME but it is done only once unless you delete/rename the database.


BASIC USAGE

the s2D program can be run as a python script as it is in the s2D folder. As such it does not need installation.
One should open a terminal window and navigate to the current folder (the one where s2D.py is located).

s2D runs typing in the terminal: 

python s2D.py INPUT [a_parameter_file_with_extension_.par -plot -plotCoil -NOplot]

where INPUT should be either a protein sequence in capital 1-letter notation (only the 20 standard amino acids) or a fasta file with one or more protein sequences.
 The arguments in square brackets are optional.
 Note also that a file named s2D_parameters.txt should be in the same location as s2D.py for the program to run [or a different custom file can be provided in the command line as a_parameter_file_with_extension_.par]
  -plot   tries to generate a plot (default False). This requires the matplotlib  version >= 1.4 python module
  -NOplot   prevent from generating a plot (this option is given as well, since the default -plot behavior can be changed at the top of the s2D.py script)
  -plotCoil   plots also the coil secondary structure population (otherwise only helix and strand are plotted (and coil population is what is missing to reach a sum of 1)




THE OUTPUT OF s2D (see file s2D_out_example.txt)

if the program is run giving as input a single sequence a file s2D_out.txt is generated (and the sequence takes the default name seq1).
The file is in tab separated format (can be opened with any text editor but also with Microsoft excel and most other spread sheet editors) and contains the predicted secondary structure populations at each residues.

If the program is run giving a fasta file - e.g. my_sequences.fasta - then a file named my_sequences_s2D_out.txt is generated.
The format is the same as for the single sequence but the file may now contain multiple sequences. The beginning of each new sequence is marked with the character > and the sequence name.

if the option -plot is given (and matplotlib version >= 1.4 is installed) a plot in .pdf format is generated for each sequence given as input (sequence_name_s2D_plot.pdf).





ADVANCED USAGE AND INSTALLATION

To run s2D as a python script there is no need to perform the installation (see section 'basic usage').
However s2D can also be installed as a python module and as an excecutable (that needs to be moved in a directory in $PATH).

After customising the s2D_parameters.txt file to point to the correct location of the psiblast database (see section 'very important'),
there are two options for installing s2D.

One can type in the terminal: (install for all users, make sure before installing that the blast database is in a location accessible to every user)

sudo python setup.py 

Or, if not admin: (install only for current user, the psiblast database can be in any location accessible by the current user)

python setup.py --user


Either will install s2D as a python module. 

It will aslo generate an executable named simply 's2D' in the current folder. This can be moved in a bin directory in your $PATH (e.g. in /usr/local/bin) so that it is possible to call the s2D program from every folder just by typing in the terminal 's2D INPUT_sequence/INPUT_fasta_file [a_parameter_file_with_extension_.par -plot -plotCoil -NOplot]' (see section 'basic usage for options and input description)



USING THE PYTHON MODULE.

Explore the source file s2D_class.py for more details.
The basic usage follows. In your python script, or directly into the python console:

import s2D_class

predictor= s2D_class.s2D_method() # loads the networks and the default parameters (given in s2D_parameters.txt before installation). Different parameters can be overwritten after this call (e.g. by adding predictor.keep_pssm_profile_files=True or predictor.psiblast_ncpu=4)

predictor.run('MKAHE...G') # give a protein sequence as a str (capital letter only 20 standard amino acids) or a Biopython Seq object

predictor.print_results('my_s2D_outputfile.txt') # This will create a file in the current working directory. One can also give predictor.print_results(sys.stdout) or any other opened file stream.

# To access the output and input without printing to file:
output_numpy_ndarray = predictor.output  # a numpy ndarray object with the results of the prediction
input_numpy_ndarray_with_pssm_and_DL_results = predictor.input # DL stand for deep learning, meaning that this array, besides the input pssm matrix, contains the outputs of the first 2 networks and of the Nto1 network.

predictor.plot_results(save=True,show=False,plotCoil=plotCoil,dpi=300) # plot figure if you have matplotlib version >= 1.4. One can also give save='my_figure.png' to give a custom file name/file format.



If you run on multiple sequences init the predictor only once (by calling predictor= s2D_class.s2D_method() as loading the networks take some time) and then put predictor.run() into a loop (and possibly also predictor.print_results or save a copy of predictor.output, as the output will be overwritten at every cycle).

Also make sure that when running on multiple sequences (especially if running in parallel) each sequence has a unique name.




NB: for each sequence psiblast may generate a warning message "Warning: lcl|Query_1 seq_1: Warning: Composition-based score adjustment conditioned..." this can safely be ignored.




UNINSTALL

go to the directory where s2D has been installed (this is printed out by running 'sudo python setup.py' or 'python setup.py --user' depending on how s2D was installed in the first place. Just rerun this command if you don’t remember where s2D has been installed).
The message of setup.py will look like:

==> Installing in folder '/Library/Python/2.7/site-packages/' 

So in the example you should cd to '/Library/Python/2.7/site-packages/'

and type 'sudo rm -rf s2D*' which will remove s2D from your system.
(you don't need sudo if you previously installed with --user)

Then you might wish to deleted this folder and its content (the installation folder) and maybe also the psiblast database.
