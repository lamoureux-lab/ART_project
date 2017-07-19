import argparse
from os.path import join, abspath, dirname
from os import listdir, makedirs, getcwd
from os.path import isfile, join, exists, splitext
import glob
from subprocess import call
import re
import sys

#Program directories
parsing_scripts_directory = 'parsing_src'
project_input_directory = getcwd()
default_min_output_directory = 'min_opt'
# default_sad_output_directory = 'sad_opt'

#Loads parsing scripts
# sys.path.insert(0, join('..', parsing_scripts_directory)+'/')
sys.path.insert(0, dirname(abspath(__file__)) + "/..")
import parsing_src.parsing_gaussian_files as parsing_gaussian_files
import parsing_src.parsing_art_files as parsing_art_files
import parsing_src.parsing_submission_files as parsing_submission_files
import parsing_src.utility as utility

# TODO sort out ValueError: Attempted relative import in non-package error with the following to import a less hacky way
# from ..parsing_src import parsing_gaussian_files
# from ..parsing_src import parsing_art_files
# from ..parsing_src import parsing_submission_files
# from ..parsing_src import utility

default_gaussian_ext = '.inp'

default_input_file = '../ethane.inp'


#Handles file parameter passing
# #TODO add a good application description

parser = argparse.ArgumentParser(description = 'Create an input and submission file')
parser.add_argument('-min_opt','--min_optimization',
                    help = 'Route section optimization setting for min files (note - previous optimization will be removed from original.inp')
# parser.add_argument('-sad_opt','--sad_optimization',
#                     help = 'Route section optimization setting for sad files (note - previous optimization will be removed from original.inp')
parser.add_argument('-f', '--input_files', nargs='*',
                    help='specific input files to submit from project directory')
parser.add_argument('-i','--gaussian_input', help='Gaussian input file to extract the method and basis set from')
parser.add_argument('-out','--output_file',
                    help = 'Name of the output file')
args = parser.parse_args()

def get_atomic_coordinates(gaussian_input_params, ART_output_file):
    with open(ART_output_file) as f:
        gaussian_input_params['atom_coordinates'] = ''
        for i in range(0,3):
            _= f.readline()
        for line in f:
            gaussian_input_params['atom_coordinates'] = gaussian_input_params['atom_coordinates'] + line
        return gaussian_input_params

def create_submission_file(filename):
    submission_script = filename + '.sub'
    with open(submission_script,'w') as n:
        n.write('''
#!/bin/bash
#PBS -S /bin/bash 
#PBS -l nodes=1:ppn=4 
#PBS -l mem=1800MB 
#PBS -l walltime=5:00:00 
#PBS -N Ethane_Opt 
            
# Adjust the mem and ppn above to match the requirements of your job 
# Sample Gaussian job script 
cd $PBS_O_WORKDIR 
echo "Current working directory is `pwd`" 
echo "Running on `hostname`" 
echo "Starting run at: `date`" 
            
# Set up the Gaussian environment using the module command: 
module load gaussian \n # Run Submission 
g09 ''' + filename + '''.com''')

    return submission_script


def check_min_or_sad(logfile):
    with open(logfile)as f:
        if logfile.startswith('min'):
            frequency = []
            for line in f:
                if line.startswith(" Frequencies"):
                    frequency.append(line)
                    print(line)

            for line in frequency:
                freq = line.split()
                check = float(freq[2])
                if check < 0:
                    print("Optimization failed")


if __name__ == '__main__':

    # Creates object containing all gaussian.inp information
    input_data = parsing_gaussian_files.gaussian_input(join(project_input_directory, default_input_file))

    #Create min and sad output directories
    utility.create_directory(default_min_output_directory)
    # utility.create_directory(default_sad_output_directory)

    #for min file in directory :
    filetype = 'min1'
    # Excludes .xyz files from cluster parsing
    file_list = [fn for fn in glob.glob(filetype + '*')
             if not splitext(fn)[1] == '.xyz']

    for min_file in file_list:
        gaussian_input_params = get_atomic_coordinates(input_data.gaussian_input_params, min_file)


        # create_new_gaussian_input(gaussian_input_params, default_min_output_directory, min_file, default_gaussian_ext)
        input_data.write_header(default_min_output_directory, min_file, default_gaussian_ext)
        input_data.write_coordinates(default_min_output_directory, min_file, default_gaussian_ext)

        submission_script = create_submission_file(join(default_min_output_directory, min_file))
        #TODO uncomment qsub when running real tests
        # call(['qsub', '-N', 'gau_opt_' + min_file, submission_script], shell=False)


