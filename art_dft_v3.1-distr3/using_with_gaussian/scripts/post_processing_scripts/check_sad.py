import argparse
from os.path import join, abspath, dirname, split
from os import listdir, makedirs, getcwd
from os.path import isfile, join, exists, splitext, relpath
import glob
from subprocess import call
import random
import sys
import re

project_input_directory = getcwd()
default_sad_output_directory = 'sad_opt'
default_gaussian_ext = '.inp'


parser = argparse.ArgumentParser(description = 'Create an input and submission file')
parser.add_argument('-sad_opt','--sad_optimization',
                    help = 'Route section optimization setting for min files (note - previous optimization will be removed from original.inp')
parser.add_argument('-s', '--sad_files', nargs='*',
                    help='specific input files to submit from project directory (e.g., sad1001')
parser.add_argument('-i','--art_input', help='ART input file to extract the method and basis set from')
parser.add_argument('-out','--output_file',
                    help = 'Name of the output file')
parser.add_argument('-n','--job_name',
                    help = 'Name of the job')
parser.add_argument('-w','--wall_time',
                    help = 'Wall time in hh:mm:ss')
parser.add_argument('-c','--check_one_per_cluster', action='store_true',
                    help = 'Option to automatically check one file per cluster')


args = parser.parse_args()


def get_atomic_coordinates(ART_output_file):
    with open (ART_output_file) as f:
        atomic_coords = ''
        for i in range (0,3):
            _ = f.readline()
        for line in f:
            atomic_coords = atomic_coords + line
        return atomic_coords

def create_submission_file(filename):
    submission_script = filename + '.sub'
    with open(submission_script,'w') as n:
        n.write('''
#!/bin/bash
#PBS -S /bin/bash 
#PBS -l nodes=1:ppn=4 
#PBS -l mem=1800MB 
#PBS -l walltime=''' + wall_time + '''
#PBS -N''' + job_name + '''
            
# Adjust the mem and ppn above to match the requirements of your job 
# Sample Gaussian job script 
cd $PBS_O_WORKDIR 
echo "Current working directory is `pwd`" 
echo "Running on `hostname`" 
echo "Starting run at: `date`" 
            
# Set up the Gaussian environment using the module command: 
module load gaussian \n # Run Submission 
g09 ''' + filename + '.inp > ' + filename + '.log\n'
 + 'python ' + join(dirname(relpath(__file__)), 'check_sad_freq.py') +' < ' + filename + '.log' + ' > ' + filename + '_results.txt' '\n')

    return submission_script

def get_gaussian_header(input_file):
    i = 0
    j = 0
    with open(input_file) as f:
        gaussian_header = ''
        for line in f:
            if line.startswith('%'):
                j = j + 1
            elif line.startswith('#'):
                j = j + 1
            gaussian_header = gaussian_header + line
            i = i + 1
            if i > j + 3: #Because space-title-space-charge_multiplicity
                break
        return gaussian_header

def get_charge_multiplicity(input_file):
    j = 0
    with open(input_file) as f:
        for line in f:
            if line.startswith('%'):
                j = j + 1
            elif line.startswith('#'):
                j = j + 1
            k = j + 3
        for i in range(0,k):
            _ = f.readline()
        charge_multiplicity = f.readline()
        return charge_multiplicity


def get_file_number(filename):
    return re.split('(\d+)', filename)[1]

def get_min_sad_coordinates(saddle_file):

    file_counter = int(get_file_number(saddle_file))
    initial_min = 'min' + str(file_counter - 1)
    final_min = 'min' + str(file_counter)

    sad_coord = get_atomic_coordinates(saddle_file)
    initial_min_coord = get_atomic_coordinates(initial_min)
    final_min_coord = get_atomic_coordinates(final_min)
    return initial_min_coord + '\n' + 'Title Card Required' + '\n\n' + charge_multiplicity + final_min_coord + '\n' + 'Title Card Required' + '\n\n' + charge_multiplicity + sad_coord + '\n' + '\n'

def create_gaussian_input_file(saddle_file):
    gaussian_input_file = saddle_file + '.inp'
    with open(gaussian_input_file, 'w') as f:
        f.write(gaussian_header + coords)

    return gaussian_input_file

def create_directory(directory):
    if not exists(directory):
        makedirs(directory)
        return True
    return False


files_to_test = args.sad_files
ART_input_file = args.art_input


if __name__ == '__main__':

    gaussian_header = get_gaussian_header(ART_input_file)
    charge_multiplicity = get_charge_multiplicity(ART_input_file)
    create_directory(default_sad_output_directory)

    for sad_file in files_to_test:
        submission_script = create_submission_file(join(default_sad_output_directory, sad_file))
        coords = get_min_sad_coordinates(sad_file)
        create_gaussian_input_file(join(default_sad_output_directory, sad_file))
        #call(['qsub', '-N', 'gau_opt_' + sad_file, submission_script], shell=False)
    

    





