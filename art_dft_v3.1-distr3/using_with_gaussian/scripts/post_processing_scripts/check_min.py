import argparse
from os.path import join, abspath, dirname, split
from os import listdir, makedirs, getcwd
from os.path import isfile, join, exists, splitext, relpath
import glob
from subprocess import call
import random
import sys
import cluster_old


#Program directories
parsing_scripts_directory = 'parsing_src'
project_input_directory = getcwd()
default_min_output_directory = 'min_opt'
default_gaussian_ext = '.inp'



parser = argparse.ArgumentParser(description = 'Create an input and submission file')
parser.add_argument('-min_opt','--min_optimization',
                    help = 'Route section optimization setting for min files (note - previous optimization will be removed from original.inp')
parser.add_argument('-m', '--min_files', nargs='*',
                    help='specific input files to submit from project directory (e.g., min1000')
parser.add_argument('-i','--art_input', help='ART input file to extract the method and basis set from')
parser.add_argument('-out','--output_file',
                    help = 'Name of the output file')
parser.add_argument('-j','--job_name', default = 'min_opt',
                    help = 'Name of the job')
parser.add_argument('-w','--wall_time', default = '1:00:00',
                    help = 'Wall time in hh:mm:ss')
parser.add_argument('-mem','--memory', default = '2000MB',
                    help = 'Memory in MB eg. 2000MB')
parser.add_argument('-np','--nodes_proc', default = 'nodes=1:ppn=4',
                    help = 'number of nodes and processors, eg. nodes=1:ppn=4')
parser.add_argument('-tol','--dist_tol', type = float, default = 0.01,
                    help = 'Distance tolerance value eg. 0.01')
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
#PBS -l ''' + np + ''' 
#PBS -l mem=''' + me + ''' 
#PBS -l walltime=''' + wt + '''
#PBS -N ''' + jn + '''
            
# Adjust the mem and ppn above to match the requirements of your job 
# Sample Gaussian job script 
cd $PBS_O_WORKDIR 
echo "Current working directory is `pwd`" 
echo "Running on `hostname`" 
echo "Starting run at: `date`" 
            
# Set up the Gaussian environment using the module command: 
module load gaussian \n # Run Submission 
g09 ''' + filename + '.inp > ' + filename + '.log\n'
'module load python\n\n'
+ 'python ' + join(dirname(relpath(__file__)), 'check_min_freq.py') +' < ' + filename + '.log' + ' > ' + filename + '_results.txt' '\n\n')

    return submission_script

def get_number_of_header_lines(input_file):
    j = 0
    with open(input_file) as f:
        for line in f:
            if line.startswith('%'):
                j = j + 1
            elif line.startswith('#'):
                j = j + 1
        k = j + 3
        return k

def get_gaussian_header(input_file):
    i = 0
    with open(input_file) as f:
        gaussian_header = ''
        for line in f:
            gaussian_header = gaussian_header + line
            i = i + 1
            if i > numb_of_head:
                break
        return gaussian_header

def get_charge_multiplicity(input_file):
    with open(input_file) as f:
        for i in range(0,numb_of_head):
            _ = f.readline()
        charge_multiplicity = f.readline()
        return charge_multiplicity


def get_min_coordinates(min_file):

    min_coord = get_atomic_coordinates(min_file)
    return min_coord 

def create_gaussian_input_file(min_file):
    gaussian_input_file = min_file + '.inp'
    with open(gaussian_input_file, 'w') as f:
        f.write(gaussian_header + coords + '\n\n')

    return gaussian_input_file

def create_directory(directory):
    if not exists(directory):
        makedirs(directory)
        return True
    return False


files_to_test = args.min_files

if args.check_one_per_cluster:
	tolerance = args.dist_tol
	files = cluster_old.get_art_files('min1')
    	map1 = cluster_old.calculate_cluster_map(files, tolerance)
    	file_dict = cluster_old.make_json_list(map1)
	
	files_to_test = []
    	for key in file_dict.iterkeys():
        	files_to_test.append(key)

ART_input_file = args.art_input
numb_of_head = get_number_of_header_lines(ART_input_file)
wt = args.wall_time
jn = args.job_name
np = args.nodes_proc
me = args.memory

if __name__ == '__main__':

    gaussian_header = get_gaussian_header(ART_input_file)
    charge_multiplicity = get_charge_multiplicity(ART_input_file)
    create_directory(default_min_output_directory)

    for min_file in files_to_test:
        submission_script = create_submission_file(join(default_min_output_directory, min_file))
        coords = get_min_coordinates(min_file)
        create_gaussian_input_file(join(default_min_output_directory, min_file))
        call(['qsub', '-N' + jn + '_' + min_file, submission_script], shell=False)


