import argparse
from os.path import join, abspath, dirname, split
from os import listdir, makedirs, getcwd
from os.path import isfile, join, exists, splitext, relpath
import glob
from subprocess import call
import random
import sys
import re

#Program directories
parsing_scripts_directory = 'parsing_src'
project_input_directory = getcwd()
default_sad_output_directory = 'sad_opt'

#Loads parsing scripts
# sys.path.insert(0, join('..', parsing_scripts_directory)+'/')
sys.path.insert(0, dirname(abspath(__file__)) + "/..")
import parsing_src.parsing_gaussian_files as parsing_gaussian_files
import parsing_src.parsing_art_files as parsing_art_files
import parsing_src.parsing_submission_files as parsing_submission_files
import parsing_src.utility as utility


default_gaussian_ext = '.inp'



#Handles file parameter passing
# #TODO add a good application description

parser = argparse.ArgumentParser(description = 'Create an input and submission file')
parser.add_argument('-sad_opt','--sad_optimization',
                    help = 'Route section optimization setting for min files (note - previous optimization will be removed from original.inp')
parser.add_argument('-s', '--sad_files', nargs='*',
                    help='specific input files to submit from project directory (e.g., sad1001')
parser.add_argument('-i','--gaussian_input', help='Gaussian input file to extract the method and basis set from')
parser.add_argument('-out','--output_file',
                    help = 'Name of the output file')
parser.add_argument('-n','--job_name',
                    help = 'Name of the job')
parser.add_argument('-w','--wall_time',
                    help = 'Wall time in hh:mm:ss')
parser.add_argument('-c','--check_one_per_cluster', action='store_true',
                    help = 'Option to automatically check one file per cluster')

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


def get_directories(starts_with):
    return glob.glob(starts_with + '*')

def get_files(directory):
    return listdir(join(getcwd(), directory))

def get_representative(file_list):
    # return file_list[0]
    return file_list[random.randint(1, len(file_list))-1]

# def create_clusters(file_starts_with):
#     files = cluster.get_art_files(file_starts_with)
#     map_to_cluster = cluster.calculate_cluster_map(files)
#     cluster.organizing_clustered_files(map_to_cluster)

def get_representatives_from_clusters(files_starts_with):
    import cluster
    starts_with = 'cluster_'
    cluster_list = get_directories(starts_with)

    file_dict = {}
    representative_files = []
    for cluster in cluster_list:
        file_dict[cluster] = get_files(cluster)
        representative_files.append(get_representative(file_dict[cluster]))
    return representative_files

def missing_saddle_arguments(args):
    if not args.check_one_per_cluster and not args.sad_files:
        return True
    return False

def get_file_number(filename):
    return re.split('(\d+)', filename)[1]

def get_min_sad_coordinates(saddle_file):

    file_counter = int(get_file_number(saddle_file))
    initial_min = 'min' + str(file_counter - 1)
    final_min = 'min' + str(file_counter)

    sad_coord = parsing_art_files.get_atomic_coordinates(saddle_file)
    initial_min_coord = parsing_art_files.get_atomic_coordinates(initial_min)
    final_min_coord = parsing_art_files.get_atomic_coordinates(final_min)
    return sad_coord + '\n' + initial_min_coord + '\n' + final_min_coord

if __name__ == '__main__':

    files_starts_with = 'sad1'

    if missing_saddle_arguments(args):
        print('Missing arguments: Ensure that -s <sad_file> is used to specify file to check')
        exit()

    #Looks for gaussian_file_location based on directory name if it was not provided as an arg
    if not args.gaussian_input:
        folder_path = dirname(getcwd())
        path, folder_name = split(getcwd())
        #TODO check next line
        gaussian_file_location = join(getcwd(), join('..', folder_name + '.inp'))

    #Gets one representative from each cluster
    if args.check_one_per_cluster:
        representative_files = get_representatives_from_clusters(files_starts_with)
        files_to_test = representative_files
    if args.sad_files:
        files_to_test = args.sad_files

    # Creates object containing all gaussian.inp information

    if args.gaussian_input:
        gaussian_file_location = join(project_input_directory, args.gaussian_input)


    input_data = parsing_gaussian_files.GaussianInput(gaussian_file_location)

    #Create min and sad output directories
    utility.create_directory(default_sad_output_directory)
    # utility.create_directory(default_sad_output_directory)

    for sad_file in files_to_test:

        #Overwrites the gaussian input file coordinates with those from the min-sad-min files
        input_data.gaussian_input_params['atom_coordinates'] = get_min_sad_coordinates(sad_file)

        input_data.write_gaussian_input_file(default_sad_output_directory, sad_file, default_gaussian_ext)

        submission_script = create_submission_file(join(default_sad_output_directory, sad_file))
        #TODO uncomment qsub when running real tests
        # call(['qsub', '-N', 'gau_opt_' + sad_file, submission_script], shell=False)


