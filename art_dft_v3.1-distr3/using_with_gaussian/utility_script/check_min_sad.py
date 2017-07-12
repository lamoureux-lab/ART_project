import argparse
from os.path import join
from os import listdir, makedirs, getcwd, chdir, rename
from os.path import isfile, join, exists
import glob
from subprocess import call
import re

#Program directories
script_directory = '../scripts'
default_input_file_directory = '../work'
default_min_output_directory = 'min_opt'
default_sad_output_directory = 'sad_opt'
default_gaussian_ext = '.com'

default_input_file = 'ethane.inp'


#Handles file parameter passing
# #TODO add a good application description

parser = argparse.ArgumentParser(description = 'Create an input and submission file')
parser.add_argument('-min_opt','--min_optimization',
                    help = 'Route section optimization setting for min files (note - previous optimization will be removed from original.inp')
parser.add_argument('-sad_opt','--sad_optimization',
                    help = 'Route section optimization setting for sad files (note - previous optimization will be removed from original.inp')
parser.add_argument('-f', '--input_files', nargs='*',
                    help='specific input files to submit from project directory')
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
        frequency = []
        for line in f:
            if line.startswith(" Frequencies"):
                frequency.append(line)
                print(line)

    for line in frequency:
        freq = line.split()
        check = float(freq[2])
        if check < 0:
            print("Saddle")

def load_gaussian_input(input_file, gaussian_input_params):
    """
    Loads the values gaussian parameters into a dictionary object
    :param gaussian_input_params
    :return:
    """
    params = gaussian_input_params
    section_number = 1

    with open(input_file) as input:

        #The periodic_table_dict will be loaded once if conversion is needed from atomic symbol to number
        periodic_table_dict = None

        for line in input:
                # This will load the relevant information from a gaussian input file expecting:
                # 1. A number of % lines (or maybe none) representing the link 0 section
                # 2. A number of # line describing the job representing the route section
                # 3. Blank line, title line, blank line
                # 4. Charge + Spin multiplicity line (two integers)
                # 5. Atom coordinates

                # Remove all comments from the gaussian input file
                comment_set = False
                line = line.split("!")
                if len(line) > 1:
                    comment_set = True
                line = line[0]
                if comment_set:
                    line = line + '\n'

                #Check for an empty line to increment section
                if line.strip() == '':
                    if not comment_set:
                        section_number += 1
                    continue

                #Link0 and route section
                if section_number == 1:
                    #For case issues
                    line = line.lower()
                    #Link0
                    if line.find('%') != -1:
                        params['link0_section'] = params['link0_section'] + line + " "
                    #Route section
                    elif line.find('#') != -1:
                        # Removes opt and force, which will be put back in by the ART code at different stages

                        # Removes route parameters with and without options
                        # line = option_removal_helper('opt', line)
                        # line = option_removal_helper('force', line)
                        # line = option_removal_helper('nosymm', line)

                        if line.strip() != '#':
                            params['route_section'] = params['route_section'] + line

                #Title section
                elif section_number == 2:
                    params['title'] = line.strip()

                #Molecule specification section
                if section_number == 3:
                    string_integers = line.split()
                    params['charge'] = int(string_integers[0])
                    params['multiplicity'] = int(string_integers[1])
                    section_number += 1
                    break

        return params

def create_new_gaussian_input(gaussian_input_params, new_directory, file_name, gaussian_ext):

    params = gaussian_input_params

    header = (params['link0_section']) \
             + (params['route_section'].rstrip()) + ('\n') \
             + (params['title'] + '\n') + ('\n') \
             + (str(params['charge']) + ' ' + str(params['multiplicity']) + '\n')

    coordinates = params['atom_coordinates']

    gaussian_input_params

    with open(join(new_directory, file_name + gaussian_ext), 'w') as output:
        output.write(header)
        output.write(coordinates)

def create_output_directory(directory):
        """
        Creates a new directory if it does not already exist
        :param directory:
        :return:
        """
        if not exists(directory):
            makedirs(directory)
            return True
        return False

if __name__ == '__main__':
    #Load gaussian header file informationinto dictionary
    gaussian_input_params = {'link0_section': '', 'route_section': '',
                             'title': '', 'natoms': 0, 'charge': None, 'multiplicity': None, 'atom_coordinates': ''}
    gaussian_input_params = load_gaussian_input(join(default_input_file_directory,default_input_file), gaussian_input_params)


    #Create min and sad output directories
    create_output_directory(default_min_output_directory)
    create_output_directory(default_sad_output_directory)


    #for min file in directory :
    file_list = glob.glob('min1*')
    for min_file in file_list:
        gaussian_input_params = get_atomic_coordinates(gaussian_input_params, min_file)
        create_new_gaussian_input(gaussian_input_params, default_min_output_directory, min_file, default_gaussian_ext)
        submission_script = create_submission_file(join(default_min_output_directory, min_file))
        call(['qsub', '-N', 'gau_opt_' + min_file, submission_script], shell=False)

    file_list = glob.glob('sad1*')
    for sad_file in file_list:
        gaussian_input_params = get_atomic_coordinates(gaussian_input_params, sad_file)
        submission_script = (gaussian_input_params, default_sad_output_directory, sad_file,
                                      default_gaussian_ext)
        submission_script = create_submission_file(join(default_sad_output_directory, sad_file))
        call(['qsub', '-N', 'gau_opt_' + sad_file, submission_script], shell=False)

    create_new_gaussian_input(gaussian_input_params, default_sad_output_directory, sad_file, default_gaussian_ext)

