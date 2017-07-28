
# This is an adapter to allow the user to supply only a Gaussian formatted input file and have this run through ART
# For reference: http://gaussian.com/input/

import re
import argparse
import sys
from os import listdir, makedirs, getcwd, chdir, rename
from os.path import isfile, join, exists, dirname, relpath
from shutil import copy, rmtree
from subprocess import call

#Program directories
script_directory = '../' #Relative script directory to be used by parsing_scr
program_data_directory = 'program_data'
project_input_directory = getcwd()
project_output_directory = project_input_directory
default_submission = 'GREX' #GREX or PSI

#Data
periodic_table_data = 'periodic_table_data.csv'

#Program files

art_executable_rel_location = '../../source/artdft'
#Gets the relative path which will be used to call ART. This goes back an additional directory since new structure
#directories are being made where all scripts for a particular job are stored
art_executable_location = join('../', relpath(join(dirname(__file__), art_executable_rel_location), getcwd()))

gaussian_execution_script = 'execute_gaussian.sh'
refconfig_filename = 'refconfig.dat'
gaussian_art_filename = 'gaussian_art.sh'     #Shell script containing the configuration parameters for the ART application
grex_submission_script = 'gauss_grex.sub'
psi_submission_script =  'gauss_psi.sub'

#Loads parsing scripts
import parsing_src.parsing_submission_files as parsing_submission_files
import parsing_src.parsing_gaussian_files as parsing_gaussian_files
import parsing_src.parsing_art_files as parsing_art_files
import parsing_src.utility as utility

#Default configurations
file_counter_start = 1000

#Handles file parameter passing
def argument_parser():
    parser = argparse.ArgumentParser(description='Gaussian_ART ... ') #TODO add a good application description
    parser.add_argument('-n', '--new_gaussian_header', action='store_true',
                        help='this will apply new gaussian.inp header data while continuing from previous coordinates '
                             'obtained by ART')
    parser.add_argument('-f', '--input_files', nargs='*',
                        help='specific input files to submit from project directory')
    parser.add_argument('-s', '--submission_type', choices=['GREX', 'PSI'], default=default_submission,
                        help='GREX or PSI submissions supported. ' + default_submission + ' is the default')
    parser.add_argument('-r', '--reset', action='store_true',
                        help='resets the output directory instead of continuing from latest referenced configuration')
    parser.add_argument('-t', '--time',
                        help='changes the submissions job time (e.g., -t \'walltime=010:00:00\' ')
    parser.add_argument('-m', '--memory',
                        help='changes the submissions job memory cap (e.g., -m \'mem=8000MB\' ')
    parser.add_argument('-p', '--processor_info',
                        help='changes the submissions job processor cap (e.g., -p \'nodes=1:ppn=8\' ')
    parser.add_argument('-o', '--optimize', action='store_true',
                        help='sets submission job memory and number of processors to what is specified in the file\n '
                             '(note: specifying --memory or --processor_info will take precedence)')
    parser.add_argument('-e', '--max_number_of_events', type=int,
                        help='')
    parser.add_argument('-c', '--central_atom', type=int,
                        help='')
    parser.add_argument('-d', '--radius_initial_deformation', type=int,
                        help='')
    parser.add_argument('-xyz', '--write_xyz', action='store_true',
                        help='')

    args = parser.parse_args()

    return args


def get_gaussian_input_files(input_file_directory, select_input_files):
    # Gets file names from gaussian input files directory

    #Load all files in that directory
    files_in_directory = [f for f in listdir(input_file_directory) if isfile(join(input_file_directory, f))]

    gaussian_files_in_directory = []
    for filename in files_in_directory:
        extension_list = ['.com', '.inp', '.input', '.gjf']
        if [x for x in extension_list if filename.endswith(x)]:
            gaussian_files_in_directory.append(filename)

    #Adds only those files selected by the user
    final_file_list = []
    if (select_input_files):
        for file_name in select_input_files:
            if file_name in gaussian_files_in_directory:
                final_file_list.append(file_name)
            else:
                print 'File not found: ' + file_name
    else:
        final_file_list = gaussian_files_in_directory

    print('Preparing the following gaussian input file(s) for submission: ' )
    for file in final_file_list:
        print file

    return final_file_list


if __name__ == "__main__":

    args = argument_parser()

    # Get input file names:

    print 'Loading files from directory:  ' + project_input_directory
    input_files = get_gaussian_input_files(project_input_directory, args.input_files)


    is_new_header = args.new_gaussian_header

    # Creates an output directory if not already present
    output_directory = project_output_directory
    if utility.create_directory(output_directory):
        start_from_scratch = True

    # Initialize gaussian.inp file parameters
    input_data = {}
    for input_file_name in input_files:
        # Reset variable will be passed by user to restart an output directory, erasing all data
        start_from_scratch = args.reset
        structure = input_file_name.split('.')[0]

        # Creates object containing all gaussian.inp information
        input_data[structure] = parsing_gaussian_files.gaussian_input(join(project_input_directory, input_file_name))
        input_data[structure].remove_route_parameter(['opt', 'force', 'nosymm'])

        # Makes symbols from a gaussian input file atomic numbers
        # input_data[structure].symbol_to_atomic_number(join(program_data_directory, periodic_table_data))

        # Creates output directories for each structure if not already present
        structure_output_directory = join(output_directory, structure)
        if utility.create_directory(structure_output_directory):
            start_from_scratch = True

        # Checks that an existing structure directory has the correct files to continue running ART
        if not start_from_scratch:
            filename_list = ['.'+gaussian_execution_script, gaussian_art_filename, refconfig_filename, 'filecounter']
            message = 'Reset this directory with latest .inp header and configuration settings (y) or skip structure (n) ' \
                   '\nNote: This will conserve previous coordinates and not erase existing min and sad files. '
            #TODO modify this to start from only modify header when possible (is_new_header)
            start_from_scratch = utility.check_for_missing_files_reset(structure_output_directory, filename_list, message)

        art_environment = {}
        art_environment['natom'] = (input_data[structure].gaussian_input_params['natoms'])
        art_environment['max_number_of_events'] = args.max_number_of_events
        art_environment['central_atom'] = args.central_atom
        art_environment['write_xyz'] = args.write_xyz
        art_environment['radius_initial_deformation'] = args.radius_initial_deformation

        # This is only called when the user wants to restart the structure or it was not already set
        if start_from_scratch:

            atom_coordinates = input_data[structure].gaussian_input_params['atom_coordinates']
            parsing_art_files.create_ref_config(refconfig_filename, atom_coordinates, structure_output_directory)

            output_structure_directory = join(output_directory, structure)
            parsing_art_files.set_env_config(script_directory, gaussian_art_filename, refconfig_filename,
                                             output_structure_directory, art_environment, art_executable_location)

            input_data[structure].insert_gaussian_file_header_into_execute_gaussian(script_directory, gaussian_execution_script,
                                                                                    structure_output_directory, add_param_flag=True)

            parsing_art_files.create_file_counter(file_counter_start, structure_output_directory)

        elif is_new_header:
            # TODO make sure that new structure data can be used on old coordinates by default if not resetting
            # TODO Simply change the scripts without removing refconfig/min/sad files

            output_structure_directory = join(output_directory, structure)
            parsing_art_files.set_env_config(script_directory, gaussian_art_filename, refconfig_filename,
                                             output_structure_directory, art_environment, art_executable_location)

            input_data[structure].insert_gaussian_file_header_into_execute_gaussian(script_directory, gaussian_execution_script,
                                                                                    structure_output_directory, add_param_flag=True)


        submission_type = args.submission_type
        if submission_type == 'GREX':
            #TODO handle the case where no structure directory is found/not running submission for a particular structure
            link0_section = input_data[structure].gaussian_input_params['link0_section']
            parsing_submission_files.set_submission_script(script_directory, grex_submission_script, structure_output_directory,
                                                           link0_section, optimize = args.optimize, time = args.time,
                                                           memory = args.memory, processor_info = args.processor_info)
            print 'Running submission file for: ' + structure
            wd = getcwd()
            chdir(join(wd, structure_output_directory))
            #TODO uncomment for tests on GREX
            call(['qsub', '-N','gau_art_'+ structure, grex_submission_script], shell=False)
            chdir(wd)
        else:
            print 'Only GREX submission is currently available'
        #TODO handle other submission types