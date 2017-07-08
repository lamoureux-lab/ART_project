
# This is an adapter to allow the user to supply only a Gaussian formatted input file and have this run through ART
# For reference: http://gaussian.com/input/

import re
import argparse
import sys
from os import listdir, makedirs, getcwd, chdir, rename
from os.path import isfile, join, exists
from shutil import copy, rmtree
from subprocess import call

#Program directories
script_directory = 'scripts'
program_data_directory = 'program_data'
default_input_file_directory = 'work'
default_output_directory = default_input_file_directory
default_submission = 'GREX' #GREX or PSI

#Program files
gaussian_execution_script = 'execute_gaussian.sh'
refconfig_filename = 'refconfig.dat'
gaussian_art_filename = 'gaussian_art.sh'     #Shell script containing the configuration parameters for the ART application
periodic_table_data = 'periodic_table_data.csv'
grex_submission_script = 'gauss_grex.sub'
psi_submission_script =  'gauss_psi.sub'

#Default configurations
file_counter_start = 1000

#Handles file parameter passing
parser = argparse.ArgumentParser(description='Gaussian_ART ... ') #TODO add a good application description
parser.add_argument('-d', '--project_directory',
                    help='project directory containing gaussian input files', default=default_input_file_directory)
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
args = parser.parse_args()


def set_submission_script(submission_script, structure_output_directory, link0_section, optimize=False, time = None, memory = None, processor_info = None):
    """
        Sets submissions parameters in gaussian_art.sh which contains the general configuration for the ART environment

        :param natoms:
        :return:
        """
    global script_directory
    submission_script_template = join(script_directory, submission_script)
    link0_arr = link0_section.splitlines()

    with open(submission_script_template) as input:
        updated_text = ''
        for line in input:
            original_line = line
            updated_line = ''

            # Sets the time
            time_flag = '#PBS -l walltime'
            if time_flag in line:
                if time:
                    updated_line = '#PBS -l ' + time + '\n'

            # Sets the memory
            memory_flag = '#PBS -l mem'
            if memory_flag in line:
                if memory:
                    updated_line = '#PBS -l ' + memory + '\n'
                elif optimize:
                    for item in link0_arr:
                        link0_memory = extract_link0_parameter('mem', item)
                        #Checks if the memory is set in the gaussian input file
                        if link0_memory:
                            updated_line = '#PBS -l mem=' + link0_memory + '\n'

            # Sets the processor information
            processor_flag = '#PBS -l nodes'
            if processor_flag in line:
                if processor_info:
                    updated_line = '#PBS -l ' + processor_info + '\n'
                elif optimize:
                    # % nproc = 12
                    for item in link0_arr:
                        link0_nproc = extract_link0_parameter('nproc', item)
                        # Checks if the nproc is set in the gaussian input file
                        if link0_nproc:
                            updated_line = '#PBS -l nodes=1:ppn=' + link0_nproc + '\n'

            # updates the submission file string
            if updated_line:
                updated_text = updated_text + updated_line
            else:
                updated_text = updated_text + original_line

    # overwrites the submissions file with appropriate values from the gaussian input file
    sub = open(join(structure_output_directory, submission_script), 'w+')
    sub.write(updated_text)
    # print updated_text
    sub.close()

def extract_link0_parameter(flag, link0_line):
    """
    Gets the parameter assigned for memory or number of processors in a gaussian.inp file
    :return:
    """
    if flag in link0_line:
        return link0_line.split("=")[1].strip()
    return None



def get_gaussian_input_files(input_file_directory, select_input_files):
    # Gets file names from gaussian input files directory

    #Load all files in that directory
    files_in_directory = [f for f in listdir(input_file_directory) if isfile(join(input_file_directory, f))]

    #Adds only those files selected by the user
    final_file_list = []
    if (select_input_files):
        for file_name in select_input_files:
            if file_name in files_in_directory:
                final_file_list.append(file_name)
            else:
                print 'File not found: ' + file_name
    else:
        final_file_list = files_in_directory

    print('Preparing the following gaussian input file(s) for submission: ' )
    for file in final_file_list:
        print file

    return final_file_list

def create_file_counter(file_counter_start, output_directory):
    """
    Creates a file counter that is used by the ART application to indicate the min/saddle file event
    :param file_counter_start:
    :param output_directory:
    """
    file_counter = open(join(output_directory, 'filecounter'), 'w+')
    file_counter.write('Counter:      ' + str(file_counter_start))
    file_counter.close()

def create_ref_config(atom_coordinates, output_directory):
    """
    Sets up refconfig.dat file which will be used as a starting reference by ART for atom coordinates before it is
    overwritten and updated

    :param refconfig:
    :param gaussian_input_params:
    :return:
    """
    global refconfig_filename

    config = open(join(output_directory, refconfig_filename), 'w+')     #Overwrites or creates a new file if it doesn't exist
    config.write('run_id:         1000\n')
    config.write('total_energy:   0\n')     #Placeholder, as this will be optimized by ART to the correct value
    #TODO see about removing these as they are not necessary for gaussian
    config.write('S   100.000000000000        100.000000000000        100.000000000000\n')
    config.write(atom_coordinates)
    config.close()

def set_env_config(natoms):
    """
    Sets critical parameters in gaussian_art.sh which contains the general configuration for the ART environment

    :param natoms:
    :return:
    """
    global script_directory
    global refconfig_filename
    global gaussian_art_filename

    with open(join(script_directory, gaussian_art_filename)) as input:

        updated_text = ''
        for line in input:
            original_line = line
            updated_line = ''

            #Remove comments from the strings
            line = line.split("#")
            comment_set = False
            if len(line) > 1:
                comment_set = True

            #Checks that GAU is set for energy calculation and sets it if it is not
            if line[0].find('setenv ENERGY_CALC') != -1:
                updated_line = 'setenv ENERGY_CALC GAU      '

            #Sets the number of atoms
            elif line[0].find('setenv NATOMS') != -1:
                updated_line = 'setenv NATOMS     '+ str(natoms) + '         '

            #Sets reference configuration file that is updated throughout
            elif line[0].find('setenv REFCONFIG') != -1:
                updated_line = 'setenv REFCONFIG        ' + refconfig_filename + '             '

            # Puts back configuration comments
            if comment_set and updated_line:
                updated_line = updated_line + '#' + line[1]
            elif updated_line:
                updated_line = updated_line + '\n'

            # updates the configuration file string
            if updated_line:
                updated_text = updated_text + updated_line
            else:
                updated_text = updated_text + original_line

    # overwrites the configuration file with appropriate values from the gaussian input file
    config = open(join(script_directory, gaussian_art_filename), 'w+')
    config.write(updated_text)
    config.close()


def create_gaussian_file_header(gaussian_input_params, structure_output_directory):
    """
    Inserts a gaussian header file into the bash script that will be writing to
    the art2gaussian.inp file and calling gaussian

    :param gaussian_input_params:
    :return:
    """
    global script_directory
    global gaussian_execution_script
    gaussian_execution_script_global = join(script_directory, gaussian_execution_script)
    copy(gaussian_execution_script_global, structure_output_directory)
    gaussian_execution_script_local = join(structure_output_directory, gaussian_execution_script)

    params = gaussian_input_params

    # creates a regular expression pattern that will isolate the header insertion point
    insertion_point = re.compile('#gaussian-header-begin.*?#gaussian-header-end \(DO NOT REMOVE\)', re.DOTALL)

    # open file
    f = open(gaussian_execution_script_local, 'r')
    initial_script = f.read()
    f.close()

    output = open(gaussian_execution_script_local, 'w+')

    # Builds a header for the gaussian.inp file, <OPTION> Flag will be replaced by opt or force

    start_shell_script_marker = '#gaussian-header-begin (DO NOT REMOVE) \n'
    header = 'header=\'' \
             + (params['link0_section']) \
             + (params['route_section'].rstrip() + ' <PARAM>' + '\n') + ('\n') \
             + (params['title'] + '\n') + ('\n') \
             + (str(params['charge']) + ' ' + str(params['multiplicity']) + '\n' + '\'')

    coordinate_line_start = '\n#Coordinate line where data begins\n' \
             + 'coorLineNumber=' + str(get_coordinate_line_number(header))

    title_line_start = '\n#Title of the gaussian input file\n' \
             + 'title=' + '\'' + str(gaussian_input_params['title']) + '\''

    end_shell_script_marker = '\n\n#gaussian-header-end (DO NOT REMOVE) '

    output.close()

    # Writes the header as a string variable in the gaussian execution script
    new_script_lines = start_shell_script_marker + header + coordinate_line_start + title_line_start + end_shell_script_marker
    script_with_header = insertion_point.sub(new_script_lines, initial_script)
    f = open(gaussian_execution_script_local, 'w')
    f.write(script_with_header)
    f.close()

    # Hides script to reduce clutter
    rename(gaussian_execution_script_local, join(structure_output_directory, '.' + gaussian_execution_script))


def get_coordinate_line_number(str):
    return len(str.split('\n'))

def load_input(input_file, gaussian_input_params):
    """
    Loads the values gaussian parameters into a dictionary object, while removing certain route parameters
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
                        line = option_removal_helper('opt', line)
                        line = option_removal_helper('force', line)
                        line = option_removal_helper('nosymm', line)

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
                    continue

                #Atom coordinates section
                if section_number == 4:
                    params['natoms'] = params['natoms'] + 1

                    #Checks if the atom is in non-numerical format (e.g., H instead of 1)
                    temp_line = line.split()
                    element = temp_line[0]
                    if not element.isdigit():
                        #Loads the periodic table data if it is not already set
                        if not periodic_table_dict:
                            periodic_table_dict = load_periodic_table()
                        #Replaces atomic symbols with atomic numbers
                        temp_line[0] = periodic_table_dict[element]
                        line = '    '.join(temp_line) + '\n'

                    params['atom_coordinates'] = params['atom_coordinates'] + line

                if section_number == 5:
                    break

        return params


def load_periodic_table():
    """
    Loads a periodic table csv file to extract the atomic symbols and numbers to create a mapping between these
    :return: A dictionary mapping atomic symbol to number
    """
    global periodic_table_data
    global program_data_directory

    atomic_symbol_number_map = {}
    with open(join(program_data_directory, periodic_table_data)) as input:
        for line in input:
            line = line.split(',')
            atomic_symbol_number_map[line[1].strip()] = line[0].strip()

    return atomic_symbol_number_map


def option_removal_helper(option_type, line):
    #Handles different options to keyword structures for removal to have them dynamically assigned by ART:
    #keyword = option
    #keyword(option)
    #keyword=(option1, option2)
    #keyword(option1, option2)

    #For keyword = option
    line = re.sub(r'(\s)' + re.escape(option_type) +'=\w+', r'\1', line)
    line = re.sub(r'(#)' + re.escape(option_type) +'=\w+', r'\1', line)
    line = re.sub(r'(\s)' + re.escape(option_type) +' =\w+', r'\1', line)
    line = re.sub(r'(#)' + re.escape(option_type) +' =\w+', r'\1', line)
    line = re.sub(r'(\s)' + re.escape(option_type) +' = \w+', r'\1', line)
    line = re.sub(r'(#)' + re.escape(option_type) +' = \w+', r'\1', line)
    line = re.sub(r'(\s)' + re.escape(option_type) +' = \w+', r'\1', line)
    line = re.sub(r'(#)' + re.escape(option_type) +' = \w+', r'\1', line)

    #For keyword(option..)
    line = re.sub(r'(\s)' + re.escape(option_type) + '\([^()]*\)', '', line)
    line = re.sub(r'(#)' + re.escape(option_type) + '\([^()]*\)', '', line)

    #For keyword=(option..)
    line = re.sub(r'(\s)' + re.escape(option_type) + '=\([^()]*\)', '', line)
    line = re.sub(r'(#)' + re.escape(option_type) + '=\([^()]*\)', '', line)

    #For lone keyword
    line = re.sub(r'(\s)' + re.escape(option_type), r'\1', line)
    line = re.sub(r'(#)' + re.escape(option_type), r'\1', line)

    return line

def query_yes_no(question, default="yes"):
    """
    From Python recipe http://code.activestate.com/recipes/577058/
    Asks a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def create_directory(directory):
    """
    Creates a new directory if it does not already exist
    :param directory:
    :return:
    """
    if not exists(directory):
        makedirs(directory)
        return True
    return False

def check_file_exists(filename):
    if not isfile(filename):
        print 'Missing: ' + filename
        return False
    return True

if __name__ == "__main__":

    # Get input file names:
    project_directory = args.project_directory
    print 'Loading files from directory:  ' + project_directory
    input_files = get_gaussian_input_files(project_directory, args.input_files)

    # Reset variable will determine whether to reset existing data in a directory to start anew (e.g., if
    # previous run data was generated from an incorrect gaussian.inp file
    start_from_scratch = args.reset

    # Creates an output directory if not already present
    output_directory = default_output_directory
    if create_directory(output_directory):
        start_from_scratch = True

    # Initialize gaussian.inp file parameters
    input_data = {}
    for input_file_name in input_files:
        gaussian_input_params = {'link0_section': '', 'route_section': '',
                    'title': '', 'natoms': 0, 'charge': None, 'multiplicity': None, 'atom_coordinates' : ''}
        is_new_structure = False
        structure = input_file_name.split('.')[0]

        # Creates output directories for each structure if not already present
        structure_output_directory = join(output_directory, structure)
        if create_directory(structure_output_directory):
            is_new_structure = True

        # Checks that an existing structure directory has the correct files to continue running ART
        # The option will be given to skip directory and submit files for the remaining gaussian.inp
        # or to clear the directory and restart
        if not is_new_structure and not start_from_scratch:
            file_missing = False
            if not check_file_exists(join(structure_output_directory, '.'+gaussian_execution_script)):
                file_missing = True
            if not check_file_exists(join(structure_output_directory, gaussian_art_filename)):
                file_missing = True
            if not check_file_exists(join(structure_output_directory, refconfig_filename)):
                file_missing = True
            if not check_file_exists(join(structure_output_directory, 'filecounter')):
                file_missing = True
            if file_missing:
                print 'A critical file is missing from: ' + structure_output_directory + '\n'
                question = 'Reset this directory erasing it\'s data (y) or skip structure (n)'
                erase = query_yes_no(question)
                if erase:
                    rmtree(structure_output_directory)
                    if create_directory(structure_output_directory):
                        is_new_structure = True

        # This is only called when the user wants to restart the structure or it was not already set
        if start_from_scratch or is_new_structure:
            # Loads input data for each Gaussian input file
            input_data[structure] = load_input(join(project_directory, input_file_name), gaussian_input_params)

            create_ref_config(gaussian_input_params['atom_coordinates'],structure_output_directory)
            set_env_config(gaussian_input_params['natoms'])
            create_gaussian_file_header(gaussian_input_params, structure_output_directory)
            create_file_counter(file_counter_start, structure_output_directory)

            # TODO Temporary method to simply copy scripts to appropriate structure directories
            copy(join(script_directory, gaussian_art_filename), structure_output_directory)

        submission_type = args.submission_type
        if submission_type == 'GREX':
            set_submission_script(grex_submission_script, structure_output_directory, gaussian_input_params['link0_section'], args.optimize, time = args.time, memory = args.memory, processor_info = args.processor_info)
            print 'Running submission file for: ' + structure
            wd = getcwd()
            chdir(join(wd, structure_output_directory))
            # call(['qsub', '-N','gau_art_'+ structure, grex_submission_script], shell=False)
            chdir(wd)
        else:
            print 'Only GREX submission is currently available'
        #TODO handle other submission types