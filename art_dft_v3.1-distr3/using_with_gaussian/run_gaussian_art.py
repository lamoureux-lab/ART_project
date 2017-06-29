
# This is an adapter to allow the user to supply only a Gaussian formatted input file and have this run through ART
# For reference: http://gaussian.com/input/

import re
import os.path
import argparse

script_folder = 'scripts'

input_file = 'sample.inp'
gaussian_execution_script = 'execute_gaussian.sh'
refconfig_file = 'refconfig.dat'            #
gaussian_art_file = 'gaussian_art.sh'               #Shell script containing the configuration parameters for the ART application
periodic_table_data = os.path.join(script_folder, 'periodic_table_data.csv')



#Handles file parameter passing
parser = argparse.ArgumentParser(description='Gaussian inputfile preperation for ART')
parser.add_argument('-f', '--input_file',help='The input file location')
args = parser.parse_args()



def create_ref_config(atom_coordinates):
    """
    Sets up refconfig.dat file which will be used as a starting reference by ART for atom coordinates before it is
    overwritten and updated

    :param refconfig:
    :param gaussian_input_params:
    :return:
    """
    global refconfig_file
    global gaussian_art_file

    config = open(refconfig_file, 'w+')     #Overwrites of creates a new file if it doesn't exist
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
    global refconfig_file
    global gaussian_art_file

    with open(gaussian_art_file) as input:

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
                updated_line = 'setenv REFCONFIG        ' + refconfig_file+ '             '

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
    test = open(gaussian_art_file, 'w+')
    test.write(updated_text)


def create_gaussian_file_header(gaussian_input_params):
    """
    Inserts a gaussian header file into the bash script that will be writing to
    the art2gaussian.inp file and calling gaussian

    :param gaussian_input_params:
    :return:
    """
    global gaussian_execution_script
    params = gaussian_input_params

    # creates a regular expression pattern that will isolate the header insertion point
    insertion_point = re.compile('#gaussian-header-begin.*?#gaussian-header-end \(DO NOT REMOVE\)', re.DOTALL)

    # open file
    f = open(gaussian_execution_script, 'r')
    initial_script = f.read()
    f.close()

    output = open(gaussian_execution_script, 'w+')

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
    f = open(gaussian_execution_script, 'w')
    f.write(script_with_header)
    f.close()

def get_coordinate_line_number(str):
    return len(str.split('\n'))

def load_input(gaussian_input_params):
    """
    Loads the values gaussian parameters into a dictionary object, while removing certain route parameters
    :param gaussian_input_params
    :return:
    """
    params = gaussian_input_params
    section_number = 1

    global input_file
    if (args.input_file):
        input_file = args.input_file
        print('Loading gaussian input file: ' + args.input_file)
    else:
        print('Input file location not set, using default gaussian input file: ' + input_file)

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

    atomic_symbol_number_map = {}
    with open(periodic_table_data) as input:
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

if __name__ == "__main__":
    load_periodic_table()

    gaussian_input_params = {'link0_section': '', 'route_section': '',
                    'title': '', 'natoms': 0, 'charge': None, 'multiplicity': None, 'atom_coordinates' : ''}
    gaussian_input_params = load_input(gaussian_input_params)

    create_ref_config(gaussian_input_params['atom_coordinates'])
    set_env_config(gaussian_input_params['natoms'])
    create_gaussian_file_header(gaussian_input_params)
