
# This is an adapter to allow the user to supply only a Gaussian formatted input file and have this run through ART
# For reference: http://gaussian.com/input/

import re

input_file = 'testgaussian.inp'
gaussian_execution_script = 'execute_gaussian.sh'
refconfig_file = 'refconfig.dat'            #
gaustart_file = 'gaustart.sh'               #Shell script containing the configuration parameters for the ART application

def create_ref_config(gaussian_input_params):
    """
    Sets up refconfig.dat file which will be used as a starting reference by ART for atom coordinates before it is
    overwritten and updated

    :param refconfig:
    :param gaussian_input_params:
    :return:
    """
    global refconfig_file
    global gaustart_file

    config = open(refconfig_file, 'w+')          #Overwrites of creates a new file if it doesn't exist
    config.write('run_id:         1000\n')
    config.write('total_energy:   0\n')     #Placeholder, as this will be optimized by ART to the correct value
    config.write('S   100.000000000000        100.000000000000        100.000000000000\n')
    config.write(gaussian_input_params['atom_coordinates'])
    config.close()

def set_env_config(gaussian_input_params):
    """
    Sets gaustart.sh which contains the general configuration for the ART environment. These will eventually be written
    to the file submitted to Gaussian in the execute_gaussian.sh script

    :param gaustart_file:
    :return:
    """
    global refconfig_file
    global gaustart_file

    with open(gaustart_file) as input:

        updated_text = ''
        for line in input:
            original_line = line
            updated_line = ''
            comment_set = False

            #Remove comments from the strings
            line = line.split("#")
            try:
                if line[1].strip() != '':
                    comment_set = True
            except:
                comment_set = False

            #Checks that GAU is set for energy calculation and sets it if it is not
            if line[0].find('setenv ENERGY_CALC') != -1:
                updated_line = 'setenv ENERGY_CALC GAU      '

            #Sets the number of atoms
            elif line[0].find('setenv NATOMS') != -1:
                updated_line = 'setenv NATOMS     '+ str(gaussian_input_params['natoms']) + '         '

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
    test = open(gaustart_file, 'w+')
    test.write(updated_text)


def create_gaussian_file_header(gaussian_input_params):
    """
    This serves to insert a gaussian header file into the bash script that will be writing to
    the art2gaussian.inp file and calling gaussian

    :param gaussian_input_params:
    :return:
    """
    global gaussian_execution_script
    params = gaussian_input_params

    # creates a regular expression pattern that will isolate the header insertion point
    insertion_point = re.compile('#gaussian-header-begin.*?#gaussian-header-end', re.DOTALL)

    # open file
    f = open(gaussian_execution_script, 'r')
    data = f.read()
    f.close()

    output = open(gaussian_execution_script, 'w+')

    # Builds a header for the gaussian.inp file, <OPTION> Flag will be replaced by opt or force
    header = '#gaussian-header-begin (DO NOT REMOVE) \n' \
             + 'header=\'' \
             + (params['link0_section']) \
             + ('#' + params['route_section'] + ' <OPTION>' + '\n') + ('\n') \
             + (params['title'] + '\n') + ('\n') \
             + (str(params['charge']) + ' ' + str(params['multiplicity']) + '\n') \
             + '\''\
             + '\n#gaussian-header-end'
    output.close()

    # Writes the header as a string variable in the gaussian execution script
    script_with_header = insertion_point.sub(header, data)
    f = open(gaussian_execution_script, 'w')
    f.write(script_with_header)
    f.close()


def load_input(gaussian_input_params):
    """
    :param gaussian_input_params
    :return:
    """
    params = gaussian_input_params
    section_number = 1

    with open(input_file) as input:

        for line in input:
                # This will load the relevant information from a gaussian input file expecting:
                # 1. A number of % lines (or maybe none)
                # 2. One # line describing the job
                # 3. Blank line, title line, blank line
                # 4. Charge + Spin multiplicity line (two integers)
                # 5. Atom coordinates

                #Check for an empty line to increment section
                if line.strip() == '':
                    section_number += 1
                    continue

                #Link0 and route section
                if section_number == 1:
                    if line.find('%') != -1:
                        params['link0_section'] = params['link0_section'] + line + " "
                    elif line.find('#') != -1:
                        #TODO handle multiple instances of # lines
                        #Removes hashtag for parsing later
                        line = line.split('#')
                        #Removes opt and force, which will be put back in by the ART code at different stages
                        #TODO handle cases like Opt=ModRedun
                        line[1] = line[1].replace('opt', '')
                        line[1] = line[1].replace('force', '')
                        params['route_section'] = line[1].strip()

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
                    params['atom_coordinates'] = params['atom_coordinates'] + line

                if section_number == 5:
                    break

        return params



if __name__ == "__main__":

    gaussian_input_params = {'link0_section': '', 'route_section': '',
                    'title': '', 'natoms': 0, 'charge': None, 'multiplicity': None, 'atom_coordinates' : ''}
    print load_input(gaussian_input_params)

    create_ref_config(gaussian_input_params)
    set_env_config(gaussian_input_params)
    create_gaussian_file_header(gaussian_input_params)
