
# This is an adapter to allow the user to supply only a Gaussian formatted input file and have this run through ART

input_file = 'testgaussian.inp'
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

def set_env_config(gaussian_input_params):
    """

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

            #Remove comments from the strings
            line = line.split("#")
            #Checks that GAU is set for energy calculation and sets it if it is not
            if line[0].find('setenv ENERGY_CALC') != -1:
                updated_line = 'setenv ENERGY_CALC GAU      #' + line[1]

            #Sets the number of atoms
            elif line[0].find('setenv NATOMS') != -1:
                updated_line = 'setenv NATOMS     '+ str(gaussian_input_params['natoms']) + '         #' + line[1]

            #Sets reference configuration file that is updated throughout
            elif line[0].find('setenv REFCONFIG') != -1:
                updated_line = 'setenv REFCONFIG        ' + refconfig_file+ '             #' + line[1]
            elif line[0].find('setenv GAU_mem') != -1:
                updated_line = 'setenv GAU_mem        ' + gaussian_input_params['mem'] + '                 #' + line[1]
            elif line[0].find('setenv GAU_nproc') != -1:
                updated_line = 'setenv GAU_nproc      ' + gaussian_input_params['nproc'] + '                   #' + line[1]
            elif line[0].find('setenv GAU_desc') != -1:
                updated_line = 'setenv GAU_desc       ' + gaussian_input_params['description'] + '       #' + line[1]
            elif line[0].find('setenv GAU_charge') != -1:
                updated_line = 'setenv GAU_charge     ' + str(gaussian_input_params['charge']) + '                           #' + line[1]
            elif line[0].find('setenv GAU_multip') != -1:
                updated_line = 'setenv GAU_multip     ' + str(gaussian_input_params['multiplicity']) + '                           #' + line[
                    1]

            # updates the configuration file string
            if updated_line:
                updated_text = updated_text + updated_line
            else:
                updated_text = updated_text + original_line

    # overwrites the configuration file with appropriate values from the gaussian input file
    test = open(gaustart_file, 'w+')
    test.write(updated_text)



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

                #Check whether the line is in the first section
                #TODO find out whether to expect more than these two % lines for other_details section
                if section_number == 1:
                    if line.find('%mem') != -1:
                        params['mem'] = line.strip()
                    elif line.find('%nproc') != -1:
                        params['nproc'] = line.strip()
                    elif line.find('%') != -1:
                        params['other_details'] = params['other_details'] + line + " "
                    elif line.find('#') != -1:
                        params['description'] = line.strip()

                #Title section
                elif section_number == 2:
                    params['title'] = line.strip()

                #Charge + Spin multiplicity section
                if section_number == 3:
                    string_integers = line.split()
                    params['charge'] = int(string_integers[0])
                    params['multiplicity'] = int(string_integers[1])
                    section_number += 1
                    continue

                #Atom coordinates
                if section_number == 4:
                    params['natoms'] = params['natoms'] + 1
                    params['atom_coordinates'] = params['atom_coordinates'] + line

                if section_number == 5:
                    break

        return params



if __name__ == "__main__":


    gaussian_input_params = {'mem': '', 'nproc': '', 'other_details': '', 'description': '',
                    'title': '', 'natoms': 0, 'charge': None, 'multiplicity': None, 'atom_coordinates' : ''}
    print load_input(gaussian_input_params)

    create_ref_config(gaussian_input_params)
    set_env_config(gaussian_input_params)
