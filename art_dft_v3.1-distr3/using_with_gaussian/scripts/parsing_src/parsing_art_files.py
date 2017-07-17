
from os.path import join

def set_env_config(script_directory, gaussian_art_filename, refconfig_filename, output_directory, natoms):
    """
    Sets critical parameters in gaussian_art.sh which contains the general configuration for the ART environment

    :param natoms:
    :return:
    """
    max_number_of_events = 11
    central_atom = 30
    write_xyz = False
    radius_initial_deformation = None
    art_executable_location = '../../../source/artdft'


    with open(join(script_directory, gaussian_art_filename)) as input:
        updated_text = ''
        for line in input:
            setting_central_atom = False
            original_line = line
            updated_line = ''

            #Remove comments from the strings
            line = line.split("#")
            comment_set = False
            if len(line) > 1:
                comment_set = True

            # #Checks that GAU is set for energy calculation and sets it if it is not
            setenv, new_value = 'ENERGY_CALC', 'GAU'
            if _check_setenv(setenv, line):
                updated_line = _set_setenv(setenv, new_value)

            setenv, new_value = 'NATOMS', str(natoms)
            if _check_setenv(setenv, line):
                updated_line = _set_setenv(setenv, new_value)

            # updated_line = _check_and_set_setenv('setenv NATOMS', line, str(natoms))
            setenv, new_value = 'REFCONFIG', refconfig_filename
            if _check_setenv(setenv, line):
                updated_line = _set_setenv(setenv, new_value)

            #Optional set or unset Write_xyz
            if write_xyz is not None:
                setenv ='Write_xyz'
                if write_xyz is True:
                    new_value = '.true.'
                elif write_xyz is False:
                    new_value = '.false.'
                if _check_setenv(setenv, line):
                    updated_line = _set_setenv(setenv, new_value)

            # Optional sets maximum events that will run in a job
            setenv, new_value = 'Max_Number_Events', str(max_number_of_events)
            if _check_setenv(setenv, line):
                updated_line = _set_setenv(setenv, new_value)

            # Optional sets maximum events that will run in a job
            if central_atom is not None:
                setenv, new_value = 'Central_Atom', str(central_atom)
                if _check_setenv(setenv, line):
                    setting_central_atom = True
                    updated_line = _set_setenv(setenv, new_value)

                # Optional sets radius_initial_deformation if that parameter was added
                if radius_initial_deformation is not None:
                    setenv, new_value = 'Radius_Initial_Deformation', str(radius_initial_deformation)
                    if _check_setenv(setenv, line):
                        updated_line = _set_setenv(setenv, new_value)

                # Sets to local event to use central_atom setting
                setenv, new_value = 'Type_of_Events', 'local'
                if _check_setenv(setenv, line):
                    updated_line = _set_setenv(setenv, new_value)

            #set ART executable location from script
            var_name, new_value = 'art_location', ''
            if _check_bash_variable(var_name, line):
                _set_bash_variable(var_name, new_value)

            updated_text = _update_line(line, original_line, comment_set, updated_line, updated_text)

    # overwrites the configuration file with appropriate values from the gaussian input file
    config = open(join(output_directory, gaussian_art_filename), 'w+')
    config.write(updated_text)
    config.close()


def _update_line(line, original_line, comment_set, updated_line, updated_text):
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

    return updated_text


def _check_bash_variable(variable_name, line):

    substr = variable_name + '='
    if line[0].find(substr) != -1:
        return True
    return False


def _set_bash_variable(variable_name, new_value):
    return variable_name + '=' + new_value

def _check_setenv(setenv, line):

    substr = 'setenv ' + str(setenv)
    #Check for uncommented or commented versions of the string
    for section in line:
        if section.find(substr) != -1:
            return True

    return False

def _set_setenv(setenv, new_value):
    space_before_comment_section = '               '
    return 'setenv ' + setenv + '  ' + new_value + space_before_comment_section


def create_file_counter(file_counter_start, output_directory):
    """
    Creates a file counter that is used by the ART application to indicate the min/saddle file event
    :param file_counter_start:
    :param output_directory:
    """
    file_counter = open(join(output_directory, 'filecounter'), 'w+')
    file_counter.write('Counter:      ' + str(file_counter_start))
    file_counter.close()

def create_ref_config(refconfig_filename, atom_coordinates, output_directory):
    """
    Sets up refconfig.dat file which will be used as a starting reference by ART for atom coordinates before it is
    overwritten and updated

    :param refconfig:
    :param gaussian_input_params:
    :return:
    """

    config = open(join(output_directory, refconfig_filename),
                  'w+')  # Overwrites or creates a new file if it doesn't exist
    config.write('run_id:         1000\n')
    config.write('total_energy:   0\n')  # Placeholder, as this will be optimized by ART to the correct value
    # TODO see about removing these as they are not necessary for gaussian
    config.write('S   100.000000000000        100.000000000000        100.000000000000\n')
    config.write(atom_coordinates)
    config.close()